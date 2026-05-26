"""Flatten the raw ProteinBase export into the canonical designs table.

Reads:
    data/raw/proteinbase/proteinbase_collection_rbx1-binder-competition-results.csv

Writes:
    data/designs.csv
    data/designs.parquet
    data/designs.fasta

The raw export ships six columns (`id`, `name`, `sequence`, `author`,
`designMethod`, `evaluations`). The last is a stringified JSON list of
metric records. This script unpacks every metric into a flat column with
a stable name and writes the result to `data/designs.csv`.

The unpacking rules — which JSON `metric` becomes which column — are
documented in `docs/DATA.md`. Add a new mapping in `_METRIC_TO_COLUMN`
when ProteinBase ships a new metric and re-run this script.

Run via::

    mise run build
    # or:
    uv run python scripts/data/build_designs.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.utils.load_data import repo_root

# ---------------------------------------------------------------------------
# Constants — denominators, controls, mappings
# ---------------------------------------------------------------------------

# Author handle used in the public ProteinBase release for spike-in controls. These rows are tagged
# `is_control=True` and dropped by `load_designs(drop_controls=True)`.
# The ProteinBase release tags platform-authored spike-in controls with this
# author handle. Kept as data — not a brand reference, just the verbatim
# string the upstream export uses to mark control rows.
_CONTROL_AUTHOR = "adaptyv-bio"

# Scalar metrics — one record per design (or aggregate-ready). Anything not in
# this map is dropped on build; add a row to bring a new metric through. Keep
# names stable — analyses depend on them.
_SCALAR_METRIC_TO_COLUMN: dict[str, str] = {
    # Sequence-derived (per design)
    "foldstring":               "pb_foldstring",
    "molecular_weight":         "pb_molecular_weight",
    "isoelectric_point":        "pb_isoelectric_point",
    "classification":           "pb_classification",
    "design_class":             "pb_design_class",
    "ted_confidence":           "pb_ted_confidence",
    "novelty":                  "pb_novelty",
    # In-silico folding (per design)
    "esmfold_plddt":            "pb_esmfold_plddt",
    "proteinmpnn_score":            "pb_proteinmpnn_score",
    "proteinmpnn_seq_recovery":     "pb_proteinmpnn_seq_recovery",
    "redesigned_proteinmpnn_score": "pb_redesigned_proteinmpnn_score",
}

# Replicate-level metrics — multiple records per design. We aggregate per design.
_REPLICATE_BOOL_METRICS = {"expressed", "binding"}
_REPLICATE_STR_METRICS = {"binding_strength"}
_REPLICATE_FLOAT_METRICS = {"kd", "koff", "kon"}

# Artifact-URL metrics — value is `{"url": "..."}`. May appear multiple times
# per design (one per replicate sensorgram); we collect all URLs.
_URL_METRIC_TO_COLUMN: dict[str, str] = {
    "esmfold_structure_prediction": "pb_esmfold_cif_url",
    "esmfold_stylized_image":       "pb_stylized_png_url",
}
_MULTI_URL_METRIC_TO_COLUMN: dict[str, str] = {
    "spr_kinetic_curves": "pb_spr_curves_urls",
    "bli_kinetic_curves": "pb_bli_curves_urls",
}

# Nested-dict metrics handled by hand below.
# - `seqidentity`: dict with a `value` key (number)
# - `domainmatch`: dict with a `metrics` array (TM-score, evalue, seqIdentity vs AFDB50)
_BINDING_STRENGTH_RANK = {
    "Strong": 4,
    "Medium": 3,
    "Weak": 2,
    "Potential binder": 1,
    "Non-binder": 0,
    "No expression": -1,
    "Unknown": -2,
    "None": -3,
    None: -4,
}

# Method-family normalisation. Free-text designMethod → bucket. Add new
# entries here when a fresh method string shows up. **Order matters** —
# more specific patterns must come before more general ones (e.g.
# `bindcraft-2` before `bindcraft`, `boltz.*proteinmpnn` before plain
# `proteinmpnn`).
_METHOD_FAMILY_PATTERNS: list[tuple[str, str]] = [
    # ---- Combos (must come first; more specific than the constituent tools)
    (r"rfantibody.*ppiflow",    "RFantibody + PPIFlow"),
    (r"originflow.*prosllam",   "OriginFlow + ProSLLaM"),
    (r"originflow.*prodesign",  "OriginFlow + ProDESIGN-LE"),
    (r"boltz.*proteinmpnn",     "Boltz + ProteinMPNN"),
    (r"pepmind.*alphafold",     "PepMind + AF3"),
    (r"pepmlm.*peptiverse",     "PepMLM + Peptiverse"),
    (r"bagel.*solumpnn",        "Bagel + SoluMPNN"),
    (r"sae-steered.*esm",       "SAE-steered ESM2 + Boltz-2"),
    (r"adflip",                 "AF2 hallucination + ADFlip"),
    # ---- Standalone tools — order: most-specific-first inside each family
    (r"ligandforge",            "LigandForge"),
    (r"originflow",             "OriginFlow"),
    (r"prosllam",               "ProSLLaM"),
    (r"prodesign[-_\s]?le|prodesign",  "ProDESIGN-LE"),
    (r"rfantibody",             "RFantibody"),
    (r"rfdiffusion",            "RFdiffusion"),
    (r"rfpeptides",             "RFpeptides"),
    (r"bindcraft[-_\s]*2",      "BindCraft 2"),
    (r"bindcraft",              "BindCraft"),
    (r"boltzgen",               "BoltzGen"),
    (r"\bboltz",                "Boltz"),
    (r"lfm2",                   "LFM2"),
    (r"moppit",                 "MoPPIt"),
    (r"foldcraft",              "FoldCraft"),
    (r"hallucinat",             "Hallucination"),
    (r"proteinmpnn",            "ProteinMPNN"),
    (r"alphafold|af2|af3",      "AF-based"),
    (r"mosaic",                 "Mosaic"),
    (r"orbit",                  "ORBIT"),
    (r"ppiflow",                "PPIFlow"),
    (r"rosetta",                "Rosetta"),
    (r"proos",                  "ProOS"),
    (r"pxdesign",               "pXdesign"),
    (r"^giraf",                 "GIRAF"),
    (r"^arena",                 "Arena"),
    (r"jointdiff",              "JointDiff"),
    (r"protein-hunter|caliby",  "Protein Hunter / Caliby"),
    (r"llm-guided.*pipeline",   "LLM-guided pipeline"),
    (r"tea-leaves",             "Tea Leaves"),
    (r"iggm",                   "IgGM"),
    (r"binding-site-detection", "Binding-site detection workflow"),
    (r"flow[-_\s]*match",       "Flow-matching"),
    (r"swapped.*conformation",  "Swapped-conformation (homolog)"),
    (r"esm",                    "ESM-based"),
    (r"diffus",                 "Diffusion (other)"),
]


def _method_family_from_text(text: str | None) -> str | None:
    """Fallback: derive method family from any free-text blob.

    Used against `submission_core_models` first, then `submission_method_summary`.
    Returns None if no pattern matches.
    """
    if not isinstance(text, str) or not text.strip():
        return None
    t = text.lower()
    for pat, label in _METHOD_FAMILY_PATTERNS:
        if re.search(pat, t):
            return label
    return None


# ---------------------------------------------------------------------------
# Evaluation-record unpacking
# ---------------------------------------------------------------------------


def _unwrap_scalar(raw: Any) -> Any:
    """Pull the scalar value out of a ProteinBase metric record."""
    if isinstance(raw, dict):
        for k in ("value", "score", "distance"):
            if k in raw:
                inner = raw[k]
                if isinstance(inner, dict):
                    for kk in ("value", "score"):
                        if kk in inner:
                            return inner[kk]
                    return None
                return inner
        return None
    return raw


def _unwrap_domainmatch(value: Any) -> dict[str, Any]:
    """Pull TM-score, seqidentity, and evalue out of a `domainmatch` record."""
    if not isinstance(value, dict):
        return {}
    out: dict[str, Any] = {}
    for m in value.get("metrics", []) or []:
        if not isinstance(m, dict):
            continue
        slug = m.get("slug") or ""
        v = m.get("value")
        if isinstance(v, dict):
            v = v.get("seqIdentity") if "seqIdentity" in v else v.get("value")
        if slug == "seqidentity_afdb50":
            out["pb_seqidentity_afdb50"] = v
        elif slug == "evalue_afdb50":
            out["pb_evalue_afdb50"] = v
        elif slug in ("tm_score_afdb50", "tmscore_afdb50"):
            out["pb_tm_score_afdb50"] = v
    if "databaseId" in value:
        out["pb_afdb50_top_id"] = value.get("databaseId")
    return out


def _parse_evaluations(blob: str | None) -> dict[str, Any]:
    """Return a flat `{column: value}` dict for one row's `evaluations` JSON.

    Multi-record metrics (`expressed`, `binding`, `binding_strength`, `kd`,
    `koff`, `kon`) are aggregated per design: bools by `any()`, strings by the
    hierarchy `Strong > Medium > Weak > … > None`, floats by mean / min / max.
    """
    if not blob or not isinstance(blob, str):
        return {}
    try:
        records = json.loads(blob)
    except json.JSONDecodeError:
        return {}

    out: dict[str, Any] = {}
    bool_acc: dict[str, list[bool]] = {k: [] for k in _REPLICATE_BOOL_METRICS}
    strength_acc: list[str] = []
    float_acc: dict[str, list[float]] = {k: [] for k in _REPLICATE_FLOAT_METRICS}
    url_acc: dict[str, list[str]] = {k: [] for k in _MULTI_URL_METRIC_TO_COLUMN}
    seen_metrics: set[str] = set()

    for rec in records:
        if not isinstance(rec, dict):
            continue
        metric = rec.get("metric")
        if metric is None:
            continue
        seen_metrics.add(metric)
        value = rec.get("value")

        if metric in _SCALAR_METRIC_TO_COLUMN:
            col = _SCALAR_METRIC_TO_COLUMN[metric]
            scalar = _unwrap_scalar(value)
            if scalar is not None and col not in out:
                out[col] = scalar
        elif metric in _URL_METRIC_TO_COLUMN:
            col = _URL_METRIC_TO_COLUMN[metric]
            if isinstance(value, dict) and "url" in value and col not in out:
                out[col] = value["url"]
        elif metric in _MULTI_URL_METRIC_TO_COLUMN:
            if isinstance(value, dict) and "url" in value:
                url_acc[metric].append(value["url"])
        elif metric == "seqidentity":
            if isinstance(value, dict) and "value" in value:
                out["pb_seqidentity"] = value["value"]
        elif metric == "domainmatch":
            out.update(_unwrap_domainmatch(value))
        elif metric in _REPLICATE_BOOL_METRICS:
            if isinstance(value, bool):
                bool_acc[metric].append(value)
        elif metric in _REPLICATE_STR_METRICS:
            if isinstance(value, str):
                strength_acc.append(value)
        elif metric in _REPLICATE_FLOAT_METRICS:
            try:
                float_acc[metric].append(float(value))
            except (TypeError, ValueError):
                pass

    # Aggregate replicates
    if bool_acc["expressed"]:
        out["any_expressed"] = bool(any(bool_acc["expressed"]))
        out["n_replicates"] = len(bool_acc["expressed"])
        out["n_replicates_expressed"] = int(sum(bool_acc["expressed"]))
    if bool_acc["binding"]:
        out["any_binding"] = bool(any(bool_acc["binding"]))
        out["n_replicates_binding"] = int(sum(bool_acc["binding"]))

    if strength_acc:
        best = max(strength_acc, key=lambda s: _BINDING_STRENGTH_RANK.get(s, -5))
        out["binding_strength"] = best

    for metric, urls in url_acc.items():
        if urls:
            col = _MULTI_URL_METRIC_TO_COLUMN[metric]
            out[col] = "|".join(urls)
            out[col.replace("_urls", "_count")] = len(urls)

    for k, vals in float_acc.items():
        if not vals:
            continue
        if k == "kd":
            out["kd_M_mean"] = float(sum(vals) / len(vals))
            out["kd_M_min"] = float(min(vals))
            out["kd_M_max"] = float(max(vals))
            out["n_kd_records"] = len(vals)
        elif k == "koff":
            out["koff_mean"] = float(sum(vals) / len(vals))
        elif k == "kon":
            out["kon_mean"] = float(sum(vals) / len(vals))

    # Provenance: which upstream metric tiers populated for this design.
    out["pb_n_metric_kinds"] = len(seen_metrics)
    out["has_wetlab"] = "expressed" in seen_metrics
    out["has_predictions"] = "esmfold_plddt" in seen_metrics
    out["has_homology"] = "domainmatch" in seen_metrics
    out["has_ted_classification"] = "classification" in seen_metrics

    return out


# ---------------------------------------------------------------------------
# Local typer-output merge — picks up reruns saved by run_proteintyper.py
# ---------------------------------------------------------------------------


def _metric_value(record: dict[str, Any]) -> Any:
    """Pull the scalar/value out of a TyperJobOutput metric record."""
    value = record.get("value")
    if isinstance(value, dict):
        if "value" in value:
            return value["value"]
        if "slug" in value:
            return value["slug"]
    return value


_S3_PUB_PREFIX = "s3://proteinbase-pub/"
_S3_PUB_HTTPS = "https://proteinbase-pub.t3.storage.dev/"


def _s3_to_https(url: Any) -> str | None:
    if not isinstance(url, str):
        return None
    if url.startswith(_S3_PUB_PREFIX):
        return _S3_PUB_HTTPS + url[len(_S3_PUB_PREFIX):]
    if url.startswith("http"):
        return url
    return None


# Only keep AFDB50 columns that exist in the bulk-CSV upstream too — we
# previously also wrote `pb_rmsd_afdb50` and `pb_aligned_length_afdb50`, but
# the bulk parser doesn't extract them, so they sat at 98% NaN. The raw
# values are still in `data/metrics/proteintyper/<pb_id>.json` for anyone
# who wants them on a per-rerun basis.
_DOMAINMATCH_KEY_TO_COLUMN: dict[str, str] = {
    "seqIdentity": "pb_seqidentity_afdb50",
    "evalue": "pb_evalue_afdb50",
    "tm_score": "pb_tm_score_afdb50",
}

# `pb_domain_esmfold_plddt` and `pb_cath_detail` likewise sat at 98% NaN
# because the bulk parser doesn't extract them. Dropped to keep designs.csv
# clean — read the typer JSON directly if you need per-domain values.
_TED_METRIC_TO_COLUMN: dict[str, str] = {
    "classification": "pb_classification",
    "ted_confidence": "pb_ted_confidence",
    "foldstring": "pb_foldstring",
}


def _domainmatch_value(record: dict[str, Any]) -> Any:
    """A domainmatch metric wraps its value as `{<metric_name>: scalar}`
    (e.g. `{"seqIdentity": 28.2}`), not the usual `{"value": ...}`."""
    name = (record.get("metric_type") or {}).get("name")
    value = record.get("value")
    if isinstance(value, dict):
        if name and name in value:
            return value[name]
        if "value" in value:
            return value["value"]
    return value


def _ted_value(record: dict[str, Any]) -> Any:
    """TED metric values: classification is `{"value": "Alpha Beta"}`,
    ted_confidence is `{"confidence": 2}`, foldstring is `{"slug": "...",
    "human_readable": "..."}`."""
    name = (record.get("metric_type") or {}).get("name")
    value = record.get("value")
    if isinstance(value, dict):
        if name == "ted_confidence" and "confidence" in value:
            return value["confidence"]
        if name == "foldstring":
            return value.get("human_readable") or value.get("slug")
        if "value" in value:
            return value["value"]
        if "slug" in value:
            return value["slug"]
    return value


def _parse_local_typer_output(path: Path) -> dict[str, Any]:
    """Return a `{column: value}` dict matching the bulk-CSV `_parse_evaluations`
    convention, but read from a local TyperJobOutput JSON saved by
    `run_proteintyper.py`. Structure:

        {
          sequence: {metrics: [scalar...]},                       # design-level
          proteindomains: [
            {
              metrics: [TED scalar metrics: classification/ted_confidence/...],
              domainmatches: [{metrics: [seqIdentity/evalue/tm_score/rmsd/...]}]
            }
          ],
          structures: [{metrics: [...], file: {url}, img: {url}}]  # ESMFold output
        }
    """
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    out: dict[str, Any] = {}

    def _walk(records):
        for rec in records or []:
            if not isinstance(rec, dict):
                continue
            name = (rec.get("metric_type") or {}).get("name")
            if not name:
                continue
            v = _metric_value(rec)
            if v is None:
                continue
            col = _SCALAR_METRIC_TO_COLUMN.get(name)
            if col and col not in out:
                out[col] = v

    _walk((data.get("sequence") or {}).get("metrics"))
    for s in data.get("structures") or []:
        if not isinstance(s, dict):
            continue
        _walk(s.get("metrics"))
        # Pull CIF + stylised-image URLs from the first ESMFold structure.
        # Without this the 4 reran designs leave pb_esmfold_cif_url /
        # pb_stylized_png_url null even though they have the assets on disk.
        method_name = (s.get("method") or {}).get("name", "").lower()
        if method_name == "esmfold" or "pb_esmfold_cif_url" not in out:
            cif_url = _s3_to_https((s.get("file") or {}).get("url"))
            png_url = _s3_to_https((s.get("img") or {}).get("url"))
            if cif_url and "pb_esmfold_cif_url" not in out:
                out["pb_esmfold_cif_url"] = cif_url
            if png_url and "pb_stylized_png_url" not in out:
                out["pb_stylized_png_url"] = png_url

    # TED + AFDB50 live one level deeper. Pull from the FIRST proteindomain /
    # FIRST domainmatch — that matches the bulk-CSV upstream convention.
    domains = data.get("proteindomains") or []
    if domains and isinstance(domains[0], dict):
        d0 = domains[0]
        # TED scalar metrics on the domain (classification, ted_confidence, ...).
        for rec in d0.get("metrics") or []:
            if not isinstance(rec, dict):
                continue
            name = (rec.get("metric_type") or {}).get("name")
            col = _TED_METRIC_TO_COLUMN.get(name)
            if not col or col in out:
                continue
            v = _ted_value(rec)
            if v is not None:
                out[col] = v
        # AFDB50 domainmatch metrics (seqIdentity, evalue, tm_score, ...).
        matches = d0.get("domainmatches") or []
        if matches and isinstance(matches[0], dict):
            for rec in matches[0].get("metrics") or []:
                if not isinstance(rec, dict):
                    continue
                name = (rec.get("metric_type") or {}).get("name")
                col = _DOMAINMATCH_KEY_TO_COLUMN.get(name)
                if not col or col in out:
                    continue
                v = _domainmatch_value(rec)
                if v is not None:
                    out[col] = v

    return out


def _designs_missing_pb_predictions(root: Path) -> dict[str, dict[str, Any]]:
    """Pre-scan `data/metrics/proteintyper/*.json` so the build can backfill
    rows whose ProteinBase predictions never ran upstream.
    """
    out: dict[str, dict[str, Any]] = {}
    typer_dir = root / "data" / "metrics" / "proteintyper"
    if not typer_dir.exists():
        return out
    for path in typer_dir.glob("*.json"):
        flat = _parse_local_typer_output(path)
        if flat:
            out[path.stem] = flat
    return out


# ---------------------------------------------------------------------------
# Derived columns
# ---------------------------------------------------------------------------


def _method_family(method: str | None) -> str:
    if not method or not isinstance(method, str):
        return "Not mentioned"
    m = method.lower()
    for pat, label in _METHOD_FAMILY_PATTERNS:
        if re.search(pat, m):
            return label
    return "Other"


def _is_binder(strength: Any) -> bool:
    if not isinstance(strength, str):
        return False
    return strength.strip().lower() in {"strong", "medium", "weak"}


def _is_strong(strength: Any) -> bool:
    if not isinstance(strength, str):
        return False
    return strength.strip().lower() == "strong"


def _is_expressed(any_expressed: Any) -> bool:
    if any_expressed is None:
        return False
    if isinstance(any_expressed, bool):
        return any_expressed
    s = str(any_expressed).strip().lower()
    return s in {"true", "yes", "high", "medium", "low", "expressed"}


def _pkd(kd_m: Any) -> float | None:
    try:
        v = float(kd_m)
        if v <= 0:
            return None
        import math
        return -math.log10(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _normalize_handle(s: str) -> str:
    """Lower-case + strip non-alphanumeric. `Tom Carroll` and `tom-carroll` collide."""
    if not isinstance(s, str):
        return ""
    return re.sub(r"[^a-z0-9]", "", s.lower())


# Hand-curated mapping for cases where the ProteinBase handle and the
# submission author_name don't normalise the same way. Each entry was
# verified by cross-checking the proteinbase `designMethod` token against
# the submission's `core_models` / `method_summary` to make sure the same
# tool family appears on both sides.
_HANDLE_OVERRIDES: dict[str, str] = {
    # ProteinBase handle → submissions-CSV author_name (verified matches)
    "pacesalab":        "Pacesa Lab",         # bindcraft-2 ↔ BindCraft
    "x.rustamov":       "Khondamir Rustamov", # foldcraft ↔ FoldCraft
    "rdgonzalez":       "Ricardo González",   # latin-1 name
    "d-barradas":       "Didier Barradas Bautista",
    "t-carroll":        "Thomas Carroll",     # boltzgen ↔ Boltz-2, BoltzGen
    "luciau":           "Lucia Urcelay",
    "getu-tadesse":     "GETU TADESSE FELLEK",
    "zhimeng-zhou":     "zhimeng",
    "ievapudz":         "Ieva Pudžiuvelytė",
    "wufandi":          "Fandi Wu",
    "tanay":            "Tanay Lohia",
    "guanjiaweitaskin": "管佳威",
    "hz3519":           "Haowen Zhao",        # giraf ↔ GIRAF pipeline
    "niki-iva":         "Nikita Ivanisenko",  # boltz2-proteinmpnn ↔ ProteinMPNN, Boltz-2
    "maxk":             "Maksim Kuznetsov",   # lfm2 ↔ LFM
    "pimi":             "Pi",                 # pxdesign ↔ PXDesign, ProteinMPNN

    # Confirmed via ProteinBase per-submission collection scrape
    # (data/raw/submissions/recovered_methods.csv):
    "nanogenomic":               "Andre Watson",       # → LigandForge
    "professionalmouthpipettor": "Christian Teague",   # → RFantibody + PPIFlow + ProteinMPNN
    "zhangpeioo":                "Pei Zhang",          # → OriginFlow + ProSLLaM + ProDESIGN-LE

    # The remaining handles below have no row in the submissions CSV at all
    # (no entry-form filed). Mark with empty string so the build doesn't
    # keep flagging them.
    "drtheone":                  "",  # 1 design, custom-binding-site workflow
    "falmassen":                 "",  # 7 designs, plain `boltzgen` (33 BoltzGen submissions — too ambiguous)
}


def _submissions_path(root: Path) -> Path:
    """Prefer the private (un-sanitized) submissions CSV if it exists.

    The public-release version drops PII columns (`author_name`, `submitted_at`,
    `link`, `votes`, `comments`) and uses pseudonymous `submission_id`. The
    private side-file `*.private.csv` is gitignored and keeps the raw
    fields needed for the team-handle join.
    """
    private = root / "data" / "raw" / "submissions" / "rbx1_designated_submissions.private.csv"
    if private.exists():
        return private
    return root / "data" / "raw" / "submissions" / "rbx1_designated_submissions.csv"


def _load_submissions(root: Path) -> pd.DataFrame:
    """Read the per-submission CSV at data/raw/submissions/, redact PII columns.

    Returns one row per submission with normalised columns:
        submission_author_name, submission_link, submission_modality,
        submission_core_models, submission_method_summary,
        submission_target_region, submission_targets_idr,
        submission_binder_length, submission_team_type,
        submission_new_method, submission_design_type,
        submission_validation_stack, submission_epitope_source,
        submission_uniref_check, submission_sabdab_check,
        submission_overall_homology, submission_votes,
        submission_n_proteins.
    Empty DataFrame if the file is missing.
    """
    path = _submissions_path(root)
    if not path.exists():
        return pd.DataFrame()
    raw = pd.read_csv(path)
    # Drop every PII / internal-review column at load. Public release ships
    # the sanitized version anyway, but be defensive in case build runs
    # against the private file.
    pii_drop = [c for c in ("author_email", "submitted_at", "link", "votes", "comments")
                if c in raw.columns]
    if pii_drop:
        raw = raw.drop(columns=pii_drop)

    # Scrub email addresses embedded in free-text columns. Some submitters
    # added "contact me at foo@bar.com" inside their method writeup; that's
    # PII and stays out of the public CSVs.
    email_re = re.compile(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+")
    for prose_col in ("description", "method_summary", "parsed_pdfs",
                       "validation_stack", "epitope_source", "core_models"):
        if prose_col in raw.columns:
            raw[prose_col] = (
                raw[prose_col]
                .astype("string")
                .apply(lambda s: email_re.sub("[email-redacted]", s) if isinstance(s, str) else s)
            )
    rename = {
        "author_name":          "submission_author_name",
        "n_proteins":           "submission_n_proteins",
        "description":          "submission_description",
        "total_submissions":    "submission_total_submissions",
        "modality":             "submission_modality",
        "core_models":          "submission_core_models",
        "method_summary":       "submission_method_summary",
        "target_region":        "submission_target_region",
        "targets_idr":          "submission_targets_idr",
        "binder_length":        "submission_binder_length",
        "team_type":            "submission_team_type",
        "new_method":           "submission_new_method",
        "design_type":          "submission_design_type",
        "num_generation_methods": "submission_num_generation_methods",
        "validation_stack":     "submission_validation_stack",
        "epitope_source":       "submission_epitope_source",
        "uniref_check":         "submission_uniref_check",
        "sabdab_check":         "submission_sabdab_check",
        "overall_homology":     "submission_overall_homology",
        "votes":                "submission_votes",
    }
    raw = raw.rename(columns=rename)
    raw["_handle_key"] = raw["submission_author_name"].map(_normalize_handle)
    return raw


def _build_team_to_submission_index(
    team_handles: list[str], submissions: pd.DataFrame
) -> dict[str, pd.Series]:
    """For each ProteinBase team handle, return the best matching submission row.

    A handle may have multiple submissions; we pick the one with the largest
    `submission_n_proteins`. Missing matches return an empty row.
    """
    if submissions.empty:
        return {h: pd.Series(dtype=object) for h in team_handles}

    index: dict[str, pd.Series] = {}
    sub_by_key = (
        submissions.sort_values("submission_n_proteins", ascending=False)
        .drop_duplicates("_handle_key", keep="first")
        .set_index("_handle_key")
    )
    for handle in team_handles:
        key = _normalize_handle(handle)
        override = _HANDLE_OVERRIDES.get(handle)
        if override:
            key = _normalize_handle(override)
        index[handle] = sub_by_key.loc[key] if key in sub_by_key.index else pd.Series(dtype=object)
    return index


def build(
    raw_csv: Path | None = None,
    out_dir: Path | None = None,
) -> pd.DataFrame:
    root = repo_root()
    if raw_csv is None:
        raw_csv = (
            root
            / "data"
            / "raw"
            / "proteinbase"
            / "proteinbase_collection_rbx1-binder-competition-results.csv"
        )
    if out_dir is None:
        out_dir = root / "data"

    raw = pd.read_csv(raw_csv)
    if "id" not in raw.columns:
        raise RuntimeError(
            f"Raw export at {raw_csv} is missing the `id` column — "
            "did the schema change?"
        )

    submissions = _load_submissions(root)
    team_handles = [(getattr(r, "author", None) or "") for r in raw.itertuples(index=False)]
    sub_by_handle = _build_team_to_submission_index(team_handles, submissions)

    # Build the author_name → submission_id pseudonym map. Sorted alpha so
    # the mapping is stable across re-builds, even if submission row order
    # changes upstream.
    global _SUB_ID_MAP, _TEAM_ID_MAP
    if not submissions.empty and "submission_author_name" in submissions.columns:
        names = sorted(
            {n for n in submissions["submission_author_name"].dropna() if str(n).strip()},
            key=lambda s: str(s).lower(),
        )
        width = max(3, len(str(len(names))))
        _SUB_ID_MAP = {n: f"sub_{i:0{width}d}" for i, n in enumerate(names, start=1)}
    else:
        _SUB_ID_MAP = {}

    # Team-handle pseudonym map — applied at write time. Sorted alpha so
    # rebuilds produce identical pseudonyms.
    handles = sorted({h for h in team_handles if h}, key=lambda s: s.lower())
    width = max(3, len(str(len(handles))))
    _TEAM_ID_MAP = {h: f"team_{i:0{width}d}" for i, h in enumerate(handles, start=1)}
    unmatched_handles = {
        h for h in team_handles
        if h and (h not in sub_by_handle or sub_by_handle[h].empty)
    }
    local_typer = _designs_missing_pb_predictions(root)

    # Manual recoveries from ProteinBase per-submission collection scrapes.
    # Highest-priority method-family source — overrides the regex inference
    # for the pb_ids listed. Schema: `pb_id, recovered_method, source`.
    recovered_methods: dict[str, str] = {}
    recovered_path = root / "data" / "raw" / "submissions" / "recovered_methods.csv"
    if recovered_path.exists():
        rec_df = pd.read_csv(recovered_path)
        recovered_methods = dict(zip(rec_df["pb_id"], rec_df["recovered_method"]))

    rows: list[dict[str, Any]] = []
    for i, r in enumerate(raw.itertuples(index=False), start=1):
        flat = _parse_evaluations(getattr(r, "evaluations", None))
        method = getattr(r, "designMethod", None)
        author = getattr(r, "author", None) or ""

        row: dict[str, Any] = {
            "design_id": i,
            "pb_id": getattr(r, "id", None),
            "name": getattr(r, "name", None),
            "team": author,
            "is_control": author.strip().lower() == _CONTROL_AUTHOR,
            "sequence": getattr(r, "sequence", "") or "",
            "design_method": method or "",
            "method_family": _method_family(method),
        }
        row["sequence_length"] = len(row["sequence"])
        row.update(flat)

        # If this design got a local typer rerun (data/metrics/proteintyper/<pb_id>.json),
        # backfill any pb_* column the upstream release left blank.
        local_pb_id = row.get("pb_id")
        if local_pb_id and local_pb_id in local_typer:
            for col, value in local_typer[local_pb_id].items():
                if col not in row or row.get(col) in (None, "", float("nan")):
                    row[col] = value
            row["pb_predictions_source"] = "local_rerun"
            row["has_predictions"] = True
            row["has_ted_classification"] = row.get("has_ted_classification", False) or bool(row.get("pb_classification"))
            row["has_homology"] = row.get("has_homology", False) or bool(row.get("pb_tm_score_afdb50"))
        else:
            row["pb_predictions_source"] = "proteinbase_release" if row.get("pb_esmfold_plddt") is not None else "missing"

        sub_row = sub_by_handle.get(author, pd.Series(dtype=object))
        # PII columns (`submission_author_name`, `submission_link`) are NOT
        # copied into designs.csv — they would leak into the public release.
        # `submission_id` is set below as a stable pseudonymous pointer.
        for col in (
            "submission_modality", "submission_core_models",
            "submission_method_summary", "submission_target_region",
            "submission_targets_idr", "submission_binder_length",
            "submission_team_type", "submission_new_method",
            "submission_design_type", "submission_num_generation_methods",
            "submission_validation_stack", "submission_epitope_source",
            "submission_uniref_check", "submission_sabdab_check",
            "submission_overall_homology", "submission_n_proteins",
            "submission_total_submissions",
        ):
            if sub_row is not None and col in sub_row.index:
                value = sub_row.get(col)
                if pd.isna(value):
                    value = None
                row[col] = value

        # Stable pseudonym for the submission this design was filed under.
        # Derived deterministically by sanitize_pii.py from the original
        # author_name; we look it up here so designs.csv stays joinable to
        # the public submissions CSV without exposing the name.
        if sub_row is not None and "submission_author_name" in sub_row.index:
            row["submission_id"] = _SUB_ID_MAP.get(
                (sub_row.get("submission_author_name") or ""), ""
            )

        # Method-family resolution priority:
        #   0. Manual recovery (data/raw/submissions/recovered_methods.csv)
        #      — authoritative; comes from per-submission ProteinBase
        #      collection scrapes that disambiguate cases the regex can't.
        #   1. proteinbase `designMethod` (already set above via regex)
        #   2. submission-form `core_models` blob
        #   3. submission-form `method_summary` prose
        pb_id = row.get("pb_id")
        if pb_id in recovered_methods:
            row["method_family"] = recovered_methods[pb_id]
            row["method_family_source"] = "recovered_from_collection"
        elif row["method_family"] == "Not mentioned":
            rescued = (
                _method_family_from_text(row.get("submission_core_models"))
                or _method_family_from_text(row.get("submission_method_summary"))
            )
            if rescued:
                row["method_family"] = rescued
                row["method_family_source"] = (
                    "submission_core_models"
                    if _method_family_from_text(row.get("submission_core_models"))
                    else "submission_method_summary"
                )
            else:
                row["method_family_source"] = "missing"
        else:
            row["method_family_source"] = "proteinbase_designMethod"

        strength = row.get("binding_strength")
        row["is_binder"] = _is_binder(strength)
        row["is_strong"] = _is_strong(strength)
        row["is_expressed"] = _is_expressed(row.get("any_expressed"))
        row["pkd_arith_mean"] = _pkd(row.get("kd_M_mean"))
        # Explicit flag for designs that hit the wet-lab pipeline at all.
        # Public release has 1 design (ivory-dove-granite) that lacks any
        # replicate record — never made it onto a chip.
        row["tested_in_wet_lab"] = bool(row.get("n_replicates")) and row["n_replicates"] > 0

        if isinstance(row.get("kd_M_mean"), (int, float)):
            row["kd_nM_mean"] = float(row["kd_M_mean"]) * 1e9
        if isinstance(row.get("kd_M_min"), (int, float)):
            row["kd_nM_min"] = float(row["kd_M_min"]) * 1e9
        if isinstance(row.get("kd_M_max"), (int, float)):
            row["kd_nM_max"] = float(row["kd_M_max"]) * 1e9

        rows.append(row)

    df = pd.DataFrame(rows)
    # Replace the raw ProteinBase handle with the deterministic pseudonym.
    # The original handle is preserved out-of-band in `data/team_pseudonyms.csv`
    # (gitignored) so the team can still trace their own rows back.
    df["team_id"] = df["team"].map(_TEAM_ID_MAP).fillna("team_unknown")
    df = df.drop(columns=["team"])

    # Persist the private team-handle → team_id mapping out-of-band.
    mapping = sorted(_TEAM_ID_MAP.items(), key=lambda kv: kv[1])
    (out_dir / "team_pseudonyms.csv").write_text(
        "team_handle,team_id\n" + "\n".join(f"{h},{tid}" for h, tid in mapping) + "\n"
    )

    # Compute the method-known count BEFORE dropping the source column so we
    # can still report it at the end. The audit file we used to write was
    # redundant with this single number — removed for public-release tidiness.
    if "method_family_source" in df.columns:
        n_method_known = int((df["method_family_source"] != "missing").sum())
    else:
        n_method_known = int(df["method_family"].notna().sum())

    # Drop columns the public CSV shouldn't carry:
    #   - design_method         (raw free-text, redundant with method_family)
    #   - method_family_source  (provenance debug — kept only in git history)
    #   - n_kd_records          (only 9 binders, redundant with kd_nM_mean.notna())
    #   - pb_bli_curves_urls    (legacy BLI fallback — 3 designs only;
    #                            primary screen was SPR, see pb_spr_curves_*)
    #   - pb_bli_curves_count   (same)
    drop_cols = [c for c in (
        "design_method",
        "method_family_source",
        "n_kd_records",
        "pb_bli_curves_urls",
        "pb_bli_curves_count",
    ) if c in df.columns]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    leading = [
        "design_id", "pb_id", "name", "team_id", "submission_id",
        "is_control", "method_family",
        "sequence", "sequence_length",
        "tested_in_wet_lab",
        "pb_n_metric_kinds", "has_wetlab", "has_predictions",
        "has_homology", "has_ted_classification",
        "pb_predictions_source",
        "binding_strength", "is_binder", "is_strong",
        "any_expressed", "is_expressed",
        "n_replicates", "n_replicates_expressed", "n_replicates_binding",
        "kd_M_mean", "kd_nM_mean", "pkd_arith_mean",
        "kd_M_min", "kd_M_max", "kd_nM_min", "kd_nM_max",
        "koff_mean", "kon_mean",
        # Submission-CSV metadata (one row per submission joined onto every design).
        # submission_author_name + submission_link are DROPPED — PII.
        "submission_modality",
        "submission_target_region", "submission_targets_idr",
        "submission_team_type", "submission_new_method",
        "submission_design_type", "submission_core_models",
        "submission_method_summary",
    ]
    ordered = [c for c in leading if c in df.columns] + [
        c for c in df.columns if c not in leading
    ]
    df = df[ordered]

    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "designs.csv", index=False)
    try:
        df.to_parquet(out_dir / "designs.parquet", index=False)
    except Exception as e:
        print(f"[build] WARN: could not write parquet ({e}). CSV is the source of truth.")

    fasta_path = out_dir / "designs.fasta"
    lines: list[str] = []
    for _, r in df.iterrows():
        ident = f"{r['design_id']}|{r['pb_id']}|{r['team_id']}"
        seq = r["sequence"]
        if not isinstance(seq, str) or not seq:
            continue
        lines.append(f">{ident}")
        for k in range(0, len(seq), 60):
            lines.append(seq[k : k + 60])
    fasta_path.write_text("\n".join(lines) + "\n")

    n = len(df)
    n_ctrl = int(df["is_control"].sum())
    n_expr = int(df["is_expressed"].sum())
    n_bind = int(df["is_binder"].sum())
    n_strong = int(df["is_strong"].sum())
    n_with_submission = int(df["submission_id"].notna().sum() & (df["submission_id"] != "").sum()) if "submission_id" in df.columns else 0
    # `n_method_known` was computed upstream from `method_family_source` before
    # that column was dropped. Falls back to `n` if the column was already gone
    # (older builds).
    n_method_missing = n - n_method_known
    print(
        f"[build] wrote {n} rows ({n_ctrl} controls). "
        f"expressed={n_expr}  binders={n_bind}  strong={n_strong}\n"
        f"[build] submission_meta={n_with_submission}/{n}  "
        f"method_known={n_method_known}/{n}  method_missing={n_method_missing}"
    )
    if unmatched_handles:
        print(
            f"[build] teams with no submission-form match ({len(unmatched_handles)}): "
            f"{sorted(unmatched_handles)}"
        )
        print(
            "[build]   (these are not bugs — they're authors absent from the "
            "submissions CSV; proteinbase designMethod is used where available)"
        )
    print(f"[build] {out_dir / 'designs.csv'}")
    print(f"[build] {out_dir / 'designs.fasta'}")
    return df


if __name__ == "__main__":
    build()
