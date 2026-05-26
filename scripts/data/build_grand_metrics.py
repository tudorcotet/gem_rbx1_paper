"""Pool per-model metric JSONs into a single wide ``grand_metrics.csv``.

Phase 3 of the rerun. After ``run_rescoring.py`` has filled
``data/metrics/{boltz2,protenix,chai,af2m}/`` and ``data/structures/{…}/``
with one file per design, this script joins them by ``pb_id`` and
writes ``data/grand_metrics.csv``.

Schema convention (matches the prefixes documented in
``docs/DATA.md`` § "Re-scoring"):

* ``pb_id``                — join key
* ``pb_*``                 — ProteinBase release columns lifted from
                             ``designs.csv`` (ESMFold + ProteinMPNN +
                             novelty / TED / classification / AFDB50).
                             URL-typed columns are dropped — they are
                             mirroring pointers, not metrics.
* ``tp_*``                 — ProteinTyper monomer rerun (same panel as
                             ``pb_*`` above, but re-computed locally
                             from ``data/metrics/proteintyper/``).
* ``b2_*``                 — Boltz-2 complex (chains A=target, B=binder)
* ``px_*``                 — Protenix complex
* ``chai_*``               — Chai-1 complex
* ``af2m_*``               — AlphaFold-2 Multimer complex

Plus four scoring layers fired on top of the per-model predicted
complexes:

* ``prodigy_{b2,chai,px,af2m}_*``     — PRODIGY ΔG / KD per source model
* ``destress_{esmfold,b2,chai,px,af2m}_*`` — DE-STRESS energy / packing
                                        panel per source model
                                        (esmfold = monomer, others =
                                        complex)
* ``esm_pll_*``           — ESM pseudo-log-likelihood (sequence-only,
                            no per-model nesting)
* ``netsolp_*``           — NetSolP solubility / usability (sequence-
                            only, no per-model nesting)

Plus a handful of derived consensus columns at the end:

* ``ipsae_pass_4folders``  — # of {b2, px, chai, af2m} with d0chn_max >= 0.4
* ``iptm_pass_4folders``   — # of {b2, px, chai, af2m} with iptm >= 0.7

The script is idempotent. Missing JSONs leave nulls; missing CIFs leave
``cif_*_exists`` = False. Missing directories (any of af2m, prodigy,
destress, esm_pll, netsolp) leave the corresponding column blocks as
nulls — re-run after the Modal jobs land.

Usage::

    mise run build:grand
    # equivalently:
    uv run python scripts/data/build_grand_metrics.py
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.utils.load_data import repo_root

# Models we expect under data/metrics/<model>/ and data/structures/<model>/.
# `proteintyper` produces its own folder but those columns are already on
# designs.csv via build_designs.py — we re-export them with `tp_` here so
# `grand_metrics.csv` is fully self-contained.
COMPLEX_MODELS: tuple[str, ...] = ("boltz2", "protenix", "chai", "af2m")
PREFIX: dict[str, str] = {
    "boltz2": "b2",
    "protenix": "px",
    "chai": "chai",
    "af2m": "af2m",
    "proteintyper": "tp",
}

# Fields we lift off the typer monomer panel onto the grand row. These are
# the canonical ProteinBase column names (matching `pb_*` on designs.csv).
TYPER_FIELDS: tuple[str, ...] = (
    "esmfold_plddt",
    "proteinmpnn_score",
    "proteinmpnn_seq_recovery",
    "redesigned_proteinmpnn_score",
    "molecular_weight",
    "isoelectric_point",
    "novelty",
    "seqidentity",
    "seqidentity_afdb50",
    "evalue_afdb50",
    "tm_score_afdb50",
    "ted_confidence",
    "design_class",
    "classification",
    "foldstring",
)

# Fields we lift off each complex predictor's JSON. Names mirror the per-
# model output of `compute_ipsae` plus the native confidence scalars.
COMPLEX_FIELDS_BASE: tuple[str, ...] = (
    "iptm",
    "ptm",
    "mean_plddt",
    "ipsae_d0res_min",
    "ipsae_d0res_max",
    "ipsae_d0chn_min",
    "ipsae_d0chn_max",
    "ipsae_d0dom_min",
    "ipsae_d0dom_max",
    "iptm_d0chn_min",
    "iptm_d0chn_max",
    "iptm_af_min",
    "iptm_af_max",
    "pdockq",
    "pdockq2",
    "lis",
    "n_interface",
)

# Predictor-specific extra fields.
EXTRA_FIELDS: dict[str, tuple[str, ...]] = {
    "boltz2": (),
    "protenix": ("ranking_score", "model_name"),
    "chai": ("aggregate_score",),
    "af2m": (),
}

# ---------------------------------------------------------------------------
# pb_* lifted from designs.csv
# ---------------------------------------------------------------------------

# URL / blob-pointer columns to drop when copying pb_* off designs.csv —
# they are mirroring pointers, not metrics.
PB_DROP_COLUMNS: frozenset[str] = frozenset(
    {
        "pb_esmfold_cif_url",
        "pb_stylized_png_url",
        "pb_spr_curves_url",
        "pb_spr_curves_urls",
        "pb_spr_curves_count",
    }
)

# ---------------------------------------------------------------------------
# Scoring-layer column blocks
# ---------------------------------------------------------------------------

# Per-source-model scoring layers. The on-disk JSON shape is
# ``{<source_model>: {<field>: <value>, …}, …}`` where ``<source_model>`` is
# one of the strings in ``SCORING_SOURCE_MODELS`` below.
#
# We flatten as ``<layer>_<source_sub_prefix>_<field>``. e.g. ProtENIX
# prodigy → ``prodigy_px_dG``; ESMFold monomer DE-STRESS →
# ``destress_esmfold_rosetta_total``.
#
# DE-STRESS runs on every complex AND on the ESMFold monomer (that's why
# its source list is longer than PRODIGY's).
SCORING_SOURCE_PREFIX: dict[str, str] = {
    "boltz2": "b2",
    "protenix": "px",
    "chai": "chai",
    "af2m": "af2m",
    "esmfold": "esmfold",
}

PRODIGY_SOURCES: tuple[str, ...] = ("boltz2", "chai", "protenix", "af2m")
DESTRESS_SOURCES: tuple[str, ...] = ("esmfold", "boltz2", "chai", "protenix", "af2m")

# We do not hard-code prodigy / destress sub-field names because they may
# evolve as the Modal apps stabilise. Instead we discover them from the
# first JSON we manage to load per directory, then keep a stable union
# across all designs. If no JSON exists yet, the column block stays empty
# and downstream rows just carry nulls for that layer.


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _flatten(
    pb_id: str, model: str, data: dict[str, Any] | None, fields: Iterable[str]
) -> dict[str, Any]:
    """Return ``{prefix_field: value}`` for one (pb_id, model) row."""
    prefix = PREFIX[model]
    out: dict[str, Any] = {}
    out[f"{prefix}_status"] = (data or {}).get("status")
    for f in fields:
        out[f"{prefix}_{f}"] = (data or {}).get(f)
    return out


def _collect_typer(pb_id: str) -> dict[str, Any]:
    """Extract the same monomer panel ProteinBase ships, from the local
    typer JSON if present. Falls back to nulls when typer JSON is absent
    (the build script for designs.csv has the same data — this re-export
    is just so grand_metrics is self-contained)."""
    typer_path = repo_root() / "data" / "structures" / "typer_outputs" / f"{pb_id}.json"
    data = _read_json(typer_path)
    flat: dict[str, Any] = {"tp_status": "ok" if data else None}
    for f in TYPER_FIELDS:
        flat[f"tp_{f}"] = (data or {}).get(f) if data else None
    return flat


def _collect_complex(pb_id: str, model: str) -> dict[str, Any]:
    json_path = repo_root() / "data" / "metrics" / model / f"{pb_id}.json"
    cif_path = repo_root() / "data" / "structures" / model / f"{pb_id}.cif"
    data = _read_json(json_path)
    fields = COMPLEX_FIELDS_BASE + EXTRA_FIELDS.get(model, ())
    row = _flatten(pb_id, model, data, fields)
    row[f"{PREFIX[model]}_cif_exists"] = cif_path.exists()
    return row


# ---------------------------------------------------------------------------
# pb_* merge helper
# ---------------------------------------------------------------------------


def _pb_columns(designs: pd.DataFrame) -> list[str]:
    """Pick the ``pb_*`` columns to merge into grand_metrics.

    Drops URL / mirror-pointer columns listed in ``PB_DROP_COLUMNS``.
    ``pb_id`` is the join key — handled separately.
    """
    keep = [
        c
        for c in designs.columns
        if c.startswith("pb_") and c != "pb_id" and c not in PB_DROP_COLUMNS
    ]
    return keep


# ---------------------------------------------------------------------------
# Scoring-layer collectors
# ---------------------------------------------------------------------------


def _iter_scoring_json(
    layer_dir: Path, pb_ids: Iterable[str]
) -> Iterable[tuple[str, dict[str, Any] | None]]:
    for pb_id in pb_ids:
        yield pb_id, _read_json(layer_dir / f"{pb_id}.json")


def _discover_per_model_fields(
    layer_dir: Path,
    pb_ids: Iterable[str],
    sources: tuple[str, ...],
) -> dict[str, list[str]]:
    """Walk every JSON in ``layer_dir`` and return the union of sub-fields
    seen per source model. Stable sort order so column ordering is
    deterministic across runs."""
    fields_by_source: dict[str, set[str]] = {s: set() for s in sources}
    if not layer_dir.exists():
        return {s: [] for s in sources}
    for _pb_id, data in _iter_scoring_json(layer_dir, pb_ids):
        if not isinstance(data, dict):
            continue
        for source in sources:
            sub = data.get(source)
            if isinstance(sub, dict):
                fields_by_source[source].update(sub.keys())
    return {s: sorted(fields_by_source[s]) for s in sources}


def _discover_flat_fields(
    layer_dir: Path, pb_ids: Iterable[str]
) -> list[str]:
    """Sequence-only layers — no per-model nesting. Union of top-level keys."""
    seen: set[str] = set()
    if not layer_dir.exists():
        return []
    for _pb_id, data in _iter_scoring_json(layer_dir, pb_ids):
        if isinstance(data, dict):
            seen.update(data.keys())
    return sorted(seen)


def _collect_per_model_scoring(
    pb_id: str,
    layer: str,
    layer_dir: Path,
    sources: tuple[str, ...],
    fields_by_source: dict[str, list[str]],
) -> dict[str, Any]:
    """Flatten a per-source-model scoring JSON onto a single row.

    Column naming: ``{layer}_{source_sub_prefix}_{field}``. e.g.
    ``prodigy_b2_dG``, ``destress_esmfold_rosetta_total``.
    """
    data = _read_json(layer_dir / f"{pb_id}.json")
    row: dict[str, Any] = {}
    for source in sources:
        sub_prefix = SCORING_SOURCE_PREFIX[source]
        sub = (data or {}).get(source) if isinstance(data, dict) else None
        for field in fields_by_source.get(source, []):
            value = sub.get(field) if isinstance(sub, dict) else None
            row[f"{layer}_{sub_prefix}_{field}"] = value
    return row


def _collect_flat_scoring(
    pb_id: str,
    layer: str,
    layer_dir: Path,
    fields: list[str],
) -> dict[str, Any]:
    """Flatten a sequence-only scoring JSON — no per-model nesting."""
    data = _read_json(layer_dir / f"{pb_id}.json")
    row: dict[str, Any] = {}
    for field in fields:
        value = data.get(field) if isinstance(data, dict) else None
        row[f"{layer}_{field}"] = value
    return row


# ---------------------------------------------------------------------------
# Consensus
# ---------------------------------------------------------------------------


def _add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Cross-folder consensus columns. Match rescoring stack conventions:
    soft thresholds (ipSAE >= 0.4, iPTM >= 0.7) and a count over the
    four complex models we re-ran here."""
    ipsae_cols = [f"{PREFIX[m]}_ipsae_d0chn_max" for m in COMPLEX_MODELS]
    iptm_cols = [f"{PREFIX[m]}_iptm" for m in COMPLEX_MODELS]
    for c in ipsae_cols + iptm_cols:
        if c not in df.columns:
            df[c] = pd.NA
    df["ipsae_pass_4folders"] = (df[ipsae_cols].apply(pd.to_numeric, errors="coerce") >= 0.4).sum(
        axis=1
    )
    df["iptm_pass_4folders"] = (df[iptm_cols].apply(pd.to_numeric, errors="coerce") >= 0.7).sum(
        axis=1
    )
    return df


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default="data/grand_metrics.csv",
        help="Output CSV path, relative to repo root.",
    )
    args = parser.parse_args(argv)

    root = repo_root()
    designs = pd.read_csv(root / "data" / "designs.csv")
    pb_ids = designs["pb_id"].astype(str).tolist()

    # Discover scoring-layer fields once across all designs so the column
    # set is stable. If the directory doesn't exist yet, we get an empty
    # list and the layer's column block is just not added (idempotent).
    metrics_dir = root / "data" / "metrics"
    prodigy_dir = metrics_dir / "prodigy"
    destress_dir = metrics_dir / "destress"
    esm_pll_dir = metrics_dir / "esm_pll"
    netsolp_dir = metrics_dir / "netsolp"

    prodigy_fields = _discover_per_model_fields(prodigy_dir, pb_ids, PRODIGY_SOURCES)
    destress_fields = _discover_per_model_fields(destress_dir, pb_ids, DESTRESS_SOURCES)
    esm_pll_fields = _discover_flat_fields(esm_pll_dir, pb_ids)
    netsolp_fields = _discover_flat_fields(netsolp_dir, pb_ids)

    rows: list[dict[str, Any]] = []
    for pb_id in pb_ids:
        row: dict[str, Any] = {"pb_id": pb_id}
        row.update(_collect_typer(pb_id))
        for model in COMPLEX_MODELS:
            row.update(_collect_complex(pb_id, model))
        # Scoring layers (per-source-model nesting).
        row.update(
            _collect_per_model_scoring(
                pb_id, "prodigy", prodigy_dir, PRODIGY_SOURCES, prodigy_fields
            )
        )
        row.update(
            _collect_per_model_scoring(
                pb_id, "destress", destress_dir, DESTRESS_SOURCES, destress_fields
            )
        )
        # Sequence-only scoring layers.
        row.update(_collect_flat_scoring(pb_id, "esm_pll", esm_pll_dir, esm_pll_fields))
        row.update(_collect_flat_scoring(pb_id, "netsolp", netsolp_dir, netsolp_fields))
        rows.append(row)

    df = pd.DataFrame(rows)

    # Merge pb_* (drop URL/mirror-pointer ones) and the small identity
    # block from designs.csv.
    pb_cols = _pb_columns(designs)
    merge_cols = ["pb_id", "design_id", "sequence_length", "is_binder", "is_strong", *pb_cols]
    df = df.merge(designs[merge_cols], on="pb_id", how="left")

    df = _add_consensus(df)

    leading = ["design_id", "pb_id", "sequence_length", "is_binder", "is_strong"]
    cols = leading + [c for c in df.columns if c not in leading]
    df = df[cols]

    out_path = root / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    by_model = {m: int((df[f"{PREFIX[m]}_status"] == "ok").sum()) for m in COMPLEX_MODELS}
    print(
        f"Wrote {out_path}  ({len(df)} rows x {len(df.columns)} cols)\n"
        f"  ok per model: {by_model}\n"
        f"  ipsae_pass_4folders >= 2: {(df['ipsae_pass_4folders'] >= 2).sum()}\n"
        f"  iptm_pass_4folders  >= 2: {(df['iptm_pass_4folders'] >= 2).sum()}\n"
        f"  pb_* columns merged: {len(pb_cols)}\n"
        f"  scoring layers — prodigy fields: "
        f"{ {s: len(fs) for s, fs in prodigy_fields.items()} }\n"
        f"  scoring layers — destress fields: "
        f"{ {s: len(fs) for s, fs in destress_fields.items()} }\n"
        f"  esm_pll fields: {len(esm_pll_fields)}, netsolp fields: {len(netsolp_fields)}"
    )


if __name__ == "__main__":
    main()
