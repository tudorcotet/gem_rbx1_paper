"""Re-run ProteinTyper on designs missing predictions in the public release.

Submits each sequence to the Modal ProteinTyper endpoint with the **default
recipe** — matching exactly how ProteinBase invokes typer (no explicit
`recipe`, no `target`, no webhook). Polls the retrieve endpoint until the
job lands, then writes:

    data/metrics/proteintyper/<pb_id>.json   — full TyperJobOutput
    data/structures/esmfold/<pb_id>.cif      — ESMFold prediction
    data/images/<pb_id>.png                  — stylised render

The four targets are the 4 designs in the public RBX1 release that shipped
without ESMFold predictions: `lunar-falcon-lotus`, `violet-raven-ice`,
`radiant-vole-pearl`, `wild-zebra-cloud`.

Auth:
    Reads ``PROTEINTYPER_API_TOKEN`` from the environment — set it to the
    bearer token your ProteinTyper deployment expects.

Endpoints:
    ``PROTEINTYPER_SUBMIT_URL`` and ``PROTEINTYPER_RETRIEVE_URL`` env vars
    are mandatory — point them at your ProteinTyper service. The hosted
    version that built the public RBX1 release isn't advertised here
    because the bearer token is service-specific.

Run via::

    mise run rerun:typer
    # or:
    uv run python scripts/data/run_proteintyper.py

The script is idempotent — designs whose typer output already exists on
disk are skipped. Pass `--force` to re-submit anyway.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.utils.load_data import repo_root

# Endpoints are mandatory ENV vars so the public repo doesn't advertise any
# particular vendor's hosted typer. Point them at your ProteinTyper service
# (the original ProteinBase release used a hosted ProteinTyper Modal app —
# see the typer's docs for endpoint URLs).
SUBMIT_URL = os.environ.get("PROTEINTYPER_SUBMIT_URL")
RETRIEVE_URL = os.environ.get("PROTEINTYPER_RETRIEVE_URL")
POLL_INTERVAL_S = 30
POLL_TIMEOUT_S = 60 * 60  # 60 min hard cap per design (cold start + full_monomer recipe)
INITIAL_DELAY_S = 60      # don't poll for the first minute; typer level 0 alone is ~60-90s


def _auth_token() -> str:
    token = os.environ.get("PROTEINTYPER_API_TOKEN")
    if not token:
        sys.stderr.write(
            "[typer] PROTEINTYPER_API_TOKEN is unset. Export the bearer token "
            "your ProteinTyper service accepts before running.\n"
        )
        sys.exit(1)
    if not SUBMIT_URL or not RETRIEVE_URL:
        sys.stderr.write(
            "[typer] PROTEINTYPER_SUBMIT_URL / PROTEINTYPER_RETRIEVE_URL must be set "
            "to your ProteinTyper service endpoints.\n"
        )
        sys.exit(1)
    return token


def _post_json(url: str, payload: dict[str, Any], token: str, timeout: int = 60) -> dict[str, Any]:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {token}",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _submit(sequence: str, pb_id: str, token: str) -> str:
    """Submit a sequence to Modal proteintyper-submit. Returns the typer_key.

    Payload mirrors ProteinBase's workers/proteintyper-receiver: no `target`,
    no `webhook_*`. We pass `recipe.template = "full_monomer"` explicitly to
    request the same set of metrics ProteinBase ships in the public release
    (esmfold + proteinmpnn + novelty + classification + domainmatch). The
    OpenAPI default for `requested_outputs` is only 2 metrics; the docs claim
    `full_monomer` is the default but the schema disagrees, so we pin it.
    """
    payload = {
        "sequence": sequence,
        "client_id": pb_id,
        "start_typing_job": True,
        "recipe": {"template": "full_monomer"},
    }
    resp = _post_json(SUBMIT_URL, payload, token)
    sequence_id = resp.get("sequence_id")
    if not sequence_id:
        raise RuntimeError(f"submit returned no sequence_id: {resp}")
    print(
        f"[typer] {pb_id}: submitted  sequence_id={sequence_id}  "
        f"job_started={resp.get('job_started')}  "
        f"requested_outputs={len(resp.get('requested_outputs') or [])} metrics"
    )
    return sequence_id


_AFDB50_METRIC_NAMES = {
    "novelty", "seqidentity_afdb50", "evalue_afdb50", "tm_score_afdb50",
}


def _has_afdb50_metric(obj: Any) -> bool:
    if isinstance(obj, dict):
        if isinstance(obj.get("metric_type"), dict) and (obj["metric_type"] or {}).get(
            "name"
        ) in _AFDB50_METRIC_NAMES:
            return True
        return any(_has_afdb50_metric(v) for v in obj.values())
    if isinstance(obj, list):
        return any(_has_afdb50_metric(v) for v in obj)
    return False


def _full_monomer_complete(parsed: dict[str, Any]) -> bool:
    """Heuristic check that every async stage of the `full_monomer` recipe
    landed — not just ESMFold (level 0).

    The deployed typer streams its result in stages: ESMFold CIF first, then
    AFDB50 foldseek (novelty / seqidentity / evalue / tm_score), TED domain
    matching (`proteindomains`), and a stylised PNG render. The old exit
    condition fired the moment the ESMFold CIF was attached, before any of
    the downstream metrics were folded in. We now also require at least one
    AFDB50 metric OR a non-empty `proteindomains` OR a PNG render url.
    """
    structures = parsed.get("structures") or []
    if not structures or not isinstance(structures[0], dict):
        return False
    first = structures[0]
    if not first.get("file"):
        return False
    if _has_afdb50_metric(parsed):
        return True
    if parsed.get("proteindomains"):
        return True
    if (first.get("img") or {}).get("url"):
        return True
    return False


def _retrieve(typer_key: str, token: str) -> dict[str, Any] | None:
    """One retrieve call. Returns the parsed TyperJobOutput when ready, else None.

    Modal returns 500 during cold start and while the job is still running
    — both are transient. 404 / 425 / 500 → None (poll again). The response
    body is `{value: "<json-encoded TyperJobOutput>", schema, client_id}` —
    so `value` arrives as a string that needs an extra json.loads.
    """
    try:
        resp = _post_json(RETRIEVE_URL, {"typer_key": typer_key}, token, timeout=30)
    except urllib.error.HTTPError as e:
        if e.code in (404, 425, 500, 502, 503, 504):
            return None  # not ready
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"retrieve HTTP {e.code}: {body}") from e

    raw = resp.get("value")
    if not raw:
        return None
    parsed = json.loads(raw) if isinstance(raw, str) else raw
    if not isinstance(parsed, dict):
        return None
    if not _full_monomer_complete(parsed):
        return None
    return parsed


def _poll_until_ready(typer_key: str, token: str, pb_id: str) -> dict[str, Any]:
    deadline = time.time() + POLL_TIMEOUT_S
    if INITIAL_DELAY_S:
        time.sleep(INITIAL_DELAY_S)
    waited = INITIAL_DELAY_S
    while time.time() < deadline:
        out = _retrieve(typer_key, token)
        if out is not None:
            print(f"[typer] {pb_id}: ready after {waited}s "
                  f"({len([k for k,v in out.items() if v is not None])} fields)")
            return out
        time.sleep(POLL_INTERVAL_S)
        waited += POLL_INTERVAL_S
        if waited % 120 == 0:
            print(f"[typer] {pb_id}: still running after {waited}s")
    raise RuntimeError(f"typer timed out after {POLL_TIMEOUT_S}s for {pb_id}")


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "rbx1_gem_paper/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as out:
        out.write(r.read())


_S3_PUBLIC_BUCKET = "s3://proteinbase-pub/"
_S3_PUBLIC_HTTPS = "https://proteinbase-pub.t3.storage.dev/"


def _s3_to_https(url: str | None) -> str | None:
    """Convert `s3://proteinbase-pub/<key>` to the public HTTPS endpoint."""
    if not isinstance(url, str):
        return None
    if url.startswith(_S3_PUBLIC_BUCKET):
        return _S3_PUBLIC_HTTPS + url[len(_S3_PUBLIC_BUCKET):]
    return url if url.startswith("http") else None


def _first_esmfold_structure(typer_output: dict[str, Any]) -> dict[str, Any] | None:
    for s in typer_output.get("structures") or []:
        if isinstance(s, dict) and (s.get("method") or {}).get("name", "").lower() == "esmfold":
            return s
    # Fallback: take the first structure if none is labelled ESMFold.
    structures = typer_output.get("structures") or []
    return structures[0] if structures else None


def _save_artifacts(typer_output: dict[str, Any], pb_id: str, data_dir: Path) -> dict[str, Path]:
    """Save the TyperJobOutput JSON + download the CIF / PNG artifacts.

    Returns the local paths in a dict.
    """
    out: dict[str, Path] = {}

    out_json = data_dir / "metrics" / "proteintyper" / f"{pb_id}.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(typer_output, indent=2) + "\n")
    out["typer_output"] = out_json

    esmfold = _first_esmfold_structure(typer_output) or {}

    cif_url = _s3_to_https((esmfold.get("file") or {}).get("url"))
    if cif_url:
        dest = data_dir / "structures" / "esmfold" / f"{pb_id}.cif"
        _download(cif_url, dest)
        out["esmfold_cif"] = dest

    png_url = _s3_to_https((esmfold.get("img") or {}).get("url"))
    if png_url:
        dest = data_dir / "images" / f"{pb_id}.png"
        _download(png_url, dest)
        out["stylized_png"] = dest

    return out


def _designs_missing_predictions(df: pd.DataFrame, structures_root: Path) -> pd.DataFrame:
    """Designs that need a rerun.

    Truth source is the filesystem (`data/metrics/proteintyper/<pb_id>.json`),
    not `has_predictions` in `designs.csv` — that flag flips True as soon as a
    local rerun lands, so checking it would make the script non-recoverable
    if someone deleted the local JSON.
    """
    json_dir = structures_root / "metrics" / "proteintyper"
    have = {p.stem for p in json_dir.glob("*.json")} if json_dir.exists() else set()
    pb_release_done = df["pb_predictions_source"] == "proteinbase_release" if "pb_predictions_source" in df.columns else df["has_predictions"].astype(str).str.lower() == "true"
    need_rerun = ~pb_release_done & ~df["pb_id"].isin(have)
    return df[need_rerun].copy()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true",
        help="Resubmit even if data/structures/typer_outputs/<pb_id>.json exists",
    )
    parser.add_argument(
        "--pb-id", action="append", default=None,
        help="Run on a specific pb_id instead of the missing-predictions set. Repeatable.",
    )
    args = parser.parse_args(argv)

    token = _auth_token()
    root = repo_root()
    df = pd.read_csv(root / "data" / "designs.csv")

    if args.pb_id:
        targets = df[df["pb_id"].isin(args.pb_id)]
    else:
        targets = _designs_missing_predictions(df, root / "data")

    if targets.empty:
        print("[typer] no targets — every design has predictions on disk.")
        return

    data_dir = root / "data"
    print(f"[typer] running on {len(targets)} designs")
    for _, row in targets.iterrows():
        pb_id = row["pb_id"]
        sequence = row["sequence"]
        if not isinstance(sequence, str) or not sequence:
            print(f"[typer] {pb_id}: no sequence, skipping")
            continue

        existing = data_dir / "metrics" / "proteintyper" / f"{pb_id}.json"
        if existing.exists() and not args.force:
            print(f"[typer] {pb_id}: already has output at {existing}, skipping")
            continue

        try:
            typer_key = _submit(sequence, pb_id, token)
            output = _poll_until_ready(typer_key, token, pb_id)
            paths = _save_artifacts(output, pb_id, data_dir)
            wrote = ", ".join(p.name for p in paths.values())
            print(f"[typer] {pb_id}: saved {wrote}")
        except Exception as e:
            print(f"[typer] {pb_id}: FAILED — {e}")


if __name__ == "__main__":
    main()
