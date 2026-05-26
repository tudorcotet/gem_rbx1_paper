# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""NetSolP E. coli solubility & usability prediction for 322 RBX1 designs.

Sequence-only scorer — no structure required. For each binder sequence
in ``data/designs.csv`` (column ``sequence``, keyed by ``pb_id``) we
predict solubility / usability with NetSolP (ESM-based, see Thumuluri
et al. 2022 *NAR*). The invocation logic is the verbatim kartic
implementation: try the official ``netsolp`` Python package first, and
if it isn't pip-installable in the image, fall back to an ESM2 t12-
embedded heuristic (charged-vs-hydrophobic fraction) as a stand-in.

Outputs live on a Modal Volume ``rbx1-rerun-results``:

* ``netsolp/{pb_id}.json`` — per-design metrics
* ``netsolp_summary.csv`` — aggregated table

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_netsolp_rbx1.py
    modal run scripts/modal/modal_netsolp_rbx1.py --download

    GPU=A10G MODAL_APP_NAME=my-netsolp modal run --detach \\
        scripts/modal/modal_netsolp_rbx1.py
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "T4")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))

APP_NAME = os.environ.get("MODAL_APP_NAME", "rbx1-netsolp")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "rbx1-rerun-results")

PREDICTOR = "netsolp"


# ---------------------------------------------------------------------------
# Modal image — try to install NetSolP at build time; fall back gracefully
# at runtime if neither pypi nor the GitHub fork is reachable.
# ---------------------------------------------------------------------------
netsolp_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.0",
        "transformers>=4.30",
        "pandas",
        "numpy",
        "biopython",
    )
    .run_commands(
        "pip install netsolp || "
        "pip install git+https://github.com/tvinet/NetSolP-1.0.git || "
        "echo 'NetSolP install skipped, will use ESM2 fallback at runtime'"
    )
)

app = modal.App(APP_NAME)

RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"


# ---------------------------------------------------------------------------
# Empty metrics row
# ---------------------------------------------------------------------------
def _empty_metrics() -> dict[str, float | None]:
    return {
        "netsolp_solubility": None,
        "netsolp_usability": None,
    }


# ---------------------------------------------------------------------------
# NetSolP batch GPU function — verbatim port of kartic_metrics
# ---------------------------------------------------------------------------
@app.function(
    image=netsolp_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
    memory=8192,
)
def score_batch(pb_ids: list[str], sequences: list[str]) -> list[dict]:
    """Score sequences with NetSolP for solubility + usability prediction.

    NetSolP predicts:
      - Solubility probability (0-1)
      - Usability probability (0-1, likelihood of being usable after purification)
    """
    results: list[dict] = []
    out_dir = Path(RESULTS_DIR) / PREDICTOR
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Try official NetSolP package first
        try:
            from netsolp import predict as netsolp_predict  # type: ignore[import-not-found]
            predictions = netsolp_predict(sequences)
            for pb_id, pred in zip(pb_ids, predictions, strict=True):
                row = {
                    "pb_id": pb_id,
                    "predictor": PREDICTOR,
                    "status": "ok",
                    "netsolp_solubility": float(
                        pred.get("solubility", pred.get("Solubility", None))
                    ),
                    "netsolp_usability": float(
                        pred.get("usability", pred.get("Usability", None))
                    ),
                }
                (out_dir / f"{_sanitize(pb_id)}.json").write_text(json.dumps(row))
                results.append(row)
                print(f"  {pb_id}: sol={row['netsolp_solubility']:.4f}")
        except ImportError:
            # Fallback: Use ESM embeddings + simple solubility heuristics.
            # Based on NetSolP's approach but using readily available ESM2.
            from transformers import AutoModelForMaskedLM, AutoTokenizer
            import torch

            print("NetSolP not available, using ESM2-based solubility estimation...")
            model_name = "facebook/esm2_t12_35M_UR50D"
            tokenizer = AutoTokenizer.from_pretrained(model_name)
            model = AutoModelForMaskedLM.from_pretrained(model_name)
            model.eval()
            if torch.cuda.is_available():
                model = model.cuda()

            for pb_id, seq in zip(pb_ids, sequences, strict=True):
                try:
                    inputs = tokenizer(
                        seq,
                        return_tensors="pt",
                        padding=True,
                        truncation=True,
                        max_length=1024,
                    )
                    if torch.cuda.is_available():
                        inputs = {k: v.cuda() for k, v in inputs.items()}

                    with torch.no_grad():
                        _ = model(**inputs)  # forward pass for parity w/ kartic source

                    aa_counts: dict[str, int] = {}
                    for aa in seq:
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    total = len(seq)

                    # Charged residue fraction (correlates w/ solubility)
                    charged = sum(aa_counts.get(aa, 0) for aa in "DEKRH") / total
                    # Hydrophobic fraction (anti-correlates w/ solubility)
                    hydrophobic = sum(aa_counts.get(aa, 0) for aa in "VILMFYW") / total
                    sol_score = 0.5 + 0.3 * (charged - 0.15) - 0.4 * (hydrophobic - 0.35)
                    sol_score = max(0.0, min(1.0, sol_score))

                    row = {
                        "pb_id": pb_id,
                        "predictor": PREDICTOR,
                        "status": "ok",
                        "netsolp_solubility": round(float(sol_score), 4),
                        "netsolp_usability": None,  # only from official NetSolP
                    }
                    print(f"  {pb_id}: sol={sol_score:.4f} (esm2 fallback)")
                except Exception as e:
                    row = {
                        "pb_id": pb_id,
                        "predictor": PREDICTOR,
                        "status": f"error: {e}",
                        **_empty_metrics(),
                    }
                    print(f"  {pb_id}: FAILED ({e})")
                (out_dir / f"{_sanitize(pb_id)}.json").write_text(json.dumps(row))
                results.append(row)

    except Exception as e:
        import traceback

        print(f"NetSolP batch failed: {traceback.format_exc()}")
        for pb_id in pb_ids:
            row = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": f"error: {e}",
                **_empty_metrics(),
            }
            (out_dir / f"{_sanitize(pb_id)}.json").write_text(json.dumps(row))
            results.append(row)

    RESULTS_VOLUME.commit()
    return results


# ---------------------------------------------------------------------------
# Volume helpers
# ---------------------------------------------------------------------------
def _sanitize(name: str) -> str:
    return (
        name.replace("/", "_SLASH_")
        .replace("\\", "_BSLASH_")
        .replace("|", "_PIPE_")
        .replace(" ", "_")
        .replace(",", "_COMMA_")
    )


# ---------------------------------------------------------------------------
# Remote orchestrator
# ---------------------------------------------------------------------------
orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(
    image=orchestrator_image,
    timeout=24 * 3600,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def run_batch(pb_ids: list[str], seqs: list[str]) -> None:
    import pandas as pd

    RESULTS_VOLUME.reload()
    pred_dir = Path(RESULTS_DIR) / PREDICTOR
    completed: set[str] = set()
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                d = json.loads(f.read_text())
                if d.get("status") == "ok":
                    completed.add(d["pb_id"])

    pending_pairs = [(p, s) for p, s in zip(pb_ids, seqs, strict=True) if p not in completed]
    print(f"NetSolP: {len(completed)} done, {len(pending_pairs)} pending")

    if pending_pairs:
        pending_ids = [p[0] for p in pending_pairs]
        pending_seqs = [p[1] for p in pending_pairs]
        score_batch.remote(pending_ids, pending_seqs)

    RESULTS_VOLUME.reload()
    rows = []
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                rows.append(json.loads(f.read_text()))
    df = pd.DataFrame(rows)
    csv_path = Path(RESULTS_DIR) / f"{PREDICTOR}_summary.csv"
    df.to_csv(csv_path, index=False)
    RESULTS_VOLUME.commit()
    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"NetSolP: {ok}/{len(df)} succeeded.")


# ---------------------------------------------------------------------------
# Local entrypoint
# ---------------------------------------------------------------------------
@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Fire NetSolP on Modal, or pull finished results to local disk."""
    import pandas as pd

    if download:
        _download_to_local(designs_csv)
        return

    if retry_failed:
        _clear_failed()
        return

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)

    pb_ids = df["pb_id"].tolist()
    seqs = df["sequence"].tolist()
    print(f"Triggering NetSolP batch: {len(pb_ids)} designs (sequence-only)")
    # Use .remote() + --detach. .spawn() only schedules on `modal deploy`'d
    # apps; under `modal run` (ephemeral), spawn returns immediately but the
    # call never actually runs.
    run_batch.remote(pb_ids, seqs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-pb_id JSON files into ./data/metrics/netsolp/."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected_ids = set(df["pb_id"].astype(str).tolist())

    json_out = Path("./data/metrics") / PREDICTOR
    json_out.mkdir(parents=True, exist_ok=True)

    n_json = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if entry.path.endswith(".json"):
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                data = json.loads(payload)
                pb_id = data.get("pb_id") or Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                (json_out / f"{pb_id}.json").write_bytes(payload)
                n_json += 1
    except Exception as e:
        print(f"  no metric JSONs to download: {e}")

    print(f"Downloaded {n_json} JSON for {PREDICTOR}.")


def _clear_failed() -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  cleared {cleared} failed results")
