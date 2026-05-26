# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""ESM2 t33 sequence pseudo-log-likelihood for 322 RBX1 designs.

Sequence-only scorer — no structure required. For each binder sequence
in ``data/designs.csv`` (column ``sequence``, keyed by ``pb_id``) we
compute the masked-marginal PLL using ``fair-esm`` ESM2 t33 650M
(matches the kartic_metrics implementation 1:1).

Outputs live on a Modal Volume ``rbx1-rerun-results``:

* ``esm_pll/{pb_id}.json`` — per-design metrics
* ``esm_pll_summary.csv``  — aggregated table

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_esm_pll_rbx1.py
    modal run scripts/modal/modal_esm_pll_rbx1.py --download

    GPU=A100 MODAL_APP_NAME=my-esm-pll modal run --detach \\
        scripts/modal/modal_esm_pll_rbx1.py
"""

from __future__ import annotations

import contextlib
import json
import math
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A10G")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))

APP_NAME = os.environ.get("MODAL_APP_NAME", "rbx1-esm-pll")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "rbx1-rerun-results")

PREDICTOR = "esm_pll"


# ---------------------------------------------------------------------------
# Modal image — fair-esm + torch
# ---------------------------------------------------------------------------
esm_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("fair-esm", "torch>=2.0", "pandas", "tqdm")
)

app = modal.App(APP_NAME)

RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"


# ---------------------------------------------------------------------------
# Empty metrics row (so failures still have the expected columns)
# ---------------------------------------------------------------------------
def _empty_metrics() -> dict[str, float | None]:
    return {
        "length": None,
        "esm_pll_total": None,
        "esm_pll_avg": None,
    }


# ---------------------------------------------------------------------------
# Per-batch GPU function — load ESM2 t33 once, score all sequences.
# Since each PLL pass is sequential masked-marginal (cheap for 60-130 aa
# binders), batching into a single container is far cheaper than fanning
# out one container per design.
# ---------------------------------------------------------------------------
@app.function(
    image=esm_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def score_batch(pb_ids: list[str], sequences: list[str]) -> list[dict]:
    """ESM2 t33 PLL (avg + total) — verbatim algorithm from kartic source."""
    import esm
    import torch

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Loading ESM2 t33 650M on {device}")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    def pll(seq: str) -> tuple[float, float]:
        data = [("p", seq)]
        *_, batch_tokens = batch_converter(data)
        log_probs: list[float] = []
        for i in range(len(seq)):
            masked = batch_tokens.clone()
            masked[0, i + 1] = alphabet.mask_idx
            with torch.no_grad():
                tp = torch.log_softmax(model(masked.to(device))["logits"], dim=-1)
            log_probs.append(tp[0, i + 1, alphabet.get_idx(seq[i])].item())
        s = math.fsum(log_probs)
        return s, s / len(seq)

    out: list[dict] = []
    out_dir = Path(RESULTS_DIR) / PREDICTOR
    out_dir.mkdir(parents=True, exist_ok=True)

    for pb_id, seq in zip(pb_ids, sequences, strict=True):
        try:
            total, avg = pll(seq)
            row = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": "ok",
                "length": len(seq),
                "esm_pll_total": float(total),
                "esm_pll_avg": float(avg),
            }
            print(f"  {pb_id}: pll_avg={avg:.4f} (len={len(seq)})")
        except Exception as e:
            row = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": f"error: {e}",
                **_empty_metrics(),
            }
            print(f"  {pb_id}: FAILED ({e})")
        (out_dir / f"{_sanitize(pb_id)}.json").write_text(json.dumps(row))
        out.append(row)

    RESULTS_VOLUME.commit()
    return out


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
    print(f"ESM-PLL: {len(completed)} done, {len(pending_pairs)} pending")

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
    print(f"ESM-PLL: {ok}/{len(df)} succeeded.")


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
    """Fire ESM-PLL on Modal, or pull finished results to local disk."""
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
    print(f"Triggering ESM-PLL batch: {len(pb_ids)} designs (sequence-only)")
    # Use .remote() + --detach. .spawn() only schedules on `modal deploy`'d
    # apps; under `modal run` (ephemeral), spawn returns immediately but the
    # call never actually runs.
    run_batch.remote(pb_ids, seqs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-pb_id JSON files into ./data/metrics/esm_pll/."""
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
