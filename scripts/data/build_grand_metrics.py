"""Pool per-model metric JSONs into a single wide ``grand_metrics.csv``.

Phase 3 of the rerun. After ``run_rescoring.py`` has filled
``data/metrics/{boltz2,protenix,chai}/`` and ``data/structures/{…}/``
with one file per design, this script joins them by ``pb_id`` and
writes ``data/grand_metrics.csv``.

Schema convention (matches the prefixes documented in
``docs/DATA.md`` § "Re-scoring"):

* ``pb_id``                — join key
* ``tp_*``                 — ProteinTyper (monomer ESMFold etc, already
                             on the row in ``designs.csv`` as ``pb_*``;
                             we re-export the canonical names here)
* ``b2_*``                 — Boltz-2 complex (chains A=target, B=binder)
* ``px_*``                 — Protenix complex
* ``chai_*``               — Chai-1 complex

Plus a handful of derived consensus columns at the end:

* ``ipsae_pass_3folders``  — # of {b2, px, chai} with d0chn_max >= 0.4
* ``iptm_pass_3folders``   — # of {b2, px, chai} with iptm >= 0.7

The script is idempotent. Missing JSONs leave nulls; missing CIFs leave
``cif_*_exists`` = False. Re-run after pulling new data.

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
COMPLEX_MODELS: tuple[str, ...] = ("boltz2", "protenix", "chai")
PREFIX: dict[str, str] = {
    "boltz2": "b2",
    "protenix": "px",
    "chai": "chai",
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
}


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


def _add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Cross-folder consensus columns. Match rescoring stack conventions:
    soft thresholds (ipSAE >= 0.4, iPTM >= 0.7) and a count over the
    three complex models we re-ran here."""
    ipsae_cols = [f"{PREFIX[m]}_ipsae_d0chn_max" for m in COMPLEX_MODELS]
    iptm_cols = [f"{PREFIX[m]}_iptm" for m in COMPLEX_MODELS]
    for c in ipsae_cols + iptm_cols:
        if c not in df.columns:
            df[c] = pd.NA
    df["ipsae_pass_3folders"] = (df[ipsae_cols].apply(pd.to_numeric, errors="coerce") >= 0.4).sum(
        axis=1
    )
    df["iptm_pass_3folders"] = (df[iptm_cols].apply(pd.to_numeric, errors="coerce") >= 0.7).sum(
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

    rows: list[dict[str, Any]] = []
    for pb_id in designs["pb_id"].astype(str):
        row: dict[str, Any] = {"pb_id": pb_id}
        row.update(_collect_typer(pb_id))
        for model in COMPLEX_MODELS:
            row.update(_collect_complex(pb_id, model))
        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.merge(
        designs[["pb_id", "design_id", "sequence_length", "is_binder", "is_strong"]],
        on="pb_id",
        how="left",
    )
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
        f"  ipsae_pass_3folders >= 2: {(df['ipsae_pass_3folders'] >= 2).sum()}\n"
        f"  iptm_pass_3folders  >= 2: {(df['iptm_pass_3folders'] >= 2).sum()}"
    )


if __name__ == "__main__":
    main()
