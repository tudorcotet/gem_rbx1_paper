"""Canonical data loaders. Every analysis script imports from here.

Why this exists: every analysis used to hard-code a `Path(__file__).parents[N]`
prefix. When the repo layout drifts (someone moves a file, renames a folder),
every script breaks at once. Funnel all loads through this module instead.

Usage::

    from scripts.utils import load_designs
    df = load_designs()                       # 322 rows
    df = load_designs(only_binders=True)      #   9 rows
    df = load_designs(only_expressed=True)    # 254 rows
    df = load_designs(only_strong=True)       #   1 row
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    """Walk up from this file to the directory that owns `pyproject.toml`."""
    here = Path(__file__).resolve()
    for parent in [here, *here.parents]:
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError(f"Could not locate repo root from {here}")


def _designs_path() -> Path:
    root = repo_root()
    parquet = root / "data" / "designs.parquet"
    if parquet.exists():
        return parquet
    return root / "data" / "designs.csv"


_BOOL_COLUMNS = (
    "is_binder",
    "is_strong",
    "is_expressed",
    "is_control",
    "any_binding",
)


def _coerce_bools(df: pd.DataFrame) -> pd.DataFrame:
    for col in _BOOL_COLUMNS:
        if col not in df.columns:
            continue
        s = df[col]
        if s.dtype == bool:
            df[col] = df[col].astype("boolean")
            continue
        df[col] = s.map(
            {
                True: True,
                False: False,
                "True": True,
                "False": False,
                "true": True,
                "false": False,
                1: True,
                0: False,
                "1": True,
                "0": False,
            }
        ).astype("boolean")
    return df


@lru_cache(maxsize=8)
def load_designs(
    *,
    only_expressed: bool = False,
    only_binders: bool = False,
    only_strong: bool = False,
    drop_controls: bool = True,
) -> pd.DataFrame:
    """Load the canonical 322-row design table.

    Bool-typed columns (`is_binder`, `is_expressed`, …) are cast to pandas
    nullable `boolean` so you can write `df.is_binder & df.is_expressed`
    without object-dtype warnings.

    Args:
        only_expressed: keep only the designs that expressed (`is_expressed`).
        only_binders:   keep only confirmed binders (`is_binder`). Implies
                        `only_expressed`.
        only_strong:    keep only Strong-bin binders (KD < 100 nM,
                        `is_strong`). Implies the two above.
        drop_controls:  drop platform-authored controls (rows where
                        `is_control=True`). On by default — the 322-design
                        headline number is controls-out.
    """
    path = _designs_path()
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)

    df = _coerce_bools(df)

    if drop_controls and "is_control" in df.columns:
        df = df[~df["is_control"].fillna(False)]

    if only_strong:
        df = df[df["is_strong"].fillna(False)]
    elif only_binders:
        df = df[df["is_binder"].fillna(False)]
    elif only_expressed:
        df = df[df["is_expressed"].fillna(False)]

    return df.reset_index(drop=True)


@lru_cache(maxsize=4)
def load_designs_fasta() -> dict[str, str]:
    """Return a `{design_id: sequence}` dict.

    Reads from `data/designs.fasta` if present (canonical source), otherwise
    falls back to the `sequence` column of `designs.csv`.
    """
    root = repo_root()
    fasta = root / "data" / "designs.fasta"
    if fasta.exists():
        out: dict[str, str] = {}
        ident, buf = None, []
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                if ident is not None:
                    out[ident] = "".join(buf)
                ident = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if ident is not None:
            out[ident] = "".join(buf)
        return out

    df = load_designs(drop_controls=False)
    return dict(zip(df["design_id"].astype(str), df["sequence"]))


def structure_path(pb_id: str, *, model: str = "esmfold") -> Path | None:
    """Return the local CIF path for `pb_id` under the given model, if it exists.

    Layout convention:
        data/structures/esmfold/<pb_id>.cif      monomer (ProteinTyper)
        data/structures/boltz2/<pb_id>.cif       complex (rescoring stack)
        data/structures/protenix/<pb_id>.cif     complex (Protenix v2)
        data/structures/chai/<pb_id>.cif         complex (Chai)
    """
    p = repo_root() / "data" / "structures" / model / f"{pb_id}.cif"
    return p if p.exists() else None


def image_path(pb_id: str) -> Path | None:
    """Return the stylised render PNG path for `pb_id`, if it exists."""
    p = repo_root() / "data" / "images" / f"{pb_id}.png"
    return p if p.exists() else None


def sensorgram_paths(pb_id: str) -> list[Path]:
    """Return all SPR / BLI sensorgram JSONs for `pb_id`.

    Pattern: `data/sensorgrams/<pb_id>_rep<NN>_{spr,bli}.json`.
    Returns an empty list when nothing has been mirrored.
    """
    sgs = sorted((repo_root() / "data" / "sensorgrams").glob(f"{pb_id}_rep*.json"))
    return list(sgs)


def metrics_path(pb_id: str, *, model: str = "proteintyper") -> Path | None:
    """Return the raw per-model output JSON path for `pb_id` (e.g. full TyperJobOutput)."""
    p = repo_root() / "data" / "metrics" / model / f"{pb_id}.json"
    return p if p.exists() else None


def binding_strength_palette() -> dict[str, str]:
    """Brand palette for the binding-strength hierarchy.

    Loaded from `theme/palettes.json` so figures stay in sync with the rest
    of the visual system.
    """
    import json

    pal_file = repo_root() / "theme" / "palettes.json"
    if pal_file.exists():
        pal = json.loads(pal_file.read_text())
        if "binding_strength_2026" in pal.get("palettes", {}):
            return pal["palettes"]["binding_strength_2026"]
        if "binding_strength_2026" in pal:
            return pal["binding_strength_2026"]
    return {
        "Binder": "#00D9FF",
        "Strong": "#30C5F5",
        "Medium": "#33C4FF",
        "Weak": "#36B7F6",
        "Non-binder": "#5C6773",
        "No expression": "#9EA2AF",
        "Unknown": "#98A2AE",
    }
