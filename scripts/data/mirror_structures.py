"""Mirror every ProteinBase structure artifact into the repo.

Reads `data/designs.csv`, follows `pb_esmfold_cif_url`, `pb_stylized_png_url`,
`pb_spr_curves_urls`, `pb_bli_curves_urls` for every row that has them, and
writes the bytes to:

    data/structures/esmfold/<pb_id>.cif         monomer (ESMFold) — under structures/
    data/images/<pb_id>.png                     stylised renders
    data/sensorgrams/<pb_id>_rep<NN>_{spr,bli}.json   SPR (primary) + BLI (fallback) curves

Sensorgrams + images sit alongside `structures/`, not inside it: a
sensorgram is a kinetic trace, not a structure, and a stylised PNG is
a render, not a model output.

Idempotent: skips files already on disk. Safe to re-run after a fresh build.

Run via::

    mise run mirror:structures
    # or:
    uv run python scripts/data/mirror_structures.py
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import pandas as pd

from scripts.utils.load_data import repo_root


def _download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "rbx1_gem_paper/0.1"})
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with urllib.request.urlopen(req, timeout=60) as r, tmp.open("wb") as out:
        out.write(r.read())
    tmp.replace(dest)
    return True


def main(only_binders: bool = False) -> None:
    root = repo_root()
    df = pd.read_csv(root / "data" / "designs.csv")
    if only_binders:
        df = df[df["is_binder"].astype(str).str.lower().eq("true")]

    data_dir = root / "data"
    counts = {"cif": [0, 0], "png": [0, 0], "spr": [0, 0], "bli": [0, 0]}

    for _, row in df.iterrows():
        pb_id = row.get("pb_id")
        if not isinstance(pb_id, str) or not pb_id:
            continue

        cif_url = row.get("pb_esmfold_cif_url")
        if isinstance(cif_url, str) and cif_url:
            try:
                wrote = _download(cif_url, data_dir / "structures" / "esmfold" / f"{pb_id}.cif")
                counts["cif"][0 if wrote else 1] += 1
            except Exception as e:
                print(f"[mirror] FAILED CIF {pb_id}: {e}")

        png_url = row.get("pb_stylized_png_url")
        if isinstance(png_url, str) and png_url:
            try:
                wrote = _download(png_url, data_dir / "images" / f"{pb_id}.png")
                counts["png"][0 if wrote else 1] += 1
            except Exception as e:
                print(f"[mirror] FAILED PNG {pb_id}: {e}")

        for col, assay in (("pb_spr_curves_urls", "spr"), ("pb_bli_curves_urls", "bli")):
            urls = row.get(col)
            if not isinstance(urls, str) or not urls:
                continue
            for i, url in enumerate(urls.split("|"), start=1):
                dest = data_dir / "sensorgrams" / f"{pb_id}_rep{i:02d}_{assay}.json"
                try:
                    wrote = _download(url, dest)
                    counts[assay][0 if wrote else 1] += 1
                except Exception as e:
                    print(f"[mirror] FAILED {assay} {pb_id}: {e}")

    for k, (new, skip) in counts.items():
        print(f"[mirror] {k}: wrote={new}  already_there={skip}")
    print(f"[mirror] root: {data_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binders", action="store_true",
                        help="Only mirror the 9 binders (default: all 322)")
    args = parser.parse_args()
    main(only_binders=args.binders)
