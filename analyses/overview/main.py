"""Overview analysis — headline counts, hit rates, modality x method cross-tab.

Outputs:
- analyses/overview/summary.json
- analyses/overview/report.md
- analyses/overview/modality_x_method.csv
- analyses/overview/hit_rate_by_method.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.utils import load_designs


def _to_md(df: pd.DataFrame, index: bool = False) -> str:
    """Render a DataFrame as a markdown table, falling back if tabulate is missing."""
    try:
        return df.to_markdown(index=index)
    except ImportError:
        cols = (([df.index.name or ""] if index else []) + list(df.columns))
        rows = []
        for idx, row in df.iterrows():
            cells = ([str(idx)] if index else []) + [str(row[c]) for c in df.columns]
            rows.append(cells)
        head = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        body = "\n".join("| " + " | ".join(r) + " |" for r in rows)
        return "\n".join([head, sep, body])

OUT = Path(__file__).resolve().parent


def main() -> None:
    df = load_designs()
    is_binder = df["is_binder"].astype("boolean").fillna(False)
    is_expressed = df["is_expressed"].astype("boolean").fillna(False)
    is_strong = df["is_strong"].astype("boolean").fillna(False)

    summary = {
        "n_designs": int(len(df)),
        "n_expressed": int(is_expressed.sum()),
        "expression_rate": float(is_expressed.mean()),
        "n_binders": int(is_binder.sum()),
        "n_strong": int(is_strong.sum()),
        "hit_rate_overall": float(is_binder.mean()),
        "hit_rate_among_expressed": float(is_binder[is_expressed].mean()) if is_expressed.any() else 0.0,
        "best_kd_nM": float(df.loc[is_binder, "kd_nM_mean"].min()) if is_binder.any() else None,
        "best_team": df.loc[df["kd_nM_mean"].idxmin(), "team_id"]
            if df["kd_nM_mean"].notna().any() else None,
        "n_teams": int(df["team_id"].nunique()),
        "n_methods_named": int((df["method_family"] != "Not mentioned").sum()),
    }

    # Modality x method cross-tab.
    xtab = (
        df.assign(
            modality=df["pb_design_class"].fillna("Unknown"),
            method=df["method_family"].fillna("Not mentioned"),
        )
        .groupby(["method", "modality"])
        .size()
        .unstack(fill_value=0)
        .sort_index()
    )
    xtab.to_csv(OUT / "modality_x_method.csv")

    # Hit rate by method family.
    by_method = (
        df.assign(
            method=df["method_family"].fillna("Not mentioned"),
            is_binder=is_binder.fillna(False),
            is_expressed=is_expressed.fillna(False),
        )
        .groupby("method")
        .agg(
            n=("design_id", "count"),
            n_expressed=("is_expressed", "sum"),
            n_binders=("is_binder", "sum"),
        )
        .assign(
            expression_rate=lambda d: d["n_expressed"] / d["n"],
            hit_rate=lambda d: d["n_binders"] / d["n"],
            hit_rate_among_expressed=lambda d: d["n_binders"] / d["n_expressed"].replace(0, pd.NA),
        )
        .sort_values("n", ascending=False)
    )
    by_method.to_csv(OUT / "hit_rate_by_method.csv")

    # Top binders for the report.
    top = (
        df[is_binder]
        .sort_values("kd_nM_mean", na_position="last")
        [["design_id", "pb_id", "team_id", "binding_strength",
          "kd_nM_mean", "method_family", "pb_design_class", "sequence_length"]]
        .head(10)
    )

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# RBX1 — overview",
        "",
        f"**Headline:** {summary['n_binders']} binders / {summary['n_designs']} "
        f"designs ({summary['hit_rate_overall']:.1%}). "
        f"Expression {summary['expression_rate']:.0%}. "
        f"Tightest KD = {summary['best_kd_nM']:.1f} nM ({summary['best_team']}).",
        "",
        "## Method",
        "",
        "Read every row of `data/designs.csv`. Counted by `is_binder`, `is_strong`, "
        "`is_expressed`. Cross-tabbed `method_family` × `pb_design_class`. Hit rates "
        "are reported both raw and conditional on expression.",
        "",
        "## Headline numbers",
        "",
        "```",
        json.dumps(summary, indent=2),
        "```",
        "",
        "## Hit rate by method family",
        "",
        _to_md(by_method.round(3).reset_index()),
        "",
        "## Top binders",
        "",
        _to_md(top),
        "",
        "## Caveats",
        "",
        "- 9 binders is anecdote territory for any per-method or per-modality split.",
        "  Use the raw `n` / `n_binders` columns above before quoting a rate.",
        "- `method_family` includes `Not mentioned` (93) and `Other` (58) — many of",
        "  the unnamed entries are top-ranked teams who left `designMethod` blank.",
        "  Add patterns to `_METHOD_FAMILY_PATTERNS` in",
        "  `scripts/data/build_designs.py` as the right labels surface.",
        "- The cross-tab is dominated by `pb_design_class=Other` (164/322). The",
        "  modality breakdown is informative but not balanced — comparison across",
        "  modalities inherits that imbalance.",
        "- The Zn²⁺-buffer rerun and the ovalbumin specificity panel are in flight.",
        "  Numbers here may shift when those land.",
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines))
    print(f"[overview] wrote {OUT/'summary.json'}, {OUT/'report.md'}")
    print(f"[overview] wrote {OUT/'modality_x_method.csv'}, {OUT/'hit_rate_by_method.csv'}")


if __name__ == "__main__":
    main()
