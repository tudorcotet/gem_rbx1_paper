"""Leaderboard — top binders by KD, per team, per method family.

Outputs:
- analyses/leaderboard/top_binders.csv
- analyses/leaderboard/team_winners.csv
- analyses/leaderboard/method_winners.csv
- analyses/leaderboard/summary.json
- analyses/leaderboard/report.md
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.utils import load_designs

OUT = Path(__file__).resolve().parent


def _to_md(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except ImportError:
        cols = list(df.columns)
        head = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        body = "\n".join(
            "| " + " | ".join(str(row[c]) for c in cols) + " |"
            for _, row in df.iterrows()
        )
        return "\n".join([head, sep, body])

_BINDER_COLS = [
    "design_id", "pb_id", "name", "team_id", "binding_strength",
    "kd_nM_mean", "kd_nM_min", "kd_nM_max", "pkd_arith_mean",
    "koff_mean", "kon_mean",
    "method_family",
    "pb_design_class", "sequence_length",
    "pb_seqidentity_afdb50", "pb_tm_score_afdb50",
    "pb_esmfold_plddt", "pb_novelty",
]


def main() -> None:
    df = load_designs(only_binders=True)
    cols = [c for c in _BINDER_COLS if c in df.columns]
    binders = df[cols].sort_values("kd_nM_mean", na_position="last").reset_index(drop=True)
    binders.to_csv(OUT / "top_binders.csv", index=False)

    team = (
        binders.groupby("team_id", as_index=False)
        .agg(
            best_kd_nM=("kd_nM_mean", "min"),
            n_binders=("design_id", "count"),
            best_design_id=("design_id", "first"),
            best_pb_id=("pb_id", "first"),
            best_method=("method_family", "first"),
        )
        .sort_values("best_kd_nM")
    )
    team.to_csv(OUT / "team_winners.csv", index=False)

    method = (
        binders.groupby("method_family", as_index=False)
        .agg(
            best_kd_nM=("kd_nM_mean", "min"),
            n_binders=("design_id", "count"),
            best_design_id=("design_id", "first"),
            best_pb_id=("pb_id", "first"),
        )
        .sort_values("best_kd_nM")
    )
    method.to_csv(OUT / "method_winners.csv", index=False)

    summary = {
        "n_binders": int(len(binders)),
        "best_kd_nM": float(binders.iloc[0]["kd_nM_mean"]) if len(binders) else None,
        "best_team": str(binders.iloc[0]["team_id"]) if len(binders) else None,
        "best_method_family": str(binders.iloc[0]["method_family"]) if len(binders) else None,
        "best_pb_id": str(binders.iloc[0]["pb_id"]) if len(binders) else None,
        "median_kd_nM_binders": float(binders["kd_nM_mean"].median()) if len(binders) else None,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    leader_md = binders[[
        "design_id", "pb_id", "team_id", "binding_strength",
        "kd_nM_mean", "method_family", "pb_design_class", "sequence_length",
    ]].rename(columns={
        "design_id": "id",
        "pb_id": "slug",
        "binding_strength": "strength",
        "kd_nM_mean": "KD (nM)",
        "method_family": "method",
        "pb_design_class": "modality",
        "sequence_length": "len",
    }).round(1)

    lines = [
        "# RBX1 — leaderboard",
        "",
        f"**Headline:** {summary['n_binders']} binders. "
        f"Tightest is {summary['best_team']} via {summary['best_method_family']} "
        f"at {summary['best_kd_nM']:.1f} nM (slug `{summary['best_pb_id']}`).",
        "",
        "## Top binders, ordered by KD",
        "",
        _to_md(leader_md),
        "",
        "## Best KD per team",
        "",
        _to_md(team.round(1)),
        "",
        "## Best KD per method family",
        "",
        _to_md(method.round(1)),
        "",
        "## Caveats",
        "",
        "- 9 binders. Per-team and per-method ranks are anecdotal.",
        "- One Weak binder has a fittable curve; it sits at the bottom of the KD",
        "  list. Filter on `binding_strength == 'Strong'` for an upper-tier slice.",
        "- The Zn²⁺-buffer rerun may add binders. Re-run after the next data drop.",
        "- KD values are mean across replicates with `fixed=True AND excluded=False`",
        "  applied upstream by ProteinBase. Re-aggregate only if you document a",
        "  different filter.",
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines))
    print(f"[leaderboard] wrote {OUT/'top_binders.csv'}, {OUT/'team_winners.csv'}, "
          f"{OUT/'method_winners.csv'}")
    print(f"[leaderboard] wrote {OUT/'summary.json'}, {OUT/'report.md'}")


if __name__ == "__main__":
    main()
