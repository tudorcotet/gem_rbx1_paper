"""
TEMPLATE — copy this into analyses/<your_handle>/ and edit.

What this analysis answers (replace this line): one-sentence framing.

Run:
    mise run analysis:<your_handle>
    # or:
    uv run python analyses/<your_handle>/main.py

Inputs:
- data/designs.csv via scripts.utils.load_designs() — ALWAYS go through
  the loader, never a hard-coded path.

Outputs (written inside this folder only):
- analyses/<your_handle>/report.md    human-readable summary, ≤1 page
- analyses/<your_handle>/summary.json  machine-readable headline numbers
"""
from __future__ import annotations

import json
from pathlib import Path

from scripts.utils import load_designs

OUT = Path(__file__).resolve().parent


def main() -> None:
    # 1. Load the canonical data. ALWAYS via load_designs() — never a
    #    hard-coded path. Filter at load time when you can.
    df = load_designs()

    # 2. Compute headline numbers.
    is_binder = df["is_binder"].astype("boolean").fillna(False)
    is_expressed = df["is_expressed"].astype("boolean").fillna(False)
    summary = {
        "n_designs": int(len(df)),
        "n_expressed": int(is_expressed.sum()),
        "n_binders": int(is_binder.sum()),
        "hit_rate_overall": float(is_binder.mean()),
    }

    # 3. Write the machine-readable summary.
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # 4. Write a human-readable report (markdown).
    lines = [
        "# <Replace with your analysis title>",
        "",
        f"**Headline:** {summary['n_binders']} binders / {summary['n_designs']} designs "
        f"({summary['hit_rate_overall']:.1%} hit rate).",
        "",
        "## Method",
        "",
        "Briefly: what was computed, which filter was applied, which stat test was used.",
        "Cite the column names from `data/designs.csv` so the reader can re-derive.",
        "",
        "## Results",
        "",
        "```",
        json.dumps(summary, indent=2),
        "```",
        "",
        "## Caveats",
        "",
        "- Small N: 9 binders across 322 designs. Per-method or per-modality splits",
        "  collapse to anecdotes fast.",
        "- Modality is imbalanced: `Other` is 51% of the cohort. Don't claim a",
        "  modality trend without addressing it.",
        "- Default-buffer screen only. The Zn²⁺ rerun lives in a separate Foundry",
        "  experiment; numbers may shift when it lands.",
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines))
    print(f"[template] wrote {OUT/'report.md'} and {OUT/'summary.json'}")


if __name__ == "__main__":
    main()
