"""Render every figure referenced in paper/main.tex.

Outputs land at ``figures/paper/<name>.{png,pdf,svg}``. Filenames match
the ``\\includegraphics{...}`` calls in ``paper/main.tex``:

    fig1_workflow                — schematic of the competition pipeline
    fig2_results                 — 4-panel results landscape
    figS1_modality_method_crosstab — heatmap
    figS2_target_region          — declared target region vs outcome
    figS3_kd_replicates          — per-binder KD with replicate spread
    figS4_complex_grand_metrics  — Boltz-2 / Chai-1 / Protenix scatter (requires
                                   ``data/metrics/{boltz2,chai,protenix}/*``;
                                   skipped with a warning if missing)

Run::

    mise run figures
    # or:
    uv run python scripts/plotting/figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.plotting import apply_theme, save_figure
from scripts.utils import binding_strength_palette, load_designs, repo_root

PAL = binding_strength_palette()
CYAN = PAL.get("Strong", "#30C5F5")
SLATE = PAL.get("Non-binder", "#5C6773")
INK = "#0F1419"


# ----------------------------------------------------------------------
# Figure 1 — schematic of the competition workflow
# ----------------------------------------------------------------------
def fig1_workflow() -> None:
    """Block-diagram schematic of the pipeline.

    Five stages in a single row. Real renders for the paper should be
    composed in Illustrator/Figma; this matplotlib version is a
    publication-ready fallback that compiles in CI.
    """
    apply_theme()
    fig, ax = plt.subplots(figsize=(7.2, 1.8))
    stages = [
        ("Target\nRBX1 (108 aa)", "RBX1"),
        ("Community\ndesign", "198 teams\n12,707 designs"),
        ("Organiser\nselection", "322 selected"),
        ("Expression\n+ SPR", "255 expressed\n9 binders"),
        ("Top binder", "26 nM\n(ORBIT)"),
    ]
    n = len(stages)
    box_w = 1.2
    gap = 0.45
    total = n * box_w + (n - 1) * gap
    x0 = -total / 2 + box_w / 2
    for i, (title, sub) in enumerate(stages):
        cx = x0 + i * (box_w + gap)
        rect = plt.Rectangle((cx - box_w / 2, -0.55), box_w, 1.1,
                              facecolor=CYAN if i == n - 1 else "white",
                              edgecolor=INK, linewidth=0.7)
        ax.add_patch(rect)
        ax.text(cx, 0.2, title, ha="center", va="center",
                fontsize=8, fontweight="bold")
        ax.text(cx, -0.25, sub, ha="center", va="center",
                fontsize=7, color=INK if i == n - 1 else SLATE)
        if i < n - 1:
            arrow_x0 = cx + box_w / 2
            arrow_x1 = cx + box_w / 2 + gap
            ax.annotate("", xy=(arrow_x1, 0), xytext=(arrow_x0, 0),
                        arrowprops=dict(arrowstyle="->", lw=0.6, color=INK))
    ax.set_xlim(-total / 2 - 0.3, total / 2 + 0.3)
    ax.set_ylim(-1, 1)
    ax.set_axis_off()
    save_figure(fig, "paper/fig1_workflow")
    plt.close(fig)


# ----------------------------------------------------------------------
# Figure 2 — 4-panel results landscape
# ----------------------------------------------------------------------
def fig2_results() -> None:
    apply_theme()
    df = load_designs()
    is_binder = df["is_binder"].astype("boolean").fillna(False)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))

    # A: counts per method family (top 10)
    top10 = (
        df["method_family"].value_counts().head(10).iloc[::-1]
    )
    ax = axes[0, 0]
    ax.barh(top10.index, top10.values, color=CYAN, edgecolor=INK)
    ax.set_xlabel("Designs")
    ax.set_title("A  Top 10 method families", loc="left", fontsize=9)

    # B: hit rate per method family (same top-10 set as panel A)
    ax = axes[0, 1]
    by_method = (
        df.assign(is_binder=is_binder)
        .groupby("method_family")
        .agg(n=("design_id", "count"), n_binders=("is_binder", "sum"))
        .loc[top10.index[::-1]]
    )
    by_method["hit_rate"] = by_method["n_binders"] / by_method["n"]
    ax.barh(by_method.index[::-1], by_method["hit_rate"].iloc[::-1],
            color=CYAN, edgecolor=INK)
    for i, (_, row) in enumerate(by_method[::-1].iterrows()):
        ax.text(row["hit_rate"] + 0.005, i,
                f"{int(row['n_binders'])}/{int(row['n'])}",
                va="center", fontsize=6, color=INK)
    ax.set_xlabel("Hit rate")
    ax.set_xlim(0, max(0.05, by_method["hit_rate"].max() * 1.4))
    ax.set_title("B  Hit rate by method family (top 10 by n)", loc="left", fontsize=9)

    # C: KD ladder for the 9 binders
    ax = axes[1, 0]
    binders = df[is_binder].sort_values("kd_nM_mean").copy()
    y = np.arange(len(binders))
    colors = [PAL.get(s, CYAN) for s in binders["binding_strength"]]
    ax.barh(y, binders["kd_nM_mean"],
            color=colors, edgecolor=INK)
    if "kd_nM_min" in binders.columns and "kd_nM_max" in binders.columns:
        ax.errorbar(binders["kd_nM_mean"], y,
                    xerr=[binders["kd_nM_mean"] - binders["kd_nM_min"],
                          binders["kd_nM_max"] - binders["kd_nM_mean"]],
                    fmt="none", ecolor=INK, capsize=2, lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(binders["pb_id"], fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("KD (nM)")
    ax.set_title(f"C  KD for the {len(binders)} binders", loc="left", fontsize=9)

    # D: AFDB50 identity vs KD (binders only with identity)
    ax = axes[1, 1]
    plot_df = binders[binders["pb_seqidentity_afdb50"].notna()].copy()
    if len(plot_df):
        ax.scatter(plot_df["pb_seqidentity_afdb50"], plot_df["kd_nM_mean"],
                   color=CYAN, edgecolor=INK, s=40, zorder=3)
        # Stagger label offsets so the tight cluster around 150-300 nM doesn't
        # collide. Sort by KD so adjacent labels alternate above/below.
        plot_df = plot_df.sort_values("kd_nM_mean").reset_index(drop=True)
        for i, r in plot_df.iterrows():
            dy = 8 if i % 2 == 0 else -10
            dx = 6 if i % 2 == 0 else 6
            ax.annotate(r["pb_id"], (r["pb_seqidentity_afdb50"], r["kd_nM_mean"]),
                        fontsize=5, color=INK, alpha=0.85,
                        xytext=(dx, dy), textcoords="offset points")
    ax.set_xlabel("Sequence identity vs AFDB50 (%)")
    ax.set_ylabel("KD (nM)")
    ax.set_yscale("log")
    ax.set_title("D  Novelty vs affinity", loc="left", fontsize=9)

    fig.tight_layout()
    save_figure(fig, "paper/fig2_results")
    plt.close(fig)


# ----------------------------------------------------------------------
# Supplementary S1 — modality × method-family crosstab
# ----------------------------------------------------------------------
def figS1_modality_method() -> None:
    apply_theme()
    df = load_designs()
    xtab = (
        df.assign(
            modality=df["pb_design_class"].fillna("Unknown"),
            method=df["method_family"].fillna("Not mentioned"),
        )
        .groupby(["method", "modality"])
        .size()
        .unstack(fill_value=0)
    )
    # keep the top 15 methods by total to fit
    xtab = xtab.loc[xtab.sum(axis=1).sort_values(ascending=False).head(15).index]
    xtab = xtab.loc[:, xtab.sum(axis=0).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    im = ax.imshow(xtab.values, cmap="cividis", aspect="auto")
    ax.set_xticks(range(len(xtab.columns)))
    ax.set_xticklabels(xtab.columns, rotation=35, ha="right", fontsize=7)
    ax.set_yticks(range(len(xtab.index)))
    ax.set_yticklabels(xtab.index, fontsize=7)
    for i in range(len(xtab.index)):
        for j in range(len(xtab.columns)):
            v = xtab.values[i, j]
            if v:
                ax.text(j, i, str(v), ha="center", va="center",
                        fontsize=6,
                        color="white" if v < xtab.values.max() * 0.55 else INK)
    fig.colorbar(im, ax=ax, label="designs", fraction=0.04, pad=0.02)
    ax.set_title("Modality × method family (top 15 methods)", fontsize=9)
    fig.tight_layout()
    save_figure(fig, "paper/figS1_modality_method_crosstab")
    plt.close(fig)


# ----------------------------------------------------------------------
# Supplementary S2 — declared target region vs binding outcome
# ----------------------------------------------------------------------
def figS2_target_region() -> None:
    apply_theme()
    df = load_designs()
    df = df[df["submission_target_region"].notna()].copy()
    df["outcome"] = np.where(df["is_binder"].fillna(False), "Binder", "Non-binder")
    xtab = (
        df.groupby(["submission_target_region", "outcome"])
        .size()
        .unstack(fill_value=0)
        .sort_values("Binder" if "Binder" in df["outcome"].unique() else "Non-binder",
                     ascending=False)
    )
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    xs = np.arange(len(xtab))
    bottom = np.zeros(len(xtab))
    for outcome, colour in (("Non-binder", SLATE), ("Binder", CYAN)):
        if outcome in xtab.columns:
            ax.bar(xs, xtab[outcome], bottom=bottom, label=outcome,
                   color=colour, edgecolor=INK, width=0.8)
            # annotate binder counts
            if outcome == "Binder":
                for x, v in zip(xs, xtab[outcome]):
                    if v > 0:
                        ax.text(x, bottom[x] + v + 1, str(int(v)),
                                ha="center", fontsize=7, color=INK)
            bottom += xtab[outcome].values
    ax.set_xticks(xs)
    ax.set_xticklabels(xtab.index, rotation=25, ha="right", fontsize=7)
    ax.set_ylabel("Designs")
    ax.set_title("Submitter-declared target region vs outcome", fontsize=9)
    ax.legend(loc="upper right", frameon=False, fontsize=7)
    fig.tight_layout()
    save_figure(fig, "paper/figS2_target_region")
    plt.close(fig)


# ----------------------------------------------------------------------
# Supplementary S3 — KD with per-replicate spread
# ----------------------------------------------------------------------
def figS3_kd_replicates() -> None:
    apply_theme()
    binders = load_designs(only_binders=True).sort_values("kd_nM_mean").copy()
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(binders))
    colors = [PAL.get(s, CYAN) for s in binders["binding_strength"]]
    ax.barh(y, binders["kd_nM_mean"],
            color=colors, edgecolor=INK)
    if "kd_nM_min" in binders.columns and "kd_nM_max" in binders.columns:
        ax.errorbar(binders["kd_nM_mean"], y,
                    xerr=[binders["kd_nM_mean"] - binders["kd_nM_min"],
                          binders["kd_nM_max"] - binders["kd_nM_mean"]],
                    fmt="none", ecolor=INK, capsize=2, lw=0.5)
        for yi, (_, r) in enumerate(binders.iterrows()):
            label = f"{r['kd_nM_mean']:.0f} nM  [{r['kd_nM_min']:.0f}, {r['kd_nM_max']:.0f}]"
            ax.text(r["kd_nM_max"] * 1.05, yi, label, va="center", fontsize=6, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(binders["pb_id"], fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("KD (nM)")
    ax.set_title("Per-replicate KD spread for the 9 binders", fontsize=9)
    fig.tight_layout()
    save_figure(fig, "paper/figS3_kd_replicates")
    plt.close(fig)


# ----------------------------------------------------------------------
# Supplementary S4 — multi-model complex re-scoring (requires grand_metrics.csv)
# ----------------------------------------------------------------------
def figS4_complex_grand_metrics() -> None:
    apply_theme()
    root = repo_root()
    grand = root / "data" / "grand_metrics.csv"
    if not grand.exists():
        print(f"[figS4] skipped — {grand} missing. Run `mise run build:grand` after all 3 models download.")
        return
    df = pd.read_csv(grand)
    # Look for canonical ipSAE columns. The exact column names depend on
    # build_grand_metrics.py's schema — pick the first that exists per model.
    candidates = {
        "boltz2":   ["b2_ipsae_d0chn_max", "b2_ipsae", "b2_iptm"],
        "chai":     ["chai_ipsae_d0chn_max", "chai_aggregate_score", "chai_iptm"],
        "protenix": ["px_ipsae_d0chn_max", "px_iptm", "px_ranking_score"],
    }
    cols = {m: next((c for c in cs if c in df.columns), None) for m, cs in candidates.items()}
    if not all(cols.values()):
        print(f"[figS4] skipped — missing model columns in grand_metrics.csv: {cols}")
        return
    is_binder = df.get("is_binder", pd.Series([False] * len(df))).astype(str).str.lower().eq("true")

    # Skip pairs where the y-axis model has all-null values (e.g. Protenix
    # rerun still in flight). Drop down to one panel rather than show empty axes.
    candidate_pairs = [("boltz2", "chai"), ("boltz2", "protenix")]
    pairs = [
        (xm, ym) for xm, ym in candidate_pairs
        if pd.to_numeric(df[cols[ym]], errors="coerce").notna().any()
    ]
    if not pairs:
        print("[figS4] skipped — every model column is null in grand_metrics.csv.")
        return

    fig, axes = plt.subplots(1, len(pairs), figsize=(3.5 * len(pairs) + 0.2, 3.4), squeeze=False)
    for ax, (xm, ym) in zip(axes[0], pairs):
        x = df[cols[xm]]; y = df[cols[ym]]
        ax.scatter(x[~is_binder], y[~is_binder], color=SLATE, edgecolor=INK,
                   s=14, alpha=0.5, label="Non-binder")
        ax.scatter(x[is_binder], y[is_binder], color=CYAN, edgecolor=INK,
                   s=40, zorder=3, label="Binder")
        ax.set_xlabel(f"{xm}  {cols[xm]}")
        ax.set_ylabel(f"{ym}  {cols[ym]}")
        ax.legend(loc="lower right", fontsize=7, frameon=False)
    fig.suptitle("Complex re-scoring across models", fontsize=9)
    fig.tight_layout()
    save_figure(fig, "paper/figS4_complex_grand_metrics")
    plt.close(fig)


def main() -> None:
    out = repo_root() / "figures" / "paper"
    out.mkdir(parents=True, exist_ok=True)
    fig1_workflow();              print("[figures] wrote fig1_workflow")
    fig2_results();               print("[figures] wrote fig2_results")
    figS1_modality_method();      print("[figures] wrote figS1_modality_method_crosstab")
    figS2_target_region();        print("[figures] wrote figS2_target_region")
    figS3_kd_replicates();        print("[figures] wrote figS3_kd_replicates")
    figS4_complex_grand_metrics();print("[figures] wrote figS4_complex_grand_metrics")


if __name__ == "__main__":
    main()
