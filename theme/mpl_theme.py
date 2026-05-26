# pyright: reportMissingTypeStubs=false
"""
Brand matplotlib theme — bioarena-trem2 (May 2026).

This is the canonical Python plotting style for this repo. Palettes are
loaded from `theme/palettes.json` so figures stay 1:1 aligned with the
ggplot R theme (`theme/ggplot_theme.R`) and the CSS tokens
(`theme/tokens.css`).

Public API
----------
``set_brand_style()``
    Apply rcParams. Call once at the top of any analysis script.

``BRAND_COLORS``
    A flat ``Dict[str, str]`` of every named brand color
    (``"cyan"``, ``"cyan_soft"``, ``"good"``, ``"ink"``, …) — the
    fastest way to grab a hex code without nested ``palettes[k][v]``
    indexing.

``get_brand_palettes()``
    Full categorical-palette dict (``binding_strength``,
    ``design_method``, etc.). Each value is a ``Dict[label, hex]``
    suitable for ``sns.color_palette`` or ``ax.bar(... color=...)``.

``apply_brand_blog_post_theme(fig, ax, ...)``
    One-shot styling helper for blog-style figures (title, subtitle,
    legend on the right, no grid, hairline axes).

Imports
-------
Used throughout the repo as::

    from theme.mpl_theme import set_brand_style, BRAND_COLORS

Drift note
----------
matplotlib cannot embed GT Pressura Extended (paid display face) inside
PNGs without bundling it into every render — so the figure font stays
Geist > Roboto > DejaVu Sans. GT Pressura is reserved for HTML/SVG hero
assets where the font can be loaded via FontFace.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, TypedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


class _ThemeSpec(TypedDict):
    style: Dict[str, object]
    core: Dict[str, str]
    palettes: Dict[str, Dict[str, str]]


# ────────────────────────────────────────────────────────────────────
# Spec loading
# ────────────────────────────────────────────────────────────────────


def _spec_path() -> Path:
    """Resolve `theme/palettes.json` relative to this file (no repo-root walk).

    Older versions walked up looking for ``pyproject.toml``; that breaks
    when the repo is unpacked outside a Python project. The theme dir is
    self-contained — its sibling JSON is always next to this file.
    """
    return Path(__file__).resolve().parent / "palettes.json"


def _load_spec() -> _ThemeSpec:
    with _spec_path().open("r", encoding="utf-8") as f:
        spec = json.load(f)
    # Strip the leading `_meta` key — it's documentation, not data.
    spec.pop("_meta", None)
    return spec  # type: ignore[return-value]


_SPEC = _load_spec()


# ────────────────────────────────────────────────────────────────────
# Public flat-color export
# ────────────────────────────────────────────────────────────────────


def _build_flat_colors() -> Dict[str, str]:
    """Flatten the `core` namespace + selected canonical names from palettes."""
    flat: Dict[str, str] = dict(_SPEC.get("core", {}))
    # Promote the canonical 2026 binding-strength colors as `binder_*` keys
    bs = _SPEC["palettes"].get("binding_strength_2026", {})
    flat.update({
        "binder":           bs.get("Binder",          flat.get("cyan_binder", "#00D9FF")),
        "binder_strong":    bs.get("Strong",          flat.get("cyan",        "#30C5F5")),
        "binder_medium":    bs.get("Medium",          flat.get("cyan_soft",   "#33C4FF")),
        "binder_weak":      bs.get("Weak",            flat.get("cyan_deep",   "#36B7F6")),
        "non_binder":       bs.get("Non-binder",      "#5C6773"),
        "missing_data":     bs.get("Missing data",    "#3E6175"),
    })
    return flat


BRAND_COLORS: Dict[str, str] = _build_flat_colors()
"""Flat dict of every named brand color (cyan/ink/good/binder/etc.)."""


# ────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────


def _as_str(value: object, default: str) -> str:
    try:
        if value is None:
            return default
        s = str(value)
        return s if s else default
    except Exception:
        return default


def _as_float(value: object, default: float) -> float:
    try:
        if value is None:
            return default
        return float(value)  # type: ignore[arg-type]
    except Exception:
        return default


def _resolve_font_family() -> str:
    """Return the first matplotlib-installed brand font, falling back gracefully.

    Order: Geist (canonical) -> Roboto (legacy 2024 EGFR theme fallback) ->
    DejaVu Sans (always present in matplotlib).
    """
    try:
        from matplotlib import font_manager
        available = {f.name for f in font_manager.fontManager.ttflist}
        for candidate in ("Geist", "Roboto"):
            if candidate in available:
                return candidate
    except Exception:
        pass
    return "DejaVu Sans"


# ────────────────────────────────────────────────────────────────────
# Palette accessors
# ────────────────────────────────────────────────────────────────────


def get_brand_palettes() -> Dict[str, Dict[str, str]]:
    """Return all categorical palettes, with an aggregated ``_all`` fallback.

    The ``_all`` palette is every unique color in insertion order — useful
    as a default cycler when a series doesn't map to a named category.
    """
    palettes: Dict[str, Dict[str, str]] = {
        k: dict(v) for k, v in _SPEC["palettes"].items()
    }
    all_colors: List[str] = []
    for pal in palettes.values():
        all_colors.extend(list(pal.values()))
    unique_colors = list(dict.fromkeys(all_colors))
    palettes["_all"] = {str(i): c for i, c in enumerate(unique_colors)}
    return palettes


# ────────────────────────────────────────────────────────────────────
# Theme application
# ────────────────────────────────────────────────────────────────────


def set_brand_style(*, dark: bool = False) -> None:
    """Apply matplotlib rcParams to match the canonical theme.

    Parameters
    ----------
    dark
        When ``True``, swap the canvas to ``--bg`` (#0F1419) with white
        ink — useful for screen-only figures meant to sit in a dark
        webpage. Defaults to the print-safe light theme.
    """
    style = _SPEC["style"]
    font_family = _resolve_font_family()

    if dark:
        text_color = _as_str(style.get("text_color_dark_canvas"), "#FFFFFF")
        background_color = _as_str(style.get("background_color_dark"), "#0F1419")
        axis_line_color = "#FFFFFF"
        legend_border_color = "#FFFFFF"
    else:
        text_color = _as_str(style.get("text_color"), "#0F1419")
        background_color = _as_str(style.get("background_color"), "white")
        axis_line_color = _as_str(style.get("axis_line_color"), "#0F1419")
        legend_border_color = _as_str(style.get("legend_border_color"), "#0F1419")

    axis_line_width = _as_float(style.get("axis_line_width"), 0.5)

    mpl.rcParams.update({
        # Base font / colors
        "font.family": font_family,
        "text.color": text_color,
        "axes.labelcolor": text_color,
        "xtick.color": text_color,
        "ytick.color": text_color,
        # Backgrounds
        "figure.facecolor": background_color,
        "axes.facecolor": background_color,
        # Grid (off in canonical brand style — let the data carry the surface)
        "axes.grid": False,
        # Spines
        "axes.edgecolor": axis_line_color,
        "axes.linewidth": axis_line_width,
        # Legend
        "legend.frameon": True,
        "legend.edgecolor": legend_border_color,
        "legend.framealpha": 1.0,
        "legend.facecolor": background_color,
        # Default sizes (per-figure overrides welcome)
        "axes.titlesize": 20,
        "axes.titleweight": "bold",
        "axes.labelsize": 12,
        "axes.labelweight": "bold",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })


# Backwards-compatible alias used by older scripts ported from egfr2024.
set_brand_matplotlib_theme = set_brand_style


def apply_brand_blog_post_theme(
    fig: Figure,
    ax: Axes,
    *,
    title: str,
    subtitle: str,
    x_label: Optional[str] = None,
    y_label: Optional[str] = "Number of designs",
    legend_title: Optional[str] = None,
) -> None:
    """Apply the standard blog-post styling to a single Axes.

    Lays out a bold title, a cyan-deep subtitle, and a right-side legend.
    Leaves data plotting to the caller — this only configures chrome.
    """
    set_brand_style()
    style = _SPEC["style"]
    subtitle_color = _as_str(style.get("subtitle_color"), "#36B7F6")
    legend_border_color = _as_str(style.get("legend_border_color"), "#0F1419")
    background_color = _as_str(style.get("background_color"), "white")

    if x_label is not None:
        ax.set_xlabel(x_label)
    if y_label is not None:
        ax.set_ylabel(y_label)

    if title:
        fig.suptitle(title, fontsize=20, fontweight="bold", y=0.98)
    if subtitle:
        ax.set_title(str(subtitle), color=subtitle_color, fontsize=16, pad=14)

    if legend_title is not None:
        leg = ax.legend(
            title=legend_title, loc="center left",
            bbox_to_anchor=(1.02, 0.5), frameon=True,
        )
    else:
        leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)

    if leg is not None:
        leg.get_frame().set_edgecolor(legend_border_color)
        leg.get_frame().set_facecolor(background_color)

    plt.subplots_adjust(right=0.78, left=0.12, top=0.90, bottom=0.22)


# ────────────────────────────────────────────────────────────────────
# Label helpers (1:1 with the R theme)
# ────────────────────────────────────────────────────────────────────


_TOKEN_MAP: Dict[str, str] = {
    "af2": "AF2",
    "ml": "ML",
    "plm": "PLM",
    "esm": "ESM",
    "esm-if": "ESM-IF",
    "rfdiffusion": "RFdiffusion",
    "bindcraft": "BindCraft",
    "timed": "TIMED",
    "proteinmpnn": "ProteinMPNN",
    "trem2": "TREM2",
    "r1": "R1",
    "r2": "R2",
}


def _sentence_case(label: str) -> str:
    """Convert snake_case to sentence case without destroying acronyms."""
    s = " ".join(label.replace("_", " ").split())
    if not s:
        return s
    return s[0].upper() + s[1:]


def _missing_like(value: object) -> bool:
    if value is None:
        return True
    s = str(value).strip()
    if not s:
        return True
    return s.lower() in {"na", "nan", "none", "null"}


def labelize_category_value(value: object) -> str:
    """Humanize category values for display (snake_case -> sentence case)."""
    if _missing_like(value):
        return "Not mentioned"
    raw = str(value).strip()
    s = raw.replace("/", " / ").replace("_", " ")
    s = " ".join(s.split())
    words = []
    for w in s.split(" "):
        lw = w.lower()
        words.append(_TOKEN_MAP.get(lw, lw))
    s2 = " ".join(words)
    return s2[:1].upper() + s2[1:] if s2 else s2


def format_column_label(column: str) -> str:
    """Format a dataframe column name into a standardized axis/legend label."""
    col = str(column).strip()
    cl = col.lower()
    fixed: Dict[str, str] = {
        "selected": "Selection status",
        "binding_strength": "Binding strength",
        "design_category": "Design category",
        "design_strategy": "Design strategy",
        "design_method": "Design method",
        "binding_method": "Binding method",
        "solubility_method": "Solubility method",
        "method_paradigm": "Method paradigm",
        "molecule_type": "Molecule type",
        "parent_molecule": "Parent molecule",
    }
    if cl in fixed:
        return fixed[cl]
    return labelize_category_value(col)


def _is_kd_metric(metric: str) -> bool:
    mm = metric.strip().lower()
    base = mm.removeprefix("-log10_").removeprefix("normalized_")
    return base == "kd" or base.endswith("_kd")


def format_metric_name(metric: str, *, for_title: bool = False) -> str:
    """Format metric identifiers into standardized display labels.

    KD is rendered as ``$K_D$``; -log10(KD) is rendered as ``$pK_D$``.
    """
    m = metric.strip()
    ml = m.lower()

    if _is_kd_metric(ml):
        if ml.startswith("-log10_"):
            return r"Binding affinity ($pK_D$)" if for_title else r"$pK_D$"
        return r"Binding affinity ($K_D$)" if for_title else r"$K_D$"

    title_map: Dict[str, str] = {
        "iptm": "ipTM",
        "ptm": "pTM",
        "plddt": "pLDDT",
        "pae_interaction": "iPAE",
        "ipae": "iPAE",
        "esm_pll": "ESM PLL",
        "normalized_esm_pll": "Normalized ESM PLL",
        "sequence_length": "Sequence length",
        "interface_nres": "Interface residues",
    }

    if for_title:
        if ml in title_map:
            return title_map[ml]
        if ml.startswith("-log10_"):
            base = ml.removeprefix("-log10_")
            return title_map.get(base, _sentence_case(base))
        return title_map.get(ml, _sentence_case(ml))

    if ml in title_map:
        return title_map[ml]
    if ml.startswith("-log10_"):
        base = ml.removeprefix("-log10_")
        return f"-log10({title_map.get(base, _sentence_case(base))})"
    return title_map.get(ml, _sentence_case(ml))


__all__ = [
    "BRAND_COLORS",
    "set_brand_style",
    "set_brand_matplotlib_theme",  # legacy alias
    "apply_brand_blog_post_theme",
    "get_brand_palettes",
    "format_metric_name",
    "format_column_label",
    "labelize_category_value",
]
