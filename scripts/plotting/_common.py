"""Plotting helpers shared by every figure script.

Apply the project theme once at the top of any plotting script::

    from scripts.plotting import apply_theme, save_figure

    apply_theme()
    fig, ax = plt.subplots()
    ...
    save_figure(fig, "paper/fig1_hit_rate")
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

from scripts.utils.load_data import repo_root


def apply_theme() -> None:
    """Load `theme/matplotlibrc` and switch matplotlib to those defaults.

    No-op if the rc file is missing — useful when collaborators are running in
    a stripped-down environment without our fonts installed.
    """
    rc = repo_root() / "theme" / "matplotlibrc"
    if rc.exists():
        mpl.rc_file(rc, use_default_template=False)
    # Also call into the canonical brand theme module if it's importable.
    try:
        from theme.mpl_theme import set_brand_style

        set_brand_style()
    except Exception:
        pass


def save_figure(
    fig: plt.Figure,
    stem: str | Path,
    *,
    formats: tuple[str, ...] = ("png", "pdf", "svg"),
    dpi: int = 300,
) -> list[Path]:
    """Save `fig` to `figures/<stem>.{png,pdf,svg}` under the repo root.

    `stem` may be relative (e.g. `"paper/fig1_hit_rate"`) or absolute. Returns
    the list of files actually written.
    """
    stem_path = Path(stem)
    if not stem_path.is_absolute():
        stem_path = repo_root() / "figures" / stem_path
    stem_path.parent.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for ext in formats:
        out = stem_path.with_suffix(f".{ext}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        written.append(out)
    return written
