# `theme/` — the brand visual system

Adaptyv brand palette, fonts, and matplotlibrc for every figure in
the paper. Pin once here so every analysis renders against the same
visual system.

## Files

| file              | role |
|---                |---|
| `palettes.json`   | Canonical colour spec. Every figure colour resolves here. |
| `matplotlibrc`    | Publication-grade rcParams (fonts, spines, ticks, legend, savefig). |
| `mpl_theme.py`    | Python helpers: `set_brand_style()`, `BRAND_COLORS`, `get_brand_palettes()`. |

## How to use

The fastest path. Two lines, one figure:

```python
from theme.mpl_theme import set_brand_style
from scripts.plotting import save_figure
import matplotlib.pyplot as plt

set_brand_style()                       # apply rcParams + brand fonts

fig, ax = plt.subplots(figsize=(3.5, 2.6))
ax.bar(["RFdiffusion", "BoltzGen", "Mosaic"], [27, 22, 25])
ax.set_ylabel("Designs")
save_figure(fig, "paper/fig_demo")        # writes png + pdf + svg
```

The matching one-liner via the plotting helpers:

```python
from scripts.plotting import apply_theme, save_figure

apply_theme()                             # calls set_brand_style + loads rc
# ...
save_figure(fig, "paper/figN_short_name")
```

## Picking colours

Always go through the named palettes — no ad hoc hex codes.

```python
from theme.mpl_theme import BRAND_COLORS, get_brand_palettes

BRAND_COLORS["cyan"]            # "#30C5F5" — brand primary
BRAND_COLORS["good"]            # "#1FE48F" — accent green
BRAND_COLORS["non_binder"]      # "#5C6773" — muted slate
BRAND_COLORS["missing_data"]    # "#3E6175"

pal = get_brand_palettes()
pal["binding_strength_2026"]      # {"Binder": "#00D9FF", "Strong": "#30C5F5", ...}
pal["method_group"]               # method-family palette ported from EGFR
pal["molecule_type"]              # miniprotein / nanobody / peptide / scFv
```

For the RBX1 paper, the load-bearing palette is
`binding_strength_2026`: cyan for binders (graded by strength), slate
for non-binders, muted grey for no expression.

## Fonts

`set_brand_style()` resolves the body font as **Geist → Roboto →
DejaVu Sans** (whichever is installed first). The publication
`matplotlibrc` widens the chain to **Helvetica → Arial → Liberation
Sans → DejaVu Sans** for journal compatibility.

To install Geist locally:

1. Download from `https://vercel.com/font` (Geist is OFL-licensed).
2. Drop the TTFs into `~/Library/Fonts/` (macOS) or
   `~/.local/share/fonts/` (Linux).
3. Restart your Python kernel so `matplotlib.font_manager` rebuilds its
   cache.

The display face (GT Pressura Extended) is a paid licence and is
**not** embedded in matplotlib output. It's only used in HTML / SVG
hero assets where the font can be loaded via FontFace. The paper
figures stay Geist / Helvetica.

## Verify the theme renders

```bash
mise run figures
```

Re-renders every figure in `figures/paper/`. If colours or fonts look
wrong, your local font install is likely the reason; see the Fonts
section above.
