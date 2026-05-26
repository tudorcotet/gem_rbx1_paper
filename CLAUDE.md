# Claude Code instructions — rbx1_gem_paper

## Read first

[`docs/DATA.md`](docs/DATA.md) for the canonical 322-row table. Then
this file for the collaboration rules.

## What this repo is

Data, code, figures, and manuscript source for the **GEM × Adaptyv
2026 RBX1 binder design competition**. 322 designs across 48 teams
through Adaptyv SPR. Outcome: **9 binders, 1 Strong (26 nM)**, 79%
expression rate.

The paper is a brief communication, two figures max. Tone is
descriptive — "we observe X" — not significance claims on top of N=322
with a 2.8% hit rate.

## Contracts

1. **`data/designs.csv` is the source of truth.** 322 rows.
   `design_id` is the join key everywhere. Load via
   `from scripts.utils import load_designs`.
2. **Raw inputs are immutable.** `data/raw/proteinbase/...csv` and
   `data/raw/submissions/...csv` are not edited. Fix `build_designs.py`
   and re-run `mise run build`.
3. **Analyses write to their own subdir only.** Copy
   `analyses/_template/` to start. Outputs (report.md, summary.json,
   derived csvs) are gitignored — they regenerate from `main.py`.
4. **Plotting:** `from scripts.plotting import apply_theme, save_figure`.
   Colours from `theme/palettes.json`. Save into `figures/paper/`.

## Adding an analysis

```bash
cp -r analyses/_template analyses/<name>
$EDITOR analyses/<name>/main.py
$EDITOR mise.toml                # add [tasks."analysis:<name>"]
mise run analysis:<name>
```

`main.py` must:

1. Load via `load_designs()` — never a hard-coded path.
2. Write only into `analyses/<name>/`.
3. Be idempotent.

## Plotting conventions

- **Palette:** `theme/palettes.json`. Binding strength is
  `pal["binding_strength_2026"]`. Same hues across every figure.
- **Fonts:** `apply_theme()` falls back Geist → Helvetica → DejaVu Sans.
  Don't hard-code a font.
- **Frame:** white canvas, 0.5pt spines, bottom + left only.
- **No 3D, no rainbow, no gridlines on top of data.**
- **Stat annotations always pair p-value with effect size.** With 9
  binders, most per-method splits are underpowered — report raw n/N
  instead.
- **Sizes:** single column 3.5×2.6 in, 1.5-column 5.0×3.4 in, full page
  7.2×4.8 in.
- **File naming:** `figures/paper/fig<N>_<short_name>.{png,pdf,svg}`.

## Voice (paper)

Clinical, methods-paper voice. Short sentences. Numbers where you'd
reach for an adjective.

- Cite the column name: "Hit rate (`is_binder`) was 2.8% (9/322)."
- State the test, pair p with effect size.
- Caveat upfront, not in the discussion.
- Avoid: *robust, leverage, delve, comprehensive, novel insights,
  stark, striking, remarkable*; the rule of three; "It's not just X —
  it's Y." parallelism.
- Use "we observe X" not "X is significant".

Blog and social copy live elsewhere — don't blend the two voices.

## Numbers come from disk, not from the model

- KDs come from `kd_nM_mean` / `pkd_arith_mean`. Cite the column.
- ESMFold / ProteinMPNN come from `pb_*` columns. Null-check before use.
- Team labels come from `team`. Method from `method_family`. Modality
  from `pb_design_class` (classifier) or `submission_modality`
  (designer-declared).
- The competition target stance is in `submission_target_region` and
  `submission_targets_idr`. 7 of 9 binders target the IDR/N-terminus.

## Domain notes

- **Target:** human RBX1 (UniProt **P62877**, 108 aa). Residues 12–39
  are disordered. Residues 40–108 are the RING-H2 finger with three
  structural Zn²⁺.
- **Default-buffer screen only.** The Zn²⁺ rerun is in flight.
- **BindCraft 2** is Pacesa Lab's unpublished branch. Cite as personal
  communication.
- **ORBIT** (Mandrake Bio) has a blog post but no preprint — the
  Strong binder. Acknowledge by name.
- The rest of the field goes into a consortium author block. Get
  written permission before naming a participant.

## What NOT to do

- Don't re-implement `load_designs()`. Import it.
- Don't read `data/raw/*` directly from an analysis — `designs.csv`
  has the canonical join.
- Don't edit `data/designs.csv` by hand. Fix the build script.
- Don't claim significance on per-method hit rates.

