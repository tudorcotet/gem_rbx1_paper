# RBX1 — the GEM × Adaptyv 2026 paper

322 designs against **RBX1** (the RING-H2 catalytic core of every
cullin-RING ligase). 48 teams, methods spanning RFdiffusion, BoltzGen,
BindCraft 2, hallucination, MoPPIt, plain Rosetta. One Adaptyv SPR
pipeline. **9 binders, 1 Strong (KD = 26 nM)**.

The competition shipped April 26 2026 at the ICLR GEM workshop. This
repo holds the data, the analyses, and the manuscript source for the
community write-up.

## The cohort

|                | designs | expressed | binders | strongest |
|---             |---      |---        |---      |---        |
| Public release | **322** | 255 (79%) | **9**   | 26 nM     |

`data/designs.csv` is the canonical table. Column reference is in
[`docs/DATA.md`](docs/DATA.md). Collaboration rules are in
[`CLAUDE.md`](CLAUDE.md).

## Quick start

```bash
mise run setup              # uv sync
mise run build              # raw csvs → data/designs.csv
mise run mirror:structures  # ESMFold CIFs, PNGs, SPR sensorgrams → data/*
mise run figures            # render every figure paper/main.tex references
mise run analysis:all       # re-run every canonical analysis
```

```python
from scripts.utils import load_designs
df = load_designs()                       # 322 rows
df = load_designs(only_binders=True)      #   9 rows
df = load_designs(only_expressed=True)    # 255 rows
```

No GPU, no credentials. The canonical CSV is self-contained.

## Layout

```
rbx1_gem_paper/
├── data/                       canonical CSV + immutable raw inputs
│   ├── designs.{csv,parquet,fasta}      ⭐ load this
│   ├── target/rbx1.fasta                UniProt P62877, 108 aa
│   └── raw/
│       ├── proteinbase/...csv           ProteinBase export — DO NOT EDIT
│       └── submissions/...csv           per-submission form export (emails redacted)
├── docs/DATA.md                every column of designs.csv
├── analyses/                   one subdir per analysis (yours goes here)
│   ├── _template/              copy this to start
│   ├── overview/               counts, hit rates, modality × method
│   └── leaderboard/            top binders, KD table
├── scripts/
│   ├── utils/load_data.py      the canonical loader
│   ├── plotting/_common.py     apply_theme + save_figure
│   └── data/
│       ├── build_designs.py    raw csvs → designs.csv
│       └── mirror_structures.py  mirror CIFs / PNGs / SPR sensorgrams
├── theme/                      Adaptyv brand kit (palette, matplotlibrc)
├── figures/{paper,exploration} rendered figures (exploration is gitignored)
└── paper/                      main.tex + references.bib + sections/*.md
```

## Adding your analysis

```bash
cp -r analyses/_template analyses/yourname
$EDITOR analyses/yourname/main.py
# add a [tasks."analysis:yourname"] block in mise.toml
mise run analysis:yourname
```

Full conventions in [`CLAUDE.md`](CLAUDE.md).

## Known caveats

- **Public release ≠ everything tested.** 322 designs hit ProteinBase;
  the wet-lab batch ran ~422 (selection cap + scratched designs).
- **Hit rate is thin** (9 / 322 = 2.8%). Per-method or per-team rates
  collapse to anecdotes. Per-cohort claims are load-bearing.
- **Modality is imbalanced.** ProteinBase tags 51% as `Other`; submitter
  declarations are richer (`submission_modality`).
- **Zinc-buffer rerun is in flight.** Default-buffer results may
  underestimate binding to the holo state.
- **N-terminal disordered tail is the most-claimed binding region.**
  7 of 9 binders explicitly target it (`submission_target_region`).
  Treat as a novelty AND a specificity concern until the Zn rerun and
  the ovalbumin panel land.

## Citation

```bibtex
@misc{rbx1_2026,
  author    = {{The GEM $\times$ Adaptyv RBX1 Consortium}},
  title     = {A community benchmark of de novo binder design against {RBX1}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {TBD},
  url       = {https://github.com/tudorcotet/gem_rbx1_paper},
}
```

## License

Code is MIT. Data and figures are CC-BY-4.0. The upstream ProteinBase
release is ODC-ODbL — derivatives inherit. See [`LICENSE`](LICENSE).
