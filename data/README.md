# `data/`

One canonical CSV, two immutable raw inputs, one target FASTA, plus
binary artifacts (structures, sensorgrams, renders, per-model metric
JSONs) mirrored locally through Git LFS.

```
data/
├── designs.{csv,parquet,fasta}     ⭐ canonical 322-row table
├── target/rbx1.fasta                UniProt P62877, 108 aa
├── raw/
│   ├── proteinbase/                 ProteinBase per-design export
│   └── submissions/                 per-submission form (emails redacted)
├── structures/                      predicted structures — one subdir per model
│   ├── esmfold/<pb_id>.cif          monomer fold (ProteinTyper / ESMFold)
│   ├── boltz2/<pb_id>.cif           complex prediction (rescoring stack)
│   ├── protenix/<pb_id>.cif         complex prediction (Protenix v2)
│   └── chai/<pb_id>.cif             complex prediction (Chai)
├── metrics/                         raw per-model output JSONs
│   ├── proteintyper/<pb_id>.json    full TyperJobOutput
│   ├── boltz2/<pb_id>.json          ipSAE / ipTM / pLDDT / PAE matrix
│   ├── protenix/<pb_id>.json        Protenix v2 raw scores
│   └── chai/<pb_id>.json            Chai raw scores
├── images/<pb_id>.png               stylised renders (one per design)
└── sensorgrams/                     SPR kinetic traces (943 spr + 3 bli fallback)
    └── <pb_id>_rep<NN>_{spr,bli}.json
```

Convention: every artifact key is the **ProteinBase slug** (`pb_id`).
Path discovery has no lookup column — `data/structures/esmfold/<pb_id>.cif`
is canonical. Use the helpers in `scripts.utils`:

```python
from scripts.utils import structure_path, image_path, sensorgram_paths, metrics_path

structure_path("small-vole-maple")              # → data/structures/esmfold/...
structure_path("small-vole-maple", model="boltz2")
image_path("small-vole-maple")
sensorgram_paths("small-vole-maple")            # list of replicate JSONs
metrics_path("small-vole-maple", model="proteintyper")
```

All four return `None` (or `[]`) when the artifact hasn't been
mirrored yet. Run `mise run mirror:structures` to populate from
ProteinBase URLs (CIFs, PNGs, and SPR sensorgrams all in one pass).

## Load

```python
from scripts.utils import load_designs
df = load_designs()                       # 322 rows
df = load_designs(only_binders=True)      #   9 rows
df = load_designs(only_expressed=True)    # 255 rows
```

Column reference: [`../docs/DATA.md`](../docs/DATA.md).

## Regenerate

```bash
mise run build              # raw csvs → designs.csv + designs.parquet + designs.fasta
mise run mirror:structures  # ESMFold CIFs + PNGs + sensorgrams → data/{structures,images,sensorgrams}
mise run rerun:typer        # re-run ProteinTyper on the 4 designs without predictions
mise run rerun:complex      # Boltz-2 + Chai-1 + Protenix-v2 complex predictions (Modal --detach)
mise run build:grand        # pool every per-model metric into data/grand_metrics.csv
```

`scripts/data/build_designs.py` is the only writer of `designs.csv`.
Deterministic. Re-running produces byte-identical output. If a row is
wrong, fix the build script.

## Provenance

| file                                                  | source                                                                                                              |
|---                                                    |---                                                                                                                  |
| `raw/proteinbase/...csv`                              | `https://proteinbase.com/api/proteins/download?collectionId=03ec16ff-…&slug=gem-x-adaptyv-rbx1-…`                   |
| `raw/submissions/...csv`                              | GEM × Adaptyv entry form export. `author_email` column redacted before commit.                                      |
| `target/rbx1.fasta`                                   | UniProt P62877.                                                                                                     |

## Reported separately

- Zn²⁺-buffer rerun on the full 322 cohort.
- Ovalbumin specificity panel on the 9 confirmed binders.
- Developability assays (HIC, HPLC-SEC, nanoDSF).

Band-level expression labels from gel imaging aren't in the
ProteinBase release — only the binary `is_expressed` flag is.
