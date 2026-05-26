# `data/designs.csv` — the canonical table

One row per submitted design. **322 rows, ~47 columns.** Every analysis
in the repo reads from here. The parquet sibling has typed columns and
is the right thing to load in code; the CSV is the readable copy for
spreadsheet eyes.

`scripts/data/build_designs.py` builds this file from the raw
ProteinBase export under `data/raw/proteinbase/`. Re-run `mise run
build` after pulling a fresh export.

## Load it

```python
from scripts.utils import load_designs

df = load_designs()                       # 322 rows, every column
df = load_designs(only_expressed=True)    # 255 rows — designs that expressed
df = load_designs(only_binders=True)      #   9 rows — confirmed binders
df = load_designs(only_strong=True)       #   1 row — KD < 100 nM
```

`load_designs()` coerces booleans (`is_binder`, `is_expressed`, …) to
pandas nullable `boolean` so `df.is_binder & df.is_expressed` works
without object-dtype warnings.

## Keys

- **`design_id`** (int, 1..322) — universal join key. Stable across
  rebuilds. **Use this everywhere.**
- **`pb_id`** (string) — ProteinBase slug, e.g. `small-vole-maple`.
  Stable across releases. Use to fetch the structure / sensorgram from
  ProteinBase.

## Column groups

### 1. Identity (322/322)

| column            | type    | notes |
|---                |---      |---|
| `design_id`       | int     | 1..322. Universal join key. |
| `pb_id`           | string  | ProteinBase slug. |
| `name`            | string  | Submission name (team's own design tag). |
| `team`            | string  | Author handle. |
| `is_control`      | bool    | True for platform-authored spike-in controls. 0 in the public RBX1 release. |
| `design_method`   | string  | Self-reported method, verbatim. |
| `method_family`   | string  | Normalised bucket: `RFdiffusion`, `Mosaic`, `BoltzGen`, `BindCraft 2`, `MoPPIt`, `LFM2`, `AF2 hallucination + ADFlip`, … See [§9](#9-method-family-normalisation). |
| `sequence`        | string  | Amino-acid sequence, uppercase. |
| `sequence_length` | int     | `len(sequence)`. Range ~25..250. |

### 2. Upstream coverage flags (322/322)

Not every design in the public ProteinBase release received every
upstream evaluation. The metrics group into tiers; some are expected
to be sparse (e.g. KD only fits for the 9 binders), some flag genuine
pipeline gaps. These boolean columns let analyses filter without
guessing what's null and why.

| column                     | type | n=True | meaning |
|---                         |---   |---:    |---|
| `tested_in_wet_lab`        | bool | 321    | Any replicate ran. False for 1 design that entered ProteinBase but never went on a chip. |
| `has_wetlab`               | bool | 321    | Mirror of `tested_in_wet_lab`; True iff `expressed` had any record. |
| `has_predictions`          | bool | 322    | ESMFold pass ran. 318 came from the original ProteinBase release; the other 4 were backfilled by a local ProteinTyper rerun (`mise run rerun:typer`). The `pb_predictions_source` column distinguishes (`proteinbase_release` / `local_rerun`). |
| `pb_predictions_source`    | str  | —      | `proteinbase_release` (318) or `local_rerun` (4). |
| `has_homology`             | bool | 264    | AFDB50 search returned a hit. The 58 designs with `False` are novel enough that no AFDB50 neighbour was found — **expected**, not a bug. |
| `has_ted_classification`   | bool | 276    | TED could call a CATH class. The 46 designs with `False` are short / unusual folds TED can't characterize — expected. |
| `pb_n_metric_kinds`        | int  | —      | Raw count of distinct metric kinds for this design (4–22). Use as a one-shot coverage signal. |

The `pb_n_metric_kinds` histogram of upstream coverage (before any
local rerun is merged in):

| metrics shipped | designs | which tier |
|---:             |---:     |---|
| 4               | 4       | Only the wet-lab tier ran upstream — backfilled locally via ProteinTyper rerun. |
| 12              | 3       | Wet-lab + partial computational tier. |
| 15              | 40      | Wet-lab + ESMFold + ProteinMPNN + novelty, no TED, no AFDB50 hit. |
| 18              | 12      | One of TED or AFDB50 missing. |
| 19              | 255     | **Modal case** — every standard metric except `kd/koff/kon` (non-binders). |
| 22              | 8       | Full bundle — includes `kd/koff/kon` (binders with a fittable curve) + `spr_kinetic_curves`. |

Backfilling those 4 designs is done by `mise run rerun:typer`, which
calls the Modal `proteintyper-submit` endpoint with the **default
recipe** — the same shape ProteinBase's `workers/proteintyper-receiver`
uses (no `target`, no `webhook_*`, `recipe.template = "full_monomer"`).
Outputs land at `data/metrics/proteintyper/<pb_id>.json` with the CIF
and stylised PNG downloaded into `data/structures/esmfold/<pb_id>.cif`
and `data/images/<pb_id>.png`. `build_designs.py` merges them into the
canonical row at build time.

### 3. Wet-lab outcome (322/322; 9 with binding)

The wet-lab pipeline ran each design through expression + SPR in default
buffer. Three to five replicates per design (`n_replicates`).
ProteinBase aggregates the per-replicate hits.

| column                   | type   | notes |
|---                       |---     |---|
| `binding_strength`       | string | Best across replicates: `Strong` (1), `Medium` (7), `Weak` (1), or null (313 non-binders). |
| `is_binder`              | bool   | `binding_strength ∈ {Strong, Medium, Weak}`. **9 binders.** |
| `is_strong`              | bool   | `binding_strength == "Strong"`. **1 design.** |
| `any_expressed`          | bool   | True if any replicate expressed. |
| `is_expressed`           | bool   | Boolean coercion of `any_expressed`. **255 / 322 ≈ 79%.** |
| `any_binding`            | bool   | True if any replicate showed binding (`binding=True` from SPR). |
| `n_replicates`           | int    | Total replicate rows from ProteinBase. |
| `n_replicates_expressed` | int    | Replicates with `expressed=True`. |
| `n_replicates_binding`   | int    | Replicates with `binding=True`. |

### 4. Affinity (9/322 — only for binders with a fit)

| column            | type  | notes |
|---                |---    |---|
| `kd_M_mean`       | float | Mean of per-replicate KD fits (in **molar**). |
| `kd_M_min`        | float | Tightest replicate. |
| `kd_M_max`        | float | Loosest replicate. |
| `kd_nM_mean`      | float | `kd_M_mean * 1e9`. **Primary KD column.** |
| `kd_nM_min`       | float | Same in nM. |
| `kd_nM_max`       | float | Same in nM. |
| `pkd_arith_mean`  | float | `-log10(kd_M_mean)`. **Use for stats.** |
| `koff_mean`       | float | Mean replicate off-rate (1/s). |
| `kon_mean`        | float | Mean replicate on-rate (1/Ms). |
| `n_kd_records`    | int   | Replicates that produced a fittable curve. |

### 5. Modality / fold class (322/322 for `pb_design_class`; 283/322 for CATH `pb_classification`)

| column            | type   | notes |
|---                |---     |---|
| `pb_design_class` | string | ProteinBase modality call: `Miniprotein` (79), `Nanobody` (33), `Peptide` (32), `scFv` (10), `Other` (164). The "Other" bucket is dominant — don't claim modality trends without controlling for it. |
| `pb_classification`| string | CATH-style class: `Mainly Alpha`, `Mainly Beta`, `Alpha Beta`. |
| `pb_foldstring`   | string | Per-residue secondary-structure string (`HHHH…`). |
| `pb_ted_confidence`| float | TED domain-call confidence. |

### 6. Sequence-derived features (322/322)

| column                  | type   | notes |
|---                      |---     |---|
| `pb_molecular_weight`   | float  | Daltons (ProteinBase recomputation). |
| `pb_isoelectric_point`  | float  | pI. |

### 7. In-silico folding metrics (322/322)

| column                              | type  | notes |
|---                                  |---    |---|
| `pb_esmfold_plddt`                  | float | Mean pLDDT from ESMFold on the binder alone (%). |
| `pb_proteinmpnn_score`              | float | ProteinMPNN log-likelihood of the submitted sequence given the fold. |
| `pb_proteinmpnn_seq_recovery`       | float | Fraction of positions where MPNN redesign agrees with the submitted sequence. |
| `pb_redesigned_proteinmpnn_score`   | float | Score of the MPNN-redesigned sequence on the same backbone. |

**No Boltz / AF3 metrics in the ProteinBase release.** Teams ran
their own in silico — the upstream release did not re-fold against
the target. Boltz-2, Chai-1, and Protenix-v2 complex predictions for
all 322 designs are computed in this repository and pooled into
`data/grand_metrics.csv` (`mise run build:grand`).

### 8. Novelty / homology (varies — filter on `has_homology`)

| column                      | type   | notes |
|---                          |---     |---|
| `pb_novelty`                | float  | ProteinBase composite novelty score. |
| `pb_seqidentity`            | float  | Sequence identity to the closest known protein (%). Method-internal. |
| `pb_seqidentity_afdb50`     | float  | Identity vs AFDB50 closest hit. From `domainmatch`. |
| `pb_evalue_afdb50`          | float  | E-value vs AFDB50 closest hit. |
| `pb_tm_score_afdb50`        | float  | TM-score vs AFDB50 closest hit. |
| `pb_afdb50_top_id`          | string | AFDB50 hit identifier. |

### 9. Artifact URLs (322/322 ESMFold CIF; 319/322 stylised PNG; subset for SPR)

| column                  | type   | notes |
|---                      |---     |---|
| `pb_esmfold_cif_url`    | string | `https://proteinbase-pub.t3.storage.dev/<ulid>.cif`. ESMFold prediction. |
| `pb_stylized_png_url`   | string | Stylised cartoon render PNG. |
| `pb_spr_curves_url`     | string | SPR sensorgram JSON URL. Only for designs with replicates. |

All artifact URLs are dereferenced into Git LFS via
`mise run mirror:structures`.

### 10. Submission metadata (293/322)

The competition entry form stored a per-submission method writeup
separate from ProteinBase. `scripts/data/build_designs.py` joins it
onto each design row by lower-cased team handle, with a hand-curated
override map in `_HANDLE_OVERRIDES` for the cases where the
ProteinBase handle doesn't normalise the same way as the form
`author_name` (e.g. `hz3519` ↔ `Haowen Zhao`, verified via the
shared method token `giraf`).

**293 of 322 designs are matched.** The remaining 29 break down as:

| reason                                                                              | designs |
|---                                                                                  |---:|
| Team filed no entry form, designMethod blank in ProteinBase too                     | 21 |
| Team filed no entry form, designMethod is `boltzgen` (too ambiguous to rescue)      |  7 |
| Team filed no entry form, custom workflow (`drtheone`, 1-off)                       |  1 |

Those 21 truly unknown designs come from three teams (`nanogenomic`,
`professionalmouthpipettor`, `zhangpeioo`); their method is captured
as `method_family = "Not mentioned"` and `method_family_source = "missing"`.

| column                          | type   | notes |
|---                              |---     |---|
| `submission_author_name`        | string | Display name from the submission form. |
| `submission_modality`           | string | Declared modality: `minibinder`, `nanobody`, `peptide`, combos. |
| `submission_target_region`      | string | Declared target face: `RING domain`, `IDR/N-terminus`, `E2-face`, or combos. |
| `submission_targets_idr`        | string | yes / no / both / unknown. **27 submissions explicitly target the disordered tail.** |
| `submission_binder_length`      | string | Declared length range. |
| `submission_team_type`          | string | `independent` (140), `academic` (43), `industry` (17). |
| `submission_new_method`         | string | yes (56) / no (144). |
| `submission_design_type`        | string | `de novo`, `motif scaffolding`, etc. |
| `submission_core_models`        | string | Comma-separated models named by the author. |
| `submission_method_summary`     | string | One-line method summary. |
| `submission_link`               | string | ProteinBase per-submission collection URL. |
| `submission_uniref_check`       | string | UniRef50 ≥25% edit-distance check pass/fail. |
| `submission_sabdab_check`       | string | SAbDab CDR ≥25% edit-distance check (antibody-only). |
| `submission_overall_homology`   | string | Overall pass / partial / fail. |
| `submission_n_proteins`         | int    | Number of sequences in the submission (1–100). |
| `submission_total_submissions`  | int    | How many separate submissions the author filed. |

### 11. Method-family normalisation

The `designMethod` field is free text. `scripts/data/build_designs.py`
runs a regex map to canonicalise it into `method_family`. Current
buckets (with first-pass counts):

| bucket                       | n   | comment |
|---                           |---  |---|
| `Not mentioned`              | 93  | Empty `designMethod`. Many top teams left it blank. |
| `Other`                      | 58  | Pattern didn't match — add to `_METHOD_FAMILY_PATTERNS`. |
| `RFdiffusion`                | 27  | Largest single named bucket. |
| `Mosaic`                     | 25  | Escalante's composite-objective wrapper. |
| `BoltzGen`                   | 22  | Boltz-2 inverse-design extension. |
| `MoPPIt`                     | 20  | PPI optimiser. |
| `LFM2`                       | 14  | Liquid AI's LFM2 with method-specific customisation. |
| `AF2 hallucination + ADFlip` | 7   | AF2 backprop + sequence editor. |
| `BindCraft 2`                | 7   | Pacesa's unpublished branch. Cite as personal comm. |
| `Bagel + SoluMPNN`           | 7   | … |
| `FoldCraft`                  | 7   | … |
| `PepMind + AF3`              | 7   | … |
| `ORBIT`                      | 1   | Mandrake's stack — the winner. |

`method_family_source` traces where each label came from:

| source                       | n   | meaning |
|---                           |---:|---|
| `proteinbase_designMethod`   | 229 | ProteinBase carried a non-empty `designMethod`. |
| `submission_core_models`     |  58 | ProteinBase was blank; the submission form's `core_models` tool list matched a pattern. |
| `submission_method_summary`  |  14 | Both above were blank; the prose summary mentioned a known tool. |
| `missing`                    |  21 | Nothing populated anywhere — three teams left both surfaces blank. |

When `Not mentioned` shows up where you expected something named, add
a pattern to `_METHOD_FAMILY_PATTERNS` in
`scripts/data/build_designs.py` and re-run `mise run build`.

## Joining with the raw export

| from               | to                                                  | key |
|---                 |---                                                  |---|
| `data/designs.csv` | `data/raw/proteinbase/...csv`                       | `pb_id` |
| `data/designs.csv` | `data/proteinbase/esmfold/<pb_id>.cif`              | `pb_id` (after LFS mirror) |
| `data/designs.csv` | `data/proteinbase/sensorgrams/<pb_id>_rep*.json`    | `pb_id` (after LFS mirror) |

## Re-scoring (rbx1 paper extras)

The public ProteinBase release only ran a monomer typer (ESMFold +
ProteinMPNN + AFDB50). The paper-extra rerun re-folds every design
against the RBX1 target with three predictors, mirroring the
complex re-scoring panel. Each model gets the
**same target MSA** (computed once for RBX1) and a single-sequence
binder — matching ProteinTyper's `msa_mode: target_only`. Outputs are
read by `scripts/data/build_grand_metrics.py` into `data/grand_metrics.csv`.

### Models

| model      | weight                                  | bundle source                                  | notes |
|---         |---                                      |---                                             |---|
| Boltz-2    | `boltz==2.2.1`                          | `scripts/modal/modal_boltz2_rbx1.py`            | Uses ColabFold's MSA server for the target. Binder pinned to `msa: empty`. |
| Protenix   | `protenix-v2` (default, +9–13 pp DockQ over v1) | `scripts/modal/modal_protenix_rbx1.py` | The v2 checkpoint is gated on ByteDance's volces.com host ([issue #295](https://github.com/bytedance/Protenix/issues/295)); the Modal app pre-fetches it from `TMF001/pxdesign-weights` on HuggingFace and pre-places it at `/protenix-models/checkpoint/protenix-v2.pt` so Protenix's loader skips the volces.com path. Override with `PROTENIX_MODEL_NAME` to use `protenix_{base,mini}_*_v0.5.0` or `pxdesign_v0.1.0` (all on the same HF mirror). |
| Chai-1     | `chai_lab>=0.6.0`                       | `scripts/modal/modal_chai1_rbx1.py`             | ESM-embedding mode (`use_esm_embeddings=True`); no external MSA. |
| ProteinTyper | Modal endpoint (default recipe)       | `scripts/data/run_proteintyper.py`              | Monomer-only — same as ProteinBase. Re-exported columns appear as `tp_*` in `grand_metrics.csv`. |

### Disk layout

```
data/
├── structures/
│   ├── esmfold/<pb_id>.cif   ← from typer (monomer)
│   ├── boltz2/<pb_id>.cif    ← Boltz-2 complex
│   ├── protenix/<pb_id>.cif  ← Protenix complex
│   └── chai/<pb_id>.cif      ← Chai-1 complex
├── metrics/
│   ├── proteintyper/<pb_id>.json  ← full TyperJobOutput
│   ├── boltz2/<pb_id>.json        ← {pb_id, predictor, status, iptm, ptm, mean_plddt, ipsae_*, pdockq, pdockq2, lis, n_interface}
│   ├── protenix/<pb_id>.json      ← same shape + ranking_score, model_name
│   └── chai/<pb_id>.json          ← same shape + aggregate_score
└── grand_metrics.csv               ← one row per pb_id, all model prefixes
```

### `data/grand_metrics.csv` columns

Prefixes match the model bucket. Generated by
`scripts/data/build_grand_metrics.py`.

| prefix | source                                  | example columns |
|---     |---                                      |---|
| `tp_`  | ProteinTyper monomer panel              | `tp_esmfold_plddt`, `tp_proteinmpnn_score`, `tp_design_class`, `tp_novelty`, `tp_tm_score_afdb50` |
| `b2_`  | Boltz-2 complex                         | `b2_iptm`, `b2_ptm`, `b2_mean_plddt`, `b2_ipsae_d0chn_max`, `b2_pdockq2`, `b2_lis`, `b2_n_interface` |
| `px_`  | Protenix complex                        | `px_iptm`, `px_ptm`, `px_ipsae_d0chn_max`, `px_ranking_score`, `px_model_name` |
| `chai_`| Chai-1 complex                          | `chai_iptm`, `chai_ptm`, `chai_ipsae_d0chn_max`, `chai_aggregate_score` |

Plus derived consensus columns at the end:

| column                  | meaning |
|---                      |---|
| `ipsae_pass_3folders`   | Count of {b2, px, chai} with `ipsae_d0chn_max >= 0.4` (soft binder threshold). |
| `iptm_pass_3folders`    | Count of {b2, px, chai} with `iptm >= 0.7` (BindCraft default). |

Missing models leave nulls; the per-model `<prefix>_status` column
distinguishes `ok` / `failed_*` / `error: …` / null (not run yet).

### Rerun commands

```bash
# All three complex predictors, all 322 designs (Modal-detached)
mise run rerun:complex            # or:  mise run rerun:rescore

# Just Protenix (defaults to protenix-v2, pre-fetched from HuggingFace)
mise run rerun:protenix
PROTENIX_MODEL_NAME=protenix_base_default_v0.5.0 mise run rerun:protenix

# Pool the per-model JSONs into one wide CSV
mise run build:grand
```

The Modal apps are idempotent — they cache per-(pb_id, model) on a Modal
Volume named `rbx1-rerun-results` and skip work that's already done.

## What's NOT in this CSV

- **Boltz-2 / AF3 confidence on the target complex.** Public release is
  ESMFold-only on the binder alone. The re-scoring rerun above produces
  `data/grand_metrics.csv` separately — `designs.csv` itself stays
  monomer-only to match ProteinBase.
- **Per-replicate SPR traces.** Aggregated in the CSV. SPR sensorgram
  URLs are in `pb_spr_curves_url`.
- **Zn²⁺-buffer rerun.** Reported separately.
- **Ovalbumin specificity panel.** Reported separately.

## Wet-lab batches

| condition                  | status |
|---                         |---|
| Default buffer             | released via ProteinBase (322 designs) |
| Zn²⁺-buffer rerun (322)    | reported separately |
| Ovalbumin specificity (9)  | reported separately |

## Reserved columns

Column slots reserved for follow-up batches; null in this release.

| column                         | source                              |
|---                             |---                                  |
| `kd_nM_mean_zinc`              | Zn²⁺-buffer rerun                   |
| `is_binder_zinc`               | Zn²⁺-buffer rerun                   |
| `ovalbumin_response`           | specificity panel                   |
| `hic_retention_time`           | developability (HIC)                |
| `sec_aggregation_pct`          | developability (HPLC-SEC)           |
| `nanodsf_tm_C`                 | developability (nanoDSF)            |

Boltz-2, Chai-1, and Protenix-v2 complex predictions land in
`data/grand_metrics.csv`, not `designs.csv` — keep `designs.csv`
monomer-only to match the ProteinBase schema.

## Regeneration

`data/designs.csv` is rebuilt by `scripts/data/build_designs.py`. The
build is deterministic — running `mise run build` twice produces
byte-identical output. If you spot a wrong row, **don't edit the CSV in
place.** Open an issue with the offending `pb_id` and the corrected
value; the fix goes into the build script (or upstream into
ProteinBase).
