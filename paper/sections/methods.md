### Target

Human RBX1 (UniProt P62877, 108 aa). The exact sequence provided to
competitors is in `data/target/rbx1.fasta`. The recombinant construct
used for expression at the Adaptyv Foundry is the full-length
P62877 with a C-terminal purification tag (vendor catalogue pending).

### Submission and selection

198 teams submitted 12,707 designs through the GEM × Adaptyv 2026
portal. Each team was capped at 100 ranked sequences (≤250 aa, ≤75%
identity to known proteins, standard amino acids only). An expert
panel from the GEM workshop selected 322 designs across 48 distinct
teams for wet-lab validation.

### Wet-lab pipeline

Designs were synthesised, expressed, and screened at the Adaptyv
Foundry. Expression was assessed on Coomassie-stained gels and BLI
loading checks. Binding was measured by surface plasmon resonance
(Carterra LSA) with 2–3 replicates per design (943 total SPR records
across 318 designs with replicates; three designs additionally have
BLI replicates). Per-replicate KD fits were averaged after curated-
replicate filtering at ProteinBase.

### Binding-strength bins

ProteinBase classifies each design into the standard hierarchy:
`Strong > Medium > Weak > Potential binder > Non-binder > No
expression > Unknown`. We treat `binding_strength ∈ {Strong, Medium,
Weak}` as a binder for headline counts (`is_binder` in
`data/designs.csv`). The *Strong* bin maps to KD < 100 nM.

### In silico

Every design was re-folded by ProteinBase using ESMFold; ProteinMPNN
log-likelihood scores and redesign recovery were computed against the
ESMFold backbone. Sequence and structural novelty were assessed via
mmseqs2 sequence identity and Foldseek TM-score against AFDB50.

No complex re-fold against the RBX1 target is included in the
ProteinBase release; Boltz-2, Chai-1, and Protenix-v2 complex
predictions are computed separately in this repository (figS4).

### Statistics

We report headline counts and rates with raw n/N. Where binary
comparisons were tested (e.g. method-family hit rates), we used
Fisher's exact (two-sided) and report the p-value alongside the odds
ratio. With 9 binders across 322 designs, almost every per-method or
per-modality split is underpowered.

### Code and data

All sequences, ESMFold predictions, and SPR sensorgrams are at
ProteinBase under ODC-ODbL. The canonical 322-row design table,
complex re-folds (Boltz-2, Chai-1, Protenix-v2), and analysis code
are released alongside this paper under MIT (code) and CC-BY-4.0
(data and figures).
