# `overview/` — counts, rates, modality × method

The headline-number generator. Reads `data/designs.csv` and writes:

- `summary.json` — every number the abstract cites.
- `report.md` — a one-page narrative for the team.
- `modality_x_method.csv` — counts cross-tabulated by `pb_design_class`
  and `method_family`.
- `hit_rate_by_method.csv` — hit rate per method family with raw n/N.

Run:

```bash
mise run analysis:overview
```

If a number in the abstract doesn't trace back here, either change the
abstract or add the number to `summary.json`.
