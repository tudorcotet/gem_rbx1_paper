# `leaderboard/` — the 9 binders, ordered

The KD table for the binder cohort plus a per-team and per-method
ranking. Reads `data/designs.csv` and writes:

- `top_binders.csv` — every confirmed binder, ordered by `kd_nM_mean`.
- `team_winners.csv` — best KD per team.
- `method_winners.csv` — best KD per method family.
- `report.md` — narrative with the leaderboard.

Run:

```bash
mise run analysis:leaderboard
```
