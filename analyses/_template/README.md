# `analyses/_template/` — copy me to start a new analysis

## 1. Copy the template into your own subdir

```bash
cp -r analyses/_template analyses/<your_handle>
$EDITOR analyses/<your_handle>/main.py
```

`<your_handle>` is whatever's easy to find: your name, your team, your
method. The subdir owns *all* your outputs. Don't write anywhere else.

## 2. Edit `main.py`

The shipped `main.py` is a working skeleton. The contract is:

1. **Load via `load_designs()`** — never read a path directly.
2. **Write to your own subdir** — `report.md`, `summary.json`, derived CSVs.
3. **Be idempotent** — running `main.py` twice gives the same output.

```python
from scripts.utils import load_designs

df = load_designs()                       # 322 rows, every column
df = load_designs(only_expressed=True)    # 255 rows that expressed
df = load_designs(only_binders=True)      #   9 confirmed binders
df = load_designs(only_strong=True)       #   1 Strong binder (KD < 100 nM)
```

Column reference: [`docs/DATA.md`](../../docs/DATA.md).

## 3. Register a `mise` task

Open the repo-root `mise.toml`, scroll to the "ADD YOURS BELOW" comment,
and add:

```toml
[tasks."analysis:<your_handle>"]
description = "One-sentence answer to: what does this analysis tell us?"
run = "uv run python analyses/<your_handle>/main.py"
```

If your analysis is load-bearing for the paper, also add
`"analysis:<your_handle>"` to the `depends = [...]` list of
`[tasks."analysis:all"]` so it runs on every `mise run analysis:all`.

## 4. Run it

```bash
mise run analysis:<your_handle>     # your analysis only
mise run analysis:all               # everyone's, end-to-end
```

## What this analysis answers

(Replace with your one-sentence framing.)

## Inputs

- `data/designs.csv` (via `scripts.utils.load_designs()`).

## Outputs

- `analyses/<your_handle>/report.md` — human-readable summary (≤1 page).
- `analyses/<your_handle>/summary.json` — machine-readable headline numbers.
- Any derived CSVs you write — keep them inside this folder.

## Don't

- Don't read `data/raw/proteinbase/*.csv` directly. The `evaluations`
  column is a nested JSON list; `scripts/data/build_designs.py` already
  unpacks it into flat columns.
- Don't modify `data/designs.csv`. Fix the build script instead.
- Don't write outside `analyses/<your_handle>/`.
