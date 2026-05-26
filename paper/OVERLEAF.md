# Overleaf sync

Git is the source of truth. Overleaf is a mirror for comments and
track-changes. The `.github/workflows/overleaf-sync.yml` workflow keeps both
sides current; you should rarely touch git or Overleaf directly to keep
them aligned.

## Setup (once)

1. Buy Overleaf Premium Standard ($16.60/mo). The Git remote, comments, and
   track-changes are all paid.
2. Create a blank Overleaf project. Copy the URL tail — that's the project
   ID, e.g. `6a159c0a442f9f862f3bf20b`.
3. Overleaf → Account Settings → Git Integration → *Create token*. Username
   is `git`, password is the token. One token works across every project.
4. GitHub repo → Settings → Secrets and variables → Actions. Add two
   secrets:
   - `OVERLEAF_PROJECT_ID` — the URL tail
   - `OVERLEAF_TOKEN` — the token
5. Actions → *Overleaf sync* → *Run workflow*. The first push seeds
   Overleaf with the current `paper/` tree.

## Edit on GitHub

```bash
$EDITOR paper/sections/results.md
mise run paper:render           # regenerates paper/_build/*.tex
git add paper/ && git commit -m "..." && git push
```

The workflow fires on every push that touches `paper/**` or
`figures/paper/**`. It re-renders pandoc, commits any `_build/*.tex`
delta back to `main`, then mirrors `paper/` into Overleaf's `master`.

## Edit on Overleaf

Collaborators open the Overleaf project, comment, and track-change.

The hourly cron (xx:17) pulls Overleaf's `master` into a `from-overleaf`
branch on GitHub. Review and merge:

```bash
gh pr create --base main --head from-overleaf
# review, merge
```

Don't push to `main` while there are pending edits on `from-overleaf` —
the next mirror push will overwrite Overleaf and erase them.

## Manual triggers

```bash
gh workflow run "Overleaf sync"                # both jobs
gh run watch                                   # follow the run
gh run view --log-failed                       # debug a failed run
```

## What gets mirrored

Overleaf sees exactly the contents of `paper/`:

- `main.tex` — entry point, references the rendered `_build/` files
- `sections/*.md` — markdown source (ignored by Overleaf's compiler;
  there for the human)
- `_build/*.tex` — rendered LaTeX, what Overleaf compiles
- `references.bib`

Figures are referenced via `\graphicspath{{../figures/paper/}}`. Overleaf
sees `figures/` only if it's copied into `paper/`; right now it isn't.
Either copy figure PNGs into `paper/figures/` for the Overleaf side, or
edit only the prose on Overleaf and let GitHub re-render the figures.

## Disabling sync

Delete `OVERLEAF_PROJECT_ID` and `OVERLEAF_TOKEN` from GitHub secrets.
The workflow becomes a no-op (the secrets check exits cleanly).
