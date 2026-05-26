# Overleaf sync

Git is the source of truth. Overleaf is a mirror for comments and
track-changes. The `.github/workflows/overleaf-sync.yml` workflow runs
the two sides asymmetrically: GitHub → Overleaf is fully automatic,
Overleaf → GitHub lands on a branch you merge by hand.

## Direction matrix

| direction | when it fires | what happens | merge step |
|---|---|---|---|
| GitHub → Overleaf | on push to `paper/**` or `figures/paper/**` (plus manual dispatch) | workflow clones Overleaf, mirrors `paper/` in, copies `figures/paper/*.pdf` into `figures/`, commits, pushes Overleaf `master` | none — lands on Overleaf within ~1 min |
| Overleaf → GitHub | hourly cron at xx:17 (plus manual dispatch) | workflow fetches Overleaf `master`, snapshots it onto branch `from-overleaf` | open PR `from-overleaf` → `main` and merge |

The asymmetry is intentional. If both sides auto-merged, a GitHub push
during in-flight Overleaf edits would overwrite the collaborator's
work on the next mirror push. The `from-overleaf` branch is a
checkpoint: anything Overleaf-side lands there first, you decide when
to fold it back.

## Setup (once)

1. Buy Overleaf Premium Standard ($16.60/mo). The Git remote, comments,
   and track-changes are all paid.
2. Create a blank Overleaf project. Copy the URL tail (the project
   ID, e.g. `6a159c0a442f9f862f3bf20b`).
3. Overleaf → Account Settings → Git Integration → *Create token*.
   Username is `git`, password is the token.
4. GitHub repo → Settings → Secrets and variables → Actions. Add:
   - `OVERLEAF_PROJECT_ID` (the URL tail)
   - `OVERLEAF_TOKEN`
5. Actions → *Overleaf sync* → *Run workflow*. The first push seeds
   Overleaf with the current `paper/` tree plus figure PDFs.

## Edit on GitHub

```bash
$EDITOR paper/sections/results.tex
git add paper/ && git commit -m "..." && git push
```

The push triggers the workflow; Overleaf updates within a minute.

## Edit on Overleaf

Collaborators open the Overleaf project and comment / track-change.

Within an hour the cron snapshots their changes onto `from-overleaf`.
To fold them into `main`:

```bash
gh pr list --head from-overleaf      # check what's pending
gh pr create --base main --head from-overleaf
# review, merge
```

If you push to `main` while `from-overleaf` has unmerged edits, the
next mirror push overwrites Overleaf and the collaborator's work is
lost. Always check `from-overleaf` first, or run:

```bash
gh pr list --head from-overleaf --json url -q '.[0].url'
```

before every push to `main`.

## Manual triggers

```bash
gh workflow run "Overleaf sync"   # fires both jobs
gh run watch                      # follow the run
gh run view --log-failed          # debug a failed run
```

## What gets mirrored

Overleaf sees:

- `main.tex`
- `sections/*.tex`
- `references.bib`
- `figures/*.pdf` (copied from `figures/paper/` by the workflow before
  push; PNG and SVG variants stay local)

`main.tex` `\graphicspath` lists `./figures/` first then
`../figures/paper/`. The same `\includegraphics{fig1_workflow}` works
on Overleaf (finds `./figures/fig1_workflow.pdf`) and locally (finds
`../figures/paper/fig1_workflow.pdf`).

## Disabling sync

Delete `OVERLEAF_PROJECT_ID` and `OVERLEAF_TOKEN` from GitHub secrets.
The workflow exits cleanly.
