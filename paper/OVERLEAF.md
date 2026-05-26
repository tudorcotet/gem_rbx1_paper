# Overleaf sync

GitHub stays canonical; Overleaf is a thin mirror where collaborators leave
comments and track-changes edits. The `.github/workflows/overleaf-sync.yml`
workflow keeps the two in sync.

## One-time setup

1. **Overleaf Premium Standard** ($16.60/mo billed annually). The Git remote +
   track-changes / comments features are Premium-only. Without Premium, this
   workflow is a no-op (the secrets check exits cleanly).

2. **Create a blank Overleaf project** named `gem_rbx1_paper`. Open it once
   and note the project ID — it's the last path segment of the URL, e.g.
   `https://www.overleaf.com/project/68472c8f1a9e3b00abcdef01` → the ID is
   `68472c8f1a9e3b00abcdef01`.

3. **Generate a Git token**: Overleaf → Account Settings → *Git Integration*
   → *Create token*. The username is the literal string `git`; the password
   is the token. The same token works across all your projects.

4. **Add two GitHub repo secrets** (Settings → Secrets and variables →
   Actions → New repository secret):
   - `OVERLEAF_PROJECT_ID` — the project ID from step 2
   - `OVERLEAF_TOKEN` — the token from step 3

5. **First push.** Trigger the workflow once manually (Actions → *Overleaf
   sync* → *Run workflow*) so the initial paper-subtree push lands on
   Overleaf as the `master` branch.

## What happens after that

| trigger | action |
|---|---|
| push to `paper/**` or `figures/paper/**` on `main` | re-render pandoc → `paper/_build/*.tex`, commit if changed, push `paper/` subtree to Overleaf |
| hourly cron (xx:17) | fetch Overleaf into a `from-overleaf` branch; you merge manually via PR if collaborator edits happened |

## Local workflow

```bash
# Edit the markdown source
$EDITOR paper/sections/results.md

# Re-render to LaTeX locally (pandoc is pinned in mise.toml [tools])
mise run paper:render

# Push to GitHub — the workflow takes over and mirrors to Overleaf
git add paper/sections/results.md paper/_build/results.tex
git commit -m "..."
git push
```

## Why not Overleaf's built-in GitHub sync?

Two killers:

- It maps the *whole* GitHub repo to a flat Overleaf project (no
  subdirectory support). Our `paper/` lives inside a larger data repo.
- You can't link an existing Overleaf project to an existing GitHub repo —
  only "import GitHub → new Overleaf" or "Overleaf → new GitHub". Neither
  matches our setup.

The Git-integration path (used here) sidesteps both.

## Troubleshooting

- **Force-push collisions.** The workflow force-pushes `paper/` to Overleaf
  on every commit. The hourly pull job catches collaborator edits *before*
  the next push overwrites them. If a collaborator edits during the gap,
  their commit is preserved on `from-overleaf` — merge that branch into
  `main` before pushing new local changes.

- **`paper/_build/*.tex` out of sync.** The workflow re-renders pandoc and
  commits before pushing. If you edited the `_build/*.tex` directly,
  changes will be overwritten — edit the `paper/sections/*.md` source.

- **Workflow no-ops with no error.** The secrets check exits cleanly if
  `OVERLEAF_PROJECT_ID` / `OVERLEAF_TOKEN` aren't set, so disabling sync is
  as simple as removing the secrets.
