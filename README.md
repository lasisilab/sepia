# SEPIA

**Systematic Evaluation of Pigmentation In Archaic hominins**

SEPIA asks what skin/hair/eye **pigmentation genetics** can — and can't — tell us about
Neanderthals and Denisovans, and how they compare to modern humans. It is a careful,
QC-first re-analysis: the archaic DNA is patchy and shallow and pigmentation genes don't
behave like ancestry markers, so the contribution is an **honest map of what the
pigmentation genetics does and doesn't support**, with the uncertainty made explicit.

This repository is the working anchor for the project. It continues and re-evaluates the
prior honors thesis **PAINT** (Lily Heald), rebuilding the analysis from clearly-traced
inputs; Lily's original thesis repository is separate and unchanged.

## Where things are

- **[`plan/plan.md`](plan/plan.md)** — the single canonical roadmap (phased; bug fixes,
  analysis redesign, QC checks `Q1`–`Q10`, infrastructure, writing). Shareable rendering:
  **[`plan/sepia-plan.html`](plan/sepia-plan.html)**.
- **[`plan/changelog.md`](plan/changelog.md)** — dated log of what was done / found / decided.
- **[`plan/docs.html`](plan/docs.html)** — hub linking every project document (data &
  pipeline record, rebuild plan, per-genome provenance & QC, bug-evidence log).
- **[`data/README.md`](data/README.md)** — per-file data provenance, and where the large /
  cluster-only data live.

## Layout

The website and the compute pipeline live in one repo; on the cluster the **repo is the
working directory**, so results flow into git while big data stays out of it.

```
*.qmd  _quarto.yml   Quarto website (source + config) → GitHub Pages
analysis/            site helper scripts + the introgression Shiny app source
code/                cluster pipeline (SLURM/bash): merge, QC, PCA
plan/                roadmap, changelog, provenance & QC docs (+ rendered HTML)
data/                small tracked inputs (see data/README.md)
output/              small tracked RESULTS that sync (PCA coords, QC tables, figures)
resources/  scratch/ gitignored, cluster-only: genomes + big intermediates (never synced)
```

On Great Lakes the checkout lives at `/nfs/turbo/lsa-tlasisi1/sepia`; reusable reference
genomes live in a shared `…/genomes` dir (referenced via `$SEPIA_GENOMES` or the gitignored
`resources/genomes` symlink) so other lab studies can call the same data without copies.

## Build

- **Website (Quarto):** `quarto render`. R-heavy analysis pages use `freeze: auto` — they
  render from committed `_freeze/` in CI (which does not run R), so they must be rendered
  once locally/on the cluster (with R + the data) to populate `_freeze/`. Prose pages render
  anywhere. `.github/workflows/publish-site.yml` deploys the rendered `_site/` to GitHub
  Pages (currently set to manual `workflow_dispatch` until the first full render).
- **Project docs under `plan/`:** regenerate the HTML with `bash plan/build_docs.sh`.
- **Pipeline:** run `code/*.slurm` from the cluster checkout (`sbatch code/…`).

## Team

Tina Lasisi (PI / advisor) · Lily Heald · Yemko Pryor. Code under the MIT License; written
content and figures under CC BY 4.0.
