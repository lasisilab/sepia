# SEPIA — collaboration worklist (internal)

**Internal working doc — not published to the site** (the publish workflow only copies
`plan/*.html` + figures, so `plan/*.md` stays off `lasisilab.github.io/sepia`). This tracks the
split of work between the cluster analyses (Tina/Claude) and Lily's pieces, with links.

Companion to [`plan.md`](plan.md) (roadmap), [`STATE.md`](STATE.md) (handoff), and
[`changelog.md`](changelog.md) (log). Last updated 2026-08-17.

---

## The split

The heavy extraction runs on Great Lakes (needs cluster access + big VCFs/BAMs Lily's laptop
can't hold), so **Tina/Claude run those**. The **modeling, interpretation, and control-set
assembly are Lily's** — they're where the science gets decided, and they run locally on the small
matrices below. Lily is *not blocked* on the cluster items; she has a full workstream from the
data already produced.

---

## Done — build on these, don't redo (links)

| Result | Where |
|---|---|
| Corrected pigmentation-panel PCA + the 3-way projection (whole-genome vs pigmentation-SNP vs pigmentation-gene) + the ascertainment finding | [`archaic_projection.qmd`](../archaic_projection.qmd) → live: https://lasisilab.github.io/sepia/archaic_projection.html |
| Extended 1,155-SNP map (confidence-tagged) + wider PCA | [`snp_panel.qmd`](../snp_panel.qmd), [`data/pigment_snp_map.tsv`](../data/pigment_snp_map.tsv) |
| **Corrected-coordinate panel matrices** (genotype + depth at right hg19 coords, 159 SGDP + 5 high-cov archaics) — *Lily's substrate* | [`output/panel_matrix/`](../output/panel_matrix/): `panel_hg19_genotypes.tsv`, `panel_hg19_depth.tsv`, `README.md` |
| Low-coverage archaic panel depth (Denisova 11, Goyet, Les Cottés) via remote ENA extraction | `output/panel_matrix/lowcov_panel_depth.tsv` (landing) |
| First-pass coverage figures (reference/validation only — Lily to own the finals) | [`plan/figures/`](figures/): `panel_atgc_heatmap.png`, `panel_depth_heatmap.png` |

---

## Lily's worklist

### 1. Genotype-confidence models on the corrected depth  *(her core piece — unstarted)*
Her binomial consensus-accuracy + Monte-Carlo damage models (currently in
[`depth.qmd`](../depth.qmd)) were built on the **wrong-coordinate** depth. Re-run them on the
corrected per-locus depth (`output/panel_matrix/panel_hg19_depth.tsv` for the high-cov archaics;
`lowcov_panel_depth.tsv` for the low-cov) to produce, **per archaic per panel SNP, the
confidence that the genotype call is correct given the observed depth**.
- *Done when:* a per-call confidence table + figure over the panel loci, so downstream we know
  which archaic calls to trust.
- *Why hers:* it's the reliability layer that makes the coverage figures actionable, and it's her
  model — nobody else should re-derive it.

### 2. Fix the missingness model  *(stats cleanup)*
In [`depth.qmd`](../depth.qmd) the `glmer(Missing ~ factor(Chromosome) + age + (1|Coverage) +
(1|Sample))` has a **degenerate `(1|Coverage)`** random effect: only 2 groups, so its variance
(≈13) is statistically unidentifiable and dwarfs the Sample variance (≈1). Make Coverage a
**fixed** effect (or drop it), and reconcile the reported coefficients (the prose numbers come
from a *numeric*-Chromosome model, not the `factor()` model actually fit). See the note in
[`pipeline_workbook.md`](pipeline_workbook.md).
- *Done when:* a non-degenerate model + corrected coefficients in the text.

### 3. Own the coverage figures + interpretation
Using `output/panel_matrix/` (genotype + depth matrices), produce the **final** variability and
depth figures your way, and — the real deliverable — **interpret** them: which panel loci / genes
are well vs poorly covered in which archaics, and what that implies for call confidence (ties to
#1). The `plan/figures/panel_*` PNGs are a first-pass reference to replace/refine.
- *Done when:* final figures + a short written interpretation (candidate text for the coverage
  section of the paper).

### 4. Assemble the immunity / neutral control gene set  *(your 3-PCA idea)*
For the "unrelated control" PCA, produce a **gene-symbol list** from a citable source
(ImmPort, InnateDB, or GO:0002376 "immune system process"), **plus** a random size-matched
control set. This feeds a control-region PCA (the cluster runs the PCA once the list exists),
built the same way as the pigmentation gene map ([`code/build_pigment_gene_map.py`](../code/build_pigment_gene_map.py)).
- *Note:* immune loci are themselves strongly selected (HLA, OAS, TLR) and introgression-rich, so
  frame immunity as a **second selected-function contrast**, not a neutral null — the random set
  is the cleaner neutral baseline. Both are worth running.
- *Done when:* `immune_genes.csv` (+ a random-control list) + a one-paragraph source note.

### 5. (Stretch) Methods note on coordinate handling
Short paragraph on reading an rsID panel at hg19 for physically-hg19 data (for the eventual paper's
methods) — documents the correction cleanly without re-litigating it.

---

## Cluster side (Tina/Claude — Lily not blocked on these)
- Low-coverage archaic panel depth (remote ENA extraction) — **running / landing.**
- Control-region PCA — runs once Lily's gene list (#4) exists.
- Chagyrskaya `.noRB` re-extraction (gets the 6th high-cov archaic into the panels).
- Full-BAM re-acquisition for the aDNA authenticity QC (ANGSD heterozygosity/contamination/sex +
  mapDamage, Q7/Q8) — the remote panel-depth trick doesn't cover these; they need whole BAMs.
