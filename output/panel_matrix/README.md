# Pigmentation-panel per-sample matrices (corrected hg19 coordinates)

The substrate for Lily's coverage/variability analyses (#1–#3), **repaired**: the thesis figures measured the 222 panel SNPs at their **hg38** positions in **hg19** data (~292 kb off-target — bug A2). These matrices read the same panel at the **correct hg19 coordinates** and span **moderns + archaics**, not archaics alone.

## Files
- `panel_hg19_genotypes.tsv` — base-level genotype (`%TGT`, e.g. `G/A`, `G/G`, `./.` for missing) per sample × locus. 161 variable panel loci × 164 samples.
- `panel_hg19_depth.tsv` — read depth (`%DP`) per sample × locus (missing call = `.`).
- `samples.txt` — sample order: 159 SGDP (all `FullyPublic`) then the 5 high-coverage archaics (`Vindija33.19`, `Mez1`, `Denisova3`, `Altai`, `Denisova25`).
- `lowcov_panel_depth.tsv` — panel depth for the 3 low-coverage archaics (`Denisova11`, `Goyet`, `LesCottes`), summed over each individual's run libraries, streamed from ENA's complete files ([`code/lowcov_panel_depth.slurm`](../../code/lowcov_panel_depth.slurm)) since the local BAMs are truncated. Mean panel depth ≈ 9.2 / low / 2.3 respectively.

## Provenance
Exported by [`code/export_panel_matrix.sh`](../../code/export_panel_matrix.sh) from `scratch/pig_panel/pig_panel.vcf.gz` (the B2a extraction, [`code/pig_panel_extract.slurm`](../../code/pig_panel_extract.slurm)). Plotted by [`plan/figures/make_panel_heatmaps.R`](../../plan/figures/make_panel_heatmaps.R) → `plan/figures/panel_{atgc,depth}_heatmap.png`.

## Known gaps
- **Chagyrskaya 8** is absent (5 high-cov archaics, not 6) — the `.noRB` VCF extracts 0 records (known bug, to fix).
- The **3 low-coverage archaics** (Denisova 11, Goyet, Les Cottés) are BAM-only, so their depth is in the separate `lowcov_panel_depth.tsv` (remote ENA extraction); a per-sample *genotype* call for them still needs genotype likelihoods (ANGSD), not covered here.
- Loci that are monomorphic-reference across all samples are not emitted as variant records, so the matrices carry the **161 variable** panel loci (of 222).
