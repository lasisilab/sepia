#!/bin/bash
# =====================================================================================
# SEPIA — export the corrected-coordinate pigmentation-panel matrices for Lily's #1/#2/#3.
#
# Source: scratch/pig_panel/pig_panel.vcf.gz — the B2a extraction (159 SGDP + 5 high-cov
# archaics) at the CORRECT hg19 panel coordinates. Produces two small, redistributable
# per-sample matrices (all 159 SGDP samples are FullyPublic; the 5 archaics are published
# snpAD calls), which drive the depth + ATGC-variability heatmaps locally:
#   panel_hg19_genotypes.tsv  — base-level genotype (%TGT: e.g. G/A, G/G, ./.) per sample/locus
#   panel_hg19_depth.tsv      — read depth (%DP) per sample/locus
#   samples.txt               — sample order (159 SGDP then the 5 archaics)
#
# NB: snpAD VCFs put a non-integer '.' in the Integer GQ FORMAT field, which bcftools 1.21
# rejects on parse. We reheader GQ as String (header-only edit; body untouched) so the
# query works — the genotypes/depths themselves are unchanged. (Same GQ issue the archaic
# projection hit; there we side-stepped it with plink2.)
#
# Run from the cluster checkout:  bash code/export_panel_matrix.sh
# =====================================================================================
set -uo pipefail
REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
SRC="$REPO/scratch/pig_panel/pig_panel.vcf.gz"
OUT="$REPO/output/panel_matrix"; mkdir -p "$OUT"

module load Bioinformatics   >/dev/null 2>&1 || true
module load bcftools/1.21     >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
command -v bcftools >/dev/null || { echo "!! bcftools not found"; exit 2; }
[ -f "$SRC" ] || { echo "!! source VCF not found: $SRC (run code/pig_panel_extract.slurm first)"; exit 2; }

# --- GQ Integer -> String so bcftools tolerates the snpAD '.' values ---
bcftools view -h "$SRC" > "$OUT/hdr.txt"
sed -i 's/##FORMAT=<ID=GQ,Number=1,Type=Integer/##FORMAT=<ID=GQ,Number=1,Type=String/' "$OUT/hdr.txt"
bcftools reheader -h "$OUT/hdr.txt" "$SRC" -o "$OUT/pig_panel.fixed.vcf.gz"
bcftools index -f "$OUT/pig_panel.fixed.vcf.gz"

bcftools query -l "$OUT/pig_panel.fixed.vcf.gz" > "$OUT/samples.txt"
{ printf 'chrom\tpos\tref\talt'; while read s; do printf '\t%s' "$s"; done < "$OUT/samples.txt"; printf '\n'; \
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' "$OUT/pig_panel.fixed.vcf.gz"; } > "$OUT/panel_hg19_genotypes.tsv"
{ printf 'chrom\tpos'; while read s; do printf '\t%s' "$s"; done < "$OUT/samples.txt"; printf '\n'; \
  bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$OUT/pig_panel.fixed.vcf.gz"; } > "$OUT/panel_hg19_depth.tsv"

echo "samples=$(wc -l < "$OUT/samples.txt")  geno_rows=$(( $(wc -l < "$OUT/panel_hg19_genotypes.tsv") - 1 ))  depth_rows=$(( $(wc -l < "$OUT/panel_hg19_depth.tsv") - 1 ))"
rm -f "$OUT/hdr.txt" "$OUT/pig_panel.fixed.vcf.gz" "$OUT/pig_panel.fixed.vcf.gz.csi"
