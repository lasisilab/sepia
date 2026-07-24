#!/bin/bash
# =====================================================================================
# SEPIA — archaic PROJECTION onto the modern PCA, for PIGMENTATION and WHOLE-GENOME.
# (This is Sanity check 2 for the WG side + the pigmentation projection to match Lily's method.)
# Builds a plink --export A dosage matrix (moderns + archaics, consistent allele coding), then
# code/project_archaics.R does prcomp(moderns) + predict(archaics). Outputs small coord TSVs.
#     cd /nfs/turbo/lsa-tlasisi1/sepia && git pull && sbatch code/project_archaics.slurm
# =====================================================================================
#SBATCH --job-name=sepia_project
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_project_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
SRC=/nfs/turbo/lsa-tlasisi1/lheald_thesis/aDNA_data
MM="$REPO/scratch/modern_merge"                              # sgdp166 bed + prune.in
WIDER="$REPO/scratch/wider_pca/wider.vcf.gz"                 # 159 SGDP + 5 archaics @ pig SNPs
AWG="$SRC/ancient_wholegenome/ancient_sgdp_wg.vcf.gz"        # archaic WG pseudo-haploid calls
SCR="$REPO/scratch/project"; OUT="$REPO/output/projections"
mkdir -p "$SCR" "$OUT"; cd "$SCR"

echo "===== SEPIA archaic projection (pig + WG) ====="; date
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21 >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
module load plink/2.0 >/dev/null 2>&1 || module load plink2 >/dev/null 2>&1 || true
module load R >/dev/null 2>&1 || true
PLINK2="$(command -v plink2 || command -v plink)"; RS="$(command -v Rscript)"
echo "plink2: $PLINK2 | Rscript: $RS | bcftools: $(command -v bcftools)"
[ -n "$PLINK2" ] && [ -n "$RS" ] || { echo "!! need plink2 + Rscript"; exit 2; }

# ---- 1. PIGMENTATION: export dosage from the wider merged VCF, project ---------------
echo "===== 1. pigmentation projection ====="
"$PLINK2" --vcf "$WIDER" --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
  --max-alleles 2 --snps-only --export A --out pig_dose
printf '%s\n' Altai Vindija33.19 Denisova3 Denisova25 Chagyrskaya8 Mez1 > pig_arc.txt
"$RS" "$REPO/code/project_archaics.R" pig_dose.raw pig_arc.txt "$OUT/pig_after"

# ---- 2. WHOLE GENOME: thin the pruned set, build moderns + archaics, merge, project --
echo "===== 2. whole-genome projection ====="
# 2a. thin ~483k pruned SNPs -> ~48k (ample for global structure); ids are chr:pos:ref:alt
awk 'NR%10==0' "$MM/sgdp166_prune.prune.in" > wg_thin.ids
awk -F: '{print $1"\t"$2}' wg_thin.ids | sort -k1,1 -k2,2n -u > wg_thin.pos
echo "  thinned WG SNPs: $(wc -l < wg_thin.ids)"
# 2b. moderns: subset sgdp166 bed to the thinned SNPs, export as VCF (consistent ref/alt)
"$PLINK2" --bfile "$MM/sgdp166" --extract wg_thin.ids --recode vcf bgz --out wg_mod
bcftools index -f -t wg_mod.vcf.gz
# 2c. archaics: clean the messy sample names to match Lily's, normalise contigs, extract
bcftools query -l "$AWG" | awk 'BEGIN{
    m["den11"]="Denisova11"; m["den5"]="Denisova5"; m["goyet"]="Goyet"; m["hst"]="Hohlenstein_Stadel";
    m["lesCottes"]="Les_Cottes"; m["mez2"]="Mezmaiskaya2"; m["Sid1253"]="El_Sidron"; m["spy"]="Spy";
    m["vindija87"]="Vindija87" }
  { old=$0; n=$0; sub(/.*\//,"",n); sub(/_sorted\.bam$/,"",n); if(n in m) n=m[n]; print old"\t"n }' > arc_rename.txt
CHRPREFIX=$(bcftools view -h "$AWG" | grep -qE '^##contig=<ID=chr' && echo yes || echo no)
if [ "$CHRPREFIX" = yes ]; then for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > chr2bare.txt; fi
( bcftools reheader -s arc_rename.txt "$AWG" \
    | bcftools view -T wg_thin.pos -Ou \
    | { [ "$CHRPREFIX" = yes ] && bcftools annotate --rename-chrs chr2bare.txt -Ou || cat; } \
    | bcftools norm -d none -Oz -o wg_arc.vcf.gz ) 2>wg_arc.err || { echo "!! archaic WG extract failed"; tail -3 wg_arc.err; }
bcftools index -f -t wg_arc.vcf.gz
echo "  archaic WG samples: $(bcftools query -l wg_arc.vcf.gz | tr '\n' ' ')"
# 2d. merge moderns + archaics (same hs37d5 ref -> matches by CHROM/POS/REF/ALT)
bcftools merge wg_mod.vcf.gz wg_arc.vcf.gz -Oz -o wg_merged.vcf.gz 2>wg_merge.err \
  || { echo "!! merge failed"; tail -3 wg_merge.err; exit 4; }
"$PLINK2" --vcf wg_merged.vcf.gz --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
  --max-alleles 2 --snps-only --export A --out wg_dose
# 2e. project (archaic list = the cleaned high-cov + all archaic names present)
bcftools query -l wg_arc.vcf.gz > wg_arc.txt
"$RS" "$REPO/code/project_archaics.R" wg_dose.raw wg_arc.txt "$OUT/wg_after"

# ---- 3. validation: high-cov archaics should sit off the modern cloud (non-African side) --
echo "===== 3. sanity: WG archaic PC1 vs modern range ====="
[ -f "$OUT/wg_after_archaic.tsv" ] && { echo "  archaic PC1 (should be extreme vs moderns):"; \
  awk -F'\t' 'NR>1{print "   "$1": PC1="$2}' "$OUT/wg_after_archaic.tsv" | head; }
{ echo "# archaic projection (pig + WG) onto the modern PCA"; date;
  echo "pig: $(wc -l < $OUT/pig_after_archaic.tsv) archaic rows; WG: $(wc -l < $OUT/wg_after_archaic.tsv) archaic rows";
  echo "WG thinned SNPs: $(wc -l < wg_thin.ids)"; } > "$OUT/projection_run.txt"
echo "===== projection done ====="; date
echo "NEXT (login node): git pull --rebase && git add output/projections && commit && push; figure laptop-side."
