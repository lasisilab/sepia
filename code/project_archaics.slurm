#!/bin/bash
# =====================================================================================
# SEPIA — archaic PROJECTION onto the modern PCA (pigmentation + whole genome), the proper
# plink2 way (supersedes the earlier hand-rolled prcomp, which mis-sided the WG archaics).
#
# Recipe (plink2 PCA projection):
#   modern:  --freq counts  +  --pca allele-wts       -> .acount + .eigenvec + .eigenvec.allele
#   archaic: --read-freq modern.acount  --score modern.eigenvec.allele 2 5 header-read
#            no-mean-imputation variance-standardize --score-col-nums 6-15
# variance-standardize uses the MODERN allele frequencies, no-mean-imputation drops missing
# archaic genotypes, and the ID+A1 columns make allele coding explicit -> correct side + scale.
# VALIDATION: we also --score the MODERNS with their own weights and confirm it reproduces
# their --pca eigenvec (r>0.99); only then are the archaic projections trustworthy.
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
MM="$REPO/scratch/modern_merge"
WIDER="$REPO/scratch/wider_pca/wider.vcf.gz"
AWG="$SRC/ancient_wholegenome/ancient_sgdp_wg.vcf.gz"
P2="$REPO/tools/plink2"                                      # newer plink2 (allele-wts + score)
SCR="$REPO/scratch/project"; OUT="$REPO/output/projections"
mkdir -p "$SCR" "$OUT"; cd "$SCR"
echo "===== SEPIA archaic projection (plink2 recipe) ====="; date
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21 >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
"$P2" --version || { echo "!! newer plink2 not found at $P2"; exit 2; }
IDARGS=(--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing --max-alleles 2 --snps-only)

# helper: pull PC1/PC2 out of an .eigenvec (modern) or .sscore (archaic) by header name -> TSV
coords_eigenvec(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1")A=i;if($i=="PC2")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }
coords_sscore(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1_AVG")A=i;if($i=="PC2_AVG")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }

# ---- 1. PIGMENTATION -----------------------------------------------------------------
echo "===== 1. pigmentation ====="
bcftools query -l "$WIDER" > all.samples
printf '%s\n' Altai Vindija33.19 Denisova3 Denisova25 Chagyrskaya8 Mez1 > arc.want
grep -xFf arc.want all.samples > pig_arc.samples; grep -vxFf arc.want all.samples > pig_mod.samples
bcftools view -S pig_mod.samples "$WIDER" -Oz -o pig_mod.vcf.gz && bcftools index -f -t pig_mod.vcf.gz
bcftools view -S pig_arc.samples "$WIDER" -Oz -o pig_arc.vcf.gz && bcftools index -f -t pig_arc.vcf.gz
"$P2" --vcf pig_mod.vcf.gz "${IDARGS[@]}" --freq counts --pca 10 allele-wts --out pig_modpca
"$P2" --vcf pig_arc.vcf.gz "${IDARGS[@]}" --read-freq pig_modpca.acount \
  --score pig_modpca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-15 --out pig_arc_proj
coords_eigenvec pig_modpca.eigenvec "$OUT/pig_after_modern.tsv"
coords_sscore  pig_arc_proj.sscore  "$OUT/pig_after_archaic.tsv"
echo "  pig: $(($(wc -l < $OUT/pig_after_modern.tsv)-1)) modern + $(($(wc -l < $OUT/pig_after_archaic.tsv)-1)) archaic"

# ---- 2. WHOLE GENOME -----------------------------------------------------------------
echo "===== 2. whole genome ====="
"$P2" --bfile "$MM/sgdp166" --extract "$MM/sgdp166_prune.prune.in" --freq counts --pca 10 allele-wts --out wg_modpca
# archaic WG genotypes at the pruned SNPs: clean names, normalise contigs, extract
awk -F: '{print $1"\t"$2}' "$MM/sgdp166_prune.prune.in" | sort -k1,1 -k2,2n -u > prune.pos
bcftools query -l "$AWG" | awk 'BEGIN{
    m["den11"]="Denisova11";m["den5"]="Denisova5";m["goyet"]="Goyet";m["hst"]="Hohlenstein_Stadel";
    m["lesCottes"]="Les_Cottes";m["mez2"]="Mezmaiskaya2";m["Sid1253"]="El_Sidron";m["spy"]="Spy";m["vindija87"]="Vindija87"}
  {old=$0;n=$0;sub(/.*\//,"",n);sub(/_sorted\.bam$/,"",n);if(n in m)n=m[n];print old"\t"n}' > arc_rename.txt
CHRPREFIX=$(bcftools view -h "$AWG" | grep -qE '^##contig=<ID=chr' && echo yes || echo no)
if [ "$CHRPREFIX" = yes ]; then for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > chr2bare.txt; fi
( bcftools reheader -s arc_rename.txt "$AWG" | bcftools view -T prune.pos -Ou \
    | { [ "$CHRPREFIX" = yes ] && bcftools annotate --rename-chrs chr2bare.txt -Ou || cat; } \
    | bcftools norm -d none -Oz -o wg_arc.vcf.gz ) 2>wg_arc.err && bcftools index -f -t wg_arc.vcf.gz
echo "  archaic WG samples: $(bcftools query -l wg_arc.vcf.gz | tr '\n' ' ')"
"$P2" --vcf wg_arc.vcf.gz "${IDARGS[@]}" --read-freq wg_modpca.acount \
  --score wg_modpca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-15 --out wg_arc_proj
coords_eigenvec wg_modpca.eigenvec "$OUT/wg_after_modern.tsv"
coords_sscore  wg_arc_proj.sscore  "$OUT/wg_after_archaic.tsv"
echo "  WG: $(($(wc -l < $OUT/wg_after_modern.tsv)-1)) modern + $(($(wc -l < $OUT/wg_after_archaic.tsv)-1)) archaic"

# ---- 3. VALIDATION: score the moderns with their own weights, must reproduce the PCA -
echo "===== 3. validation (modern self-projection r vs eigenvec) ====="
"$P2" --bfile "$MM/sgdp166" --extract "$MM/sgdp166_prune.prune.in" --read-freq wg_modpca.acount \
  --score wg_modpca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-15 --out wg_mod_selfproj 2>/dev/null
paste <(coords_eigenvec wg_modpca.eigenvec /dev/stdout) <(coords_sscore wg_mod_selfproj.sscore /dev/stdout) \
  | awk 'NR>1{print $2, $5}' | awk '{n++;sx+=$1;sy+=$2;sxx+=$1*$1;syy+=$2*$2;sxy+=$1*$2}
    END{r=(n*sxy-sx*sy)/sqrt((n*sxx-sx*sx)*(n*syy-sy*sy)); printf "  modern PC1 self-projection r = %.4f (want |r|>0.99)\n", r}'
echo "  WG archaic PC1 (should lean toward NON-African / Neanderthal-admixture side):"
awk -F'\t' 'NR>1{print "   "$1": PC1="$2}' "$OUT/wg_after_archaic.tsv" | grep -iE "Chagyr|Denisova3|Denisova5|Denisova25|Vi33"

{ echo "# archaic projection (plink2 recipe), pig + WG"; date;
  echo "pig SNPs=$(($(wc -l < pig_modpca.eigenvec.allele)-1)); WG pruned SNPs=$(wc -l < $MM/sgdp166_prune.prune.in)"; } > "$OUT/projection_run.txt"
echo "===== projection done ====="; date