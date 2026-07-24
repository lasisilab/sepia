#!/bin/bash
# =====================================================================================
# SEPIA — archaic PROJECTION onto the modern PCA (pigmentation + whole genome), plink2.
#
# Key robustness choices (learned from two failed attempts):
#  - build ONE combined plink dataset (moderns + archaics) per analysis, so variant IDs and
#    allele coding are identical -> --score matches. (Separate modern/archaic files fail:
#    the archaic WG calls carry ALT='.', so their chr:pos:ref:alt IDs never match the moderns'.)
#  - read the pigmentation VCF with plink2 directly (bcftools rejects a malformed GQ field).
#  - detect the A1 and PC columns of .eigenvec.allele by header name (a PROVISIONAL_REF?
#    column shifts them), instead of hard-coding.
# Recipe: modern subset --freq counts + --pca allele-wts; archaic subset --read-freq (modern)
#   --score <allele-wts> <idcol> <a1col> header-read no-mean-imputation variance-standardize.
# VALIDATION: modern self-projection must reproduce the eigenvec (r>0.99).
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
P2="$REPO/tools/plink2"
SCR="$REPO/scratch/project"; OUT="$REPO/output/projections"
mkdir -p "$SCR" "$OUT"; cd "$SCR"
echo "===== SEPIA archaic projection (combined-dataset plink2 recipe) ====="; date
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21 >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
"$P2" --version || { echo "!! plink2 missing"; exit 2; }
IDARGS=(--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing --max-alleles 2 --snps-only)

colnum(){ head -1 "$1" | tr '\t' '\n' | grep -nxF "$2" | cut -d: -f1; }   # 1-based col of header name
coords_eigenvec(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1")A=i;if($i=="PC2")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }
coords_sscore(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1_AVG")A=i;if($i=="PC2_AVG")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }

# project(): $1=combined bed prefix -> writes $OUT/${2}_after_{modern,archaic}.tsv  ($2 = pig|wg)
# moderns = IID starts with 'LP' (SGDP library IDs); archaics = everything else.
project(){
  local bed="$1" tag="$2"
  awk '$2 ~ /^LP/{print $1"\t"$2}' "$bed.fam" > "${tag}_mod.keep"
  awk '$2 !~ /^LP/{print $1"\t"$2}' "$bed.fam" > "${tag}_arc.keep"
  echo "  $tag: $(wc -l < ${tag}_mod.keep) modern, $(wc -l < ${tag}_arc.keep) archaic"
  "$P2" --bfile "$bed" --keep "${tag}_mod.keep" --freq counts --pca 10 allele-wts --out "${tag}_modpca"
  local A1=$(colnum "${tag}_modpca.eigenvec.allele" A1)
  local P1=$(colnum "${tag}_modpca.eigenvec.allele" PC1); local P10=$(colnum "${tag}_modpca.eigenvec.allele" PC10)
  echo "  $tag eigenvec.allele: A1 col=$A1, PC cols=${P1}-${P10}"
  "$P2" --bfile "$bed" --keep "${tag}_arc.keep" --read-freq "${tag}_modpca.acount" \
    --score "${tag}_modpca.eigenvec.allele" 2 "$A1" header-read no-mean-imputation variance-standardize \
    --score-col-nums "${P1}-${P10}" --out "${tag}_arc_proj"
  coords_eigenvec "${tag}_modpca.eigenvec" "$OUT/${tag}_after_modern.tsv"
  coords_sscore  "${tag}_arc_proj.sscore"  "$OUT/${tag}_after_archaic.tsv"
  # validation: self-project moderns, correlate with eigenvec PC1
  "$P2" --bfile "$bed" --keep "${tag}_mod.keep" --read-freq "${tag}_modpca.acount" \
    --score "${tag}_modpca.eigenvec.allele" 2 "$A1" header-read no-mean-imputation variance-standardize \
    --score-col-nums "${P1}-${P10}" --out "${tag}_mod_self" 2>/dev/null
  coords_eigenvec "${tag}_modpca.eigenvec" ev.tmp; coords_sscore "${tag}_mod_self.sscore" sp.tmp
  join -t$'\t' <(sort ev.tmp) <(sort sp.tmp) 2>/dev/null | awk 'NR>0{print $2,$4}' \
    | awk '{n++;sx+=$1;sy+=$2;sxx+=$1*$1;syy+=$2*$2;sxy+=$1*$2}
      END{if(n>2){r=(n*sxy-sx*sy)/sqrt((n*sxx-sx*sx)*(n*syy-sy*sy));printf "  %s VALIDATION: modern self-proj r=%.4f (want |r|>0.99)\n",t,r}}' t="$tag"
}

# ---- 1. PIGMENTATION: plink2 reads wider.vcf.gz directly (bcftools rejects its GQ) ----
echo "===== 1. pigmentation ====="
"$P2" --vcf "$WIDER" "${IDARGS[@]}" --make-bed --out pig_all
project pig_all pig

# ---- 2. WHOLE GENOME: build combined moderns+archaics bed at the pruned SNPs ----------
echo "===== 2. whole genome ====="
"$P2" --bfile "$MM/sgdp166" --extract "$MM/sgdp166_prune.prune.in" --export vcf bgz --out wg_mod
bcftools index -f -t wg_mod.vcf.gz
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
echo "  moderns VCF: $(bcftools index -n wg_mod.vcf.gz) sites; archaic VCF: $(bcftools index -n wg_arc.vcf.gz) sites"
# merge (archaic ALT='.' folds into the modern REF/ALT), then one bed
bcftools merge wg_mod.vcf.gz wg_arc.vcf.gz -Oz -o wg_comb.vcf.gz 2>wg_merge.err || { echo "!! merge failed"; tail -3 wg_merge.err; exit 4; }
"$P2" --vcf wg_comb.vcf.gz "${IDARGS[@]}" --make-bed --out wg_all
project wg_all wg

echo "  WG archaic PC1 (expect NON-African / Neanderthal-admixture side):"
awk -F'\t' 'NR>1{print "   "$1": PC1="$2}' "$OUT/wg_after_archaic.tsv" | grep -iE "Chagyr|Denisova3|Denisova5|Denisova25|Vi33"
{ echo "# archaic projection (combined-dataset plink2 recipe)"; date; } > "$OUT/projection_run.txt"
echo "===== projection done ====="; date