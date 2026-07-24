#!/bin/bash
# =====================================================================================
# SEPIA B6/B7 — pigmentation GENE-REGION PCA: place archaics using the variation across
# whole pigmentation-gene regions, NOT a specific ascertained SNP panel.
#
# Regions: the 1,004-gene union (data/pigment_gene_map.tsv), gene body +/-10 kb, merged.
# Variants: all SGDP variants in those regions (from the 159-sample merged set), LD-pruned.
# Method: same validated recipe as project_archaics.slurm — build the modern PCA
#   (--pca allele-wts), project the archaics (--score, read-freq + variance-standardize +
#   no-mean-imputation), on ONE combined moderns+archaics bed. Modern self-projection validates r~1.
#     cd /nfs/turbo/lsa-tlasisi1/sepia && git pull && sbatch code/gene_region_pca.slurm
# =====================================================================================
#SBATCH --job-name=sepia_generegion
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_generegion_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
SRC=/nfs/turbo/lsa-tlasisi1/lheald_thesis/aDNA_data
MM="$REPO/scratch/modern_merge"; GMAP="$REPO/data/pigment_gene_map.tsv"
AWG="$SRC/ancient_wholegenome/ancient_sgdp_wg.vcf.gz"; P2="$REPO/tools/plink2"
SCR="$REPO/scratch/generegion"; OUT="$REPO/output/pca_generegion"
mkdir -p "$SCR" "$OUT"; cd "$SCR"
echo "===== SEPIA gene-region PCA ====="; date
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21 >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
"$P2" --version || { echo "!! plink2 missing"; exit 2; }
IDARGS=(--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing --max-alleles 2 --snps-only)
colnum(){ head -1 "$1" | tr '\t' '\n' | grep -nxF "$2" | cut -d: -f1; }
coords_eigenvec(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1")A=i;if($i=="PC2")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }
coords_sscore(){ awk 'NR==1{for(i=1;i<=NF;i++){if($i=="IID")I=i;if($i=="PC1_AVG")A=i;if($i=="PC2_AVG")B=i};print "IID\tPC1\tPC2";next}{print $I"\t"$A"\t"$B}' "$1" > "$2"; }

# ---- 1. gene regions BED (hg19, gene body +/-10kb, merged) ---------------------------
echo "===== 1. build gene-region BED (+/-10kb, merged) ====="
awk -F'\t' 'NR>1 && ($14=="ok"||$14=="hg19_only") && $6!="" && $3!="" {
    s=$6-10000; if(s<0)s=0; print $3"\t"s"\t"$7+10000 }' "$GMAP" \
  | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}' > gene_regions.bed
echo "  merged intervals: $(wc -l < gene_regions.bed)   span: $(awk '{n+=$3-$2}END{printf "%.1f Mb",n/1e6}' gene_regions.bed)"

# ---- 2. modern gene-region variants (from the 159-sample bed), LD-pruned -------------
echo "===== 2. modern gene-region variants + prune ====="
"$P2" --bfile "$MM/sgdp166" --extract bed1 gene_regions.bed --make-bed --out gr_mod_all
echo "  gene-region variants in moderns: $(wc -l < gr_mod_all.bim)"
"$P2" --bfile gr_mod_all --indep-pairwise 200 25 0.4 --out gr_prune
"$P2" --bfile gr_mod_all --extract gr_prune.prune.in --recode vcf bgz --out gr_mod
bcftools index -f -t gr_mod.vcf.gz
echo "  pruned gene-region variants: $(wc -l < gr_prune.prune.in)"

# ---- 3. archaic gene-region variants; combine; project -------------------------------
echo "===== 3. archaic extraction + combine + project ====="
bcftools query -l "$AWG" | awk 'BEGIN{
    m["den11"]="Denisova11";m["den5"]="Denisova5";m["goyet"]="Goyet";m["hst"]="Hohlenstein_Stadel";
    m["lesCottes"]="Les_Cottes";m["mez2"]="Mezmaiskaya2";m["Sid1253"]="El_Sidron";m["spy"]="Spy";m["vindija87"]="Vindija87"}
  {old=$0;n=$0;sub(/.*\//,"",n);sub(/_sorted\.bam$/,"",n);if(n in m)n=m[n];print old"\t"n}' > arc_rename.txt
CHRPREFIX=$(bcftools view -h "$AWG" | grep -qE '^##contig=<ID=chr' && echo yes || echo no)
if [ "$CHRPREFIX" = yes ]; then for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > chr2bare.txt; fi
# -T (streaming targets, no index) not -R (needs an indexed file; fails on the reheader pipe)
( bcftools reheader -s arc_rename.txt "$AWG" | bcftools view -T gene_regions.bed -Ou \
    | { [ "$CHRPREFIX" = yes ] && bcftools annotate --rename-chrs chr2bare.txt -Ou || cat; } \
    | bcftools norm -d exact -Oz -o gr_arc.vcf.gz ) 2>gr_arc.err && bcftools index -f -t gr_arc.vcf.gz
bcftools merge gr_mod.vcf.gz gr_arc.vcf.gz -Oz -o gr_comb.vcf.gz 2>gr_merge.err || { echo "!! merge failed"; tail -3 gr_merge.err; exit 4; }
# --rm-dup: drop any duplicate-ID variants (archaic ALT='.' records collapse to chr:pos:ref:. dups)
"$P2" --vcf gr_comb.vcf.gz "${IDARGS[@]}" --rm-dup force-first --make-bed --out gr_all

ARC="Chagyrskaya8 Denisova3 Denisova11 Denisova25 Denisova5 Goyet Hohlenstein_Stadel Les_Cottes L4741 Mezmaiskaya2 SP4903 El_Sidron Spy Vi33.19 Vindija87"
awk -v L="$ARC" 'BEGIN{n=split(L,a," ");for(i=1;i<=n;i++)A[a[i]]=1} ($2 in A){print $1"\t"$2 > "gr_arc.keep"} !($2 in A){print $1"\t"$2 > "gr_mod.keep"}' gr_all.fam
echo "  combined: $(wc -l < gr_mod.keep) modern, $(wc -l < gr_arc.keep) archaic"
"$P2" --bfile gr_all --keep gr_mod.keep --mac 1 --freq counts --pca 10 allele-wts --out gr_modpca
A1=$(colnum gr_modpca.eigenvec.allele A1); P1=$(colnum gr_modpca.eigenvec.allele PC1); P10=$(colnum gr_modpca.eigenvec.allele PC10)
"$P2" --bfile gr_all --keep gr_arc.keep --read-freq gr_modpca.acount \
  --score gr_modpca.eigenvec.allele 2 "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${P1}-${P10}" --out gr_arc_proj
# moderns from self-projection (sign-consistent with archaics) + validation
"$P2" --bfile gr_all --keep gr_mod.keep --read-freq gr_modpca.acount \
  --score gr_modpca.eigenvec.allele 2 "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${P1}-${P10}" --out gr_mod_self 2>/dev/null
coords_sscore gr_mod_self.sscore "$OUT/generegion_after_modern.tsv"
coords_sscore gr_arc_proj.sscore "$OUT/generegion_after_archaic.tsv"
coords_eigenvec gr_modpca.eigenvec ev.tmp
join -t$'\t' <(sort ev.tmp) <(sort "$OUT/generegion_after_modern.tsv") 2>/dev/null | awk 'NR>0{print $2,$4}' \
  | awk '{n++;sx+=$1;sy+=$2;sxx+=$1*$1;syy+=$2*$2;sxy+=$1*$2}END{if(n>2)printf "  VALIDATION modern self-proj r=%.4f (want |r|>0.99)\n",(n*sxy-sx*sy)/sqrt((n*sxx-sx*sx)*(n*syy-sy*sy))}'
echo "  hi-cov archaic PC1:"; awk -F'\t' 'NR>1{print "   "$1": "$2}' "$OUT/generegion_after_archaic.tsv" | grep -iE "Chagyr|Denisova3|Denisova5|Denisova25|Vi33"
{ echo "# gene-region PCA (1004-gene union, +/-10kb)"; date;
  echo "merged intervals=$(wc -l < gene_regions.bed); gene-region variants=$(wc -l < gr_mod_all.bim); pruned=$(wc -l < gr_prune.prune.in)"; } > "$OUT/generegion_run.txt"
echo "===== gene-region PCA done ====="; date