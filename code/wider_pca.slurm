#!/bin/bash
# =====================================================================================
# SEPIA B3 — the WIDER pigmentation PCA: the extended ~1,155-SNP set (not just Lily's 222).
#
# Extracts every resolvable SNP in data/pigment_snp_map.tsv (status==ok, at its correct hg19
# coordinate) from the hg19 data (159 SGDP + high-cov archaics), then runs a JOINT PCA -- twice:
#   (a) ALL ok SNPs        -> output/pca_wider/pig_wider_pca.*
#   (b) high-confidence ONLY (confidence==high) -> output/pca_wider/pig_wider_highconf_pca.*
# so we can see how much the unreplicated single-GWAS tail moves the result.
#
# Same extraction machinery as pig_panel_extract.slurm (SGDP: stream -T; archaics: indexed -R
# seek on the 4 GB per-chr all-sites VCFs). Results are pushed from the LOGIN node afterwards
# (compute nodes have no SSH agent), so this job does NOT auto-commit.
#     cd /nfs/turbo/lsa-tlasisi1/sepia && git pull && sbatch code/wider_pca.slurm
# =====================================================================================
#SBATCH --job-name=sepia_wider_pca
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_wider_pca_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
MAP="$REPO/data/pigment_snp_map.tsv"                        # cols: 1 rsid .. 3 chrom .. 5 pos_hg19 .. 13 confidence .. 20 status
BCFLIST="$REPO/scratch/modern_merge/sgdp_bcf.list"          # 159 validated per-sample SGDP BCFs
ARCVCF="$GENOMES/archaic/vcf"
SCR="$REPO/scratch/wider_pca"; OUT="$REPO/output/pca_wider"
mkdir -p "$SCR" "$OUT"; cd "$SCR"

echo "===== SEPIA B3 — wider pigmentation PCA (extended SNP set) ====="; date
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21   >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
module load plink/2.0       >/dev/null 2>&1 || module load plink2   >/dev/null 2>&1 || true
PLINK2="$(command -v plink2 || command -v plink)"
echo "bcftools: $(command -v bcftools) | plink2: $PLINK2"
[ -n "$PLINK2" ] && [ -f "$MAP" ] && [ -f "$BCFLIST" ] || { echo "!! missing plink2/map/bcflist"; exit 2; }

# ---- 1. targets from the extended map (ok rows) + a high-confidence position list ----
echo "===== 1. targets from $MAP ====="
awk -F'\t' 'NR>1 && $20=="ok" && $5!="" && $5!="NA"{print $3"\t"$5}' "$MAP" | sort -k1,1 -k2,2n -u > targets.bare.tsv
sed 's/^/chr/' targets.bare.tsv > targets.chr.tsv
awk -F'\t' 'NR>1 && $20=="ok" && $13=="high" && $5!="" {print $3":"$5}' "$MAP" | sort -u > highconf.chrpos
NT=$(wc -l < targets.bare.tsv); NH=$(wc -l < highconf.chrpos)
echo "  target SNPs: $NT  (high-confidence: $NH)"
CHRS=$(cut -f1 targets.bare.tsv | sort -u -V | tr '\n' ' ')
[ "$NT" -ge 500 ] || { echo "!! only $NT targets; aborting"; exit 3; }

detect_style(){ bcftools view -h "$1" 2>/dev/null | grep -qE '^##contig=<ID=chr' && echo chr || echo bare; }

# ---- 2. SGDP extraction (stream -T; parallel) ---------------------------------------
echo "===== 2. extract 159 SGDP ====="
mkdir -p sgdp; : > sgdp.merge.list
SG_T=$([ "$(detect_style "$(head -1 "$BCFLIST")")" = chr ] && echo targets.chr.tsv || echo targets.bare.tsv)
i=0
while read -r f; do
  b=$(basename "$f" .bcf)
  ( bcftools view -T "$SG_T" "$f" -Ob -o "sgdp/$b.bcf" 2>/dev/null && bcftools index -f "sgdp/$b.bcf" 2>/dev/null ) &
  i=$((i+1)); [ $((i % 16)) -eq 0 ] && wait
done < "$BCFLIST"; wait
for f in sgdp/*.bcf; do [ -f "$f.csi" ] && echo "$f" >> sgdp.merge.list; done
echo "  SGDP panel files: $(wc -l < sgdp.merge.list)"

# ---- 3. high-cov archaic extraction (indexed -R seek) -------------------------------
echo "===== 3. extract high-cov archaics (-R seek) ====="
mkdir -p arc; : > arc.merge.list
declare -A ARC=( [Altai]=Altai [Vindija33.19]=Vindija33.19 [Denisova]=Denisova3 \
                 [Denisova25]=Denisova25 [Chagyrskaya8]=Chagyrskaya8 [Mez1]=Mez1 )
SG_STYLE=$(detect_style "$(head -1 "$BCFLIST")")
for d in "${!ARC[@]}"; do
  gdir="$ARCVCF/$d"; label="${ARC[$d]}"; [ -d "$gdir" ] || continue
  rep=$(ls "$gdir"/chr*.vcf.gz 2>/dev/null | head -1); [ -n "$rep" ] || continue
  st=$(detect_style "$rep"); T=$([ "$st" = chr ] && echo targets.chr.tsv || echo targets.bare.tsv)
  : > "arc/$label.parts"
  for c in $CHRS; do
    cf=$(ls "$gdir"/chr${c}[._]*.vcf.gz 2>/dev/null | grep -vE '\.(tbi|csi)$' | head -1); [ -n "$cf" ] || continue
    awk -v c="$c" -F'\t' '$1==c || $1=="chr"c' "$T" > "arc/$label.chr${c}.tgt"
    [ -s "arc/$label.chr${c}.tgt" ] || continue
    if bcftools view -R "arc/$label.chr${c}.tgt" "$cf" -Ob -o "arc/$label.chr${c}.bcf" 2>/dev/null; then
      bcftools index -f "arc/$label.chr${c}.bcf" 2>/dev/null; echo "arc/$label.chr${c}.bcf" >> "arc/$label.parts"
    fi
  done
  if [ -s "arc/$label.parts" ]; then
    bcftools concat -f "arc/$label.parts" -Ob 2>/dev/null | bcftools sort -Ob -o "arc/$label.raw.bcf" 2>/dev/null
    printf '%s %s\n' "$(bcftools query -l arc/$label.raw.bcf | head -1)" "$label" > "arc/$label.rename"
    if [ "$st" = chr ] && [ "$SG_STYLE" = bare ]; then
      for c in $CHRS; do echo "chr$c $c"; done > arc/chr2bare.txt
      bcftools annotate --rename-chrs arc/chr2bare.txt "arc/$label.raw.bcf" -Ob 2>/dev/null | bcftools reheader -s "arc/$label.rename" -o "arc/$label.bcf" 2>/dev/null
    else
      bcftools reheader -s "arc/$label.rename" "arc/$label.raw.bcf" -o "arc/$label.bcf" 2>/dev/null
    fi
    bcftools index -f "arc/$label.bcf" 2>/dev/null
    echo "  $label: $(bcftools view -H "arc/$label.bcf" 2>/dev/null | wc -l) records"
    [ -f "arc/$label.bcf.csi" ] && echo "arc/$label.bcf" >> arc.merge.list
  else echo "  $label: no records extracted"; fi
done
echo "  archaic files: $(wc -l < arc.merge.list)"

# ---- 4. merge + PLINK bed -----------------------------------------------------------
echo "===== 4. merge + bed ====="
cat sgdp.merge.list arc.merge.list > all.merge.list
bcftools merge -l all.merge.list --threads 16 -Oz -o wider.vcf.gz 2>merge.err || { echo "!! merge failed"; tail -3 merge.err; exit 4; }
bcftools index -t -f wider.vcf.gz
echo "  merged: $(bcftools query -l wider.vcf.gz | wc -l) samples, $(bcftools index -n wider.vcf.gz) sites"
"$PLINK2" --vcf wider.vcf.gz --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
  --max-alleles 2 --snps-only --make-bed --out pig_wider
echo "  bed: $(wc -l < pig_wider.fam) samples, $(wc -l < pig_wider.bim) biallelic SNPs"

# ---- 5. PCA: all SNPs, then high-confidence only ------------------------------------
echo "===== 5. PCA (all) ====="
"$PLINK2" --bfile pig_wider --freq --out pig_wider_freq
"$PLINK2" --bfile pig_wider --pca 10 --threads 8 --out pig_wider_pca
echo "  all-SNP eigenvalues:"; cat pig_wider_pca.eigenval 2>/dev/null | tr '\n' ' '; echo
# high-confidence subset: match bim variant IDs (chr:pos:..) by chr:pos
awk 'FNR==NR{k[$1]=1;next} (($1":"$4) in k){print $2}' highconf.chrpos pig_wider.bim > highconf.snplist
echo "===== 5b. PCA (high-confidence only: $(wc -l < highconf.snplist) SNPs) ====="
"$PLINK2" --bfile pig_wider --extract highconf.snplist --pca 10 --threads 8 --out pig_wider_highconf_pca 2>/dev/null \
  && { echo "  high-conf eigenvalues:"; cat pig_wider_highconf_pca.eigenval | tr '\n' ' '; echo; } || echo "  (high-conf PCA skipped)"

# ---- 6. place outputs ---------------------------------------------------------------
cp -f pig_wider_pca.eigenvec pig_wider_pca.eigenval pig_wider_freq.afreq "$OUT/" 2>/dev/null || true
[ -f pig_wider_highconf_pca.eigenvec ] && cp -f pig_wider_highconf_pca.eigenvec pig_wider_highconf_pca.eigenval "$OUT/" 2>/dev/null
{ echo "# B3 wider pigmentation PCA @ correct hg19 coords"; date;
  echo "target_snps=$NT  high_conf=$NH  biallelic_in_pca=$(wc -l < pig_wider.bim)";
  echo "samples=$(wc -l < pig_wider.fam)  archaics=$(wc -l < arc.merge.list)"; } > "$OUT/B3_run.txt"
echo "===== B3 done ====="; date
echo "NEXT (login node): git pull --rebase && git add output/pca_wider && commit && push; then figure laptop-side."
