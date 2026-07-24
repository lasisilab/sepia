#!/bin/bash
# =====================================================================================
# SEPIA — B2a: the "Lily, but correct" pigmentation-panel PCA.
#
# WHAT: extract the 222 pigmentation-panel SNPs at their CORRECT hg19 coordinates from
# the hg19 data (159 SGDP moderns + 6 high-coverage archaic snpAD VCFs), then run a
# pigmentation-panel PCA. This directly reverses bug A2: Lily read the panel's *hg38*
# positions out of *hg19* data, so 221/222 SNPs were monomorphic and the PCA was rank-1.
# Reading the SAME panel at the correct hg19 coordinates restores the real variation.
#
# WHY hg19 (not hg38): the genotype is build-independent and the data is physically hg19,
# so we extract at pos_hg19 from data/panel_hg19_hg38_map.tsv (all 222 verified by rsID).
# hg38 is only the reporting/integration layer for the later gene-network work; it is not
# needed to replicate/repair Lily's figures. No liftOver, no re-alignment.
#
# METHOD: JOINT PCA of moderns + high-cov archaics (no projection) so we don't need the
# 'allele-wts' modifier the 2019 cluster plink2 lacks. High-cov archaics are near-complete
# at well-covered pigmentation SNPs, so they don't distort the axes. The projection-based
# version (Lily's exact method) follows once a newer plink2 lands (shared with Sanity 2).
#
# LAYOUT (repo-is-workdir): run from the SEPIA cluster checkout AFTER `git pull` (needs
# data/panel_hg19_hg38_map.tsv, committed 8dfabb5). Big intermediates -> gitignored
# scratch/; the small PCA coords -> output/ (tracked, synced).
#     cd /nfs/turbo/lsa-tlasisi1/sepia && git pull && sbatch code/pig_panel_extract.slurm
# =====================================================================================
#SBATCH --job-name=sepia_pig_panel
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_pig_panel_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
SRC=/nfs/turbo/lsa-tlasisi1/lheald_thesis/aDNA_data          # Lily's data — READ-ONLY seed
MAP="$REPO/data/panel_hg19_hg38_map.tsv"                     # rsID<->hg19<->hg38 map (222 SNPs)
BCFLIST="$REPO/scratch/modern_merge/sgdp_bcf.list"           # 159 validated per-sample SGDP BCFs
ARCVCF="$GENOMES/archaic/vcf"                                # high-cov archaic snpAD VCFs (per-chr)
SCR="$REPO/scratch/pig_panel"                                # intermediates (gitignored)
OUT="$REPO/output/pca_pig"                                   # small synced results (tracked)
mkdir -p "$SCR" "$OUT"; cd "$SCR"

echo "===== SEPIA B2a — pigmentation panel @ correct hg19 coords ====="; date
echo "map: $MAP"; echo "scratch: $SCR"; echo "output: $OUT"

module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21   >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
module load plink/2.0       >/dev/null 2>&1 || module load plink2   >/dev/null 2>&1 || true
PLINK2="$(command -v plink2 || command -v plink)"
echo "--- tools ---"; command -v bcftools; echo "plink2: $PLINK2"
[ -n "$PLINK2" ] || { echo "!! plink2 not found"; exit 2; }
[ -f "$MAP" ] || { echo "!! map not found ($MAP) — did you 'git pull'?"; exit 2; }
[ -f "$BCFLIST" ] || { echo "!! SGDP BCF list not found ($BCFLIST) — run wg_modern_merge first"; exit 2; }

# ---- 1. build hg19 targets from the map (chrom=$2, pos_hg19=$4) ----------------------
# bcftools -T wants a tab-delimited CHROM<TAB>POS file (1-based). Build BOTH a bare-contig
# and a chr-prefixed version, and pick per-source by the source's own contig style.
echo "===== 1. hg19 targets from the map ====="
awk -F'\t' 'NR>1 && $4!="" && $4!="NA"{print $2"\t"$4}' "$MAP" | sort -k1,1 -k2,2n -u > targets.bare.tsv
sed 's/^/chr/' targets.bare.tsv > targets.chr.tsv
# rsID/gene label per chr:pos (for downstream annotation)
awk -F'\t' 'NR>1 && $4!="" && $4!="NA"{print $2":"$4"\t"$1"\t"$9}' "$MAP" | sort -u > pos2rsid.tsv
NT=$(wc -l < targets.bare.tsv)
echo "  panel positions with hg19 coords: $NT"
CHRS=$(cut -f1 targets.bare.tsv | sort -u -V | tr '\n' ' ')
echo "  chromosomes spanned: $CHRS"
[ "$NT" -ge 200 ] || { echo "!! only $NT targets (<200); aborting"; exit 3; }

# helper: does this file use chr-prefixed contigs internally?  echoes "chr" or "bare"
detect_style () { bcftools view -h "$1" 2>/dev/null | grep -qE '^##contig=<ID=chr' && echo chr || echo bare; }
# helper: correct targets file for a source's style
targets_for () { [ "$(detect_style "$1")" = chr ] && echo targets.chr.tsv || echo targets.bare.tsv; }

# ---- 2. extract the panel from each SGDP sample (stream -T, parallel) ----------------
echo "===== 2. extract SGDP (159) at the panel hg19 positions ====="
mkdir -p sgdp; : > sgdp.merge.list
SG_STYLE=$(detect_style "$(head -1 "$BCFLIST")"); SG_T=$([ "$SG_STYLE" = chr ] && echo targets.chr.tsv || echo targets.bare.tsv)
echo "  SGDP contig style: $SG_STYLE  (targets: $SG_T)"
i=0
while read -r f; do
  b=$(basename "$f" .bcf)
  ( bcftools view -T "$SG_T" "$f" -Ob -o "sgdp/$b.bcf" 2>/dev/null && bcftools index -f "sgdp/$b.bcf" 2>/dev/null ) &
  i=$((i+1)); [ $((i % 16)) -eq 0 ] && wait
done < "$BCFLIST"
wait
for f in sgdp/*.bcf; do [ -f "$f.csi" ] && echo "$f" >> sgdp.merge.list; done
echo "  SGDP per-sample panel files: $(wc -l < sgdp.merge.list)"

# ---- 3. extract the panel from each high-cov archaic (per-chr snpAD VCFs) ------------
# Each genome is one directory of per-chromosome VCFs with genome-specific filenames.
# For each genome: for each panel chromosome, glob its chrN file, extract the panel
# positions, then concat -> one single-sample VCF renamed to a clean archaic label.
# NB: these per-chr files are ~4 GB all-sites VCFs, so we extract by INDEXED SEEK
# (`-R`, random access) not `-T` (which streams the whole file). Only ~10 positions per
# chromosome => seeking is instant; streaming 4 GB per file would take hours. (This is
# the opposite call from the SGDP step, where -T is right: there we hit 1.4M positions,
# so -R would mean 1.4M index seeks -- death on NFS. Few positions => seek; many => stream.)
echo "===== 3. extract high-cov archaics at the panel hg19 positions ====="
mkdir -p arc
# dir_name -> clean sample label
declare -A ARC=( [Altai]=Altai [Vindija33.19]=Vindija33.19 [Denisova]=Denisova3 \
                 [Denisova25]=Denisova25 [Chagyrskaya8]=Chagyrskaya8 [Mez1]=Mez1 )
: > arc.merge.list
for d in "${!ARC[@]}"; do
  gdir="$ARCVCF/$d"; label="${ARC[$d]}"
  [ -d "$gdir" ] || { echo "  skip $label — dir missing ($gdir)"; continue; }
  # pick a representative file to detect contig style
  rep=$(ls "$gdir"/chr*.vcf.gz 2>/dev/null | head -1)
  [ -n "$rep" ] || { echo "  skip $label — no chr VCFs in $gdir"; continue; }
  st=$(detect_style "$rep"); T=$([ "$st" = chr ] && echo targets.chr.tsv || echo targets.bare.tsv)
  : > "arc/$label.parts"
  for c in $CHRS; do
    # match chrN followed by . or _ (so chr1 != chr15); take the first match
    cf=$(ls "$gdir"/chr${c}[._]*.vcf.gz 2>/dev/null | grep -vE '\.(tbi|csi)$' | head -1)
    [ -n "$cf" ] || continue
    # per-chr targets for this source's style
    grep -P "^(chr)?${c}\t" "$T" > "arc/$label.chr${c}.tgt" 2>/dev/null || awk -v c="$c" -F'\t' '$1==c || $1=="chr"c' "$T" > "arc/$label.chr${c}.tgt"
    [ -s "arc/$label.chr${c}.tgt" ] || continue
    if bcftools view -R "arc/$label.chr${c}.tgt" "$cf" -Ob -o "arc/$label.chr${c}.bcf" 2>/dev/null; then
      bcftools index -f "arc/$label.chr${c}.bcf" 2>/dev/null
      echo "arc/$label.chr${c}.bcf" >> "arc/$label.parts"
    fi
  done
  if [ -s "arc/$label.parts" ]; then
    bcftools concat -f "arc/$label.parts" -Ob 2>/dev/null | bcftools sort -Ob -o "arc/$label.raw.bcf" 2>/dev/null
    # rename the single sample to the clean label, normalize contigs to bare to match SGDP
    printf '%s\n' "$(bcftools query -l arc/$label.raw.bcf | head -1) $label" > "arc/$label.rename"
    if [ "$st" = chr ] && [ "$SG_STYLE" = bare ]; then
      for c in $CHRS; do echo "chr$c $c"; done > arc/chr2bare.txt
      bcftools annotate --rename-chrs arc/chr2bare.txt "arc/$label.raw.bcf" -Ob 2>/dev/null \
        | bcftools reheader -s "arc/$label.rename" -o "arc/$label.bcf" 2>/dev/null
    else
      bcftools reheader -s "arc/$label.rename" "arc/$label.raw.bcf" -o "arc/$label.bcf" 2>/dev/null
    fi
    bcftools index -f "arc/$label.bcf" 2>/dev/null
    n=$(bcftools view -H "arc/$label.bcf" 2>/dev/null | wc -l)
    echo "  $label: $n panel records"
    [ -f "arc/$label.bcf.csi" ] && echo "arc/$label.bcf" >> arc.merge.list
  else
    echo "  $label: no records extracted"
  fi
done
echo "  archaic panel files: $(wc -l < arc.merge.list)"

# ---- 4. merge moderns + archaics at the panel loci ----------------------------------
echo "===== 4. merge -> pig_panel.vcf.gz ====="
cat sgdp.merge.list arc.merge.list > all.merge.list
# NB: no `-0` — a genotype absent in a sample must stay ./. (missing), never be invented as
# hom-ref. plink mean-imputes missing for the PCA; fabricating ref calls would bias archaics.
bcftools merge -l all.merge.list --threads 16 -Oz -o pig_panel.vcf.gz 2>merge.err || { echo "!! merge failed"; tail -5 merge.err; exit 4; }
bcftools index -t -f pig_panel.vcf.gz
NSAMP=$(bcftools query -l pig_panel.vcf.gz | wc -l); NSITE=$(bcftools index -n pig_panel.vcf.gz)
echo "  merged samples: $NSAMP   sites: $NSITE"

# ---- 4b. TRACER VALIDATION — the whole point of the fix --------------------------------
# Vindija 33.19 at the true hg19 SLC24A5 coord (15:48426484) must return a real genotype.
# (Under the A2 bug this locus was queried at hg38 15:48134288 and came back monomorphic.)
echo "===== 4b. tracer check: Vindija33.19 @ 15:48426484 (SLC24A5) ====="
bcftools view -H -r 15:48426484 pig_panel.vcf.gz 2>/dev/null | cut -f1-5,10- | head -1 \
  || bcftools view -H -r chr15:48426484 pig_panel.vcf.gz 2>/dev/null | cut -f1-5,10- | head -1
echo "  (expect a real REF/ALT + genotypes, not an empty/monomorphic line)"

# ---- 5. PLINK: joint PCA (moderns + high-cov archaics) + modern-only scree -----------
echo "===== 5. PCA on the 222-SNP panel ====="
"$PLINK2" --vcf pig_panel.vcf.gz \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
  --max-alleles 2 --snps-only --make-bed --out pig_panel
echo "  panel: $(wc -l < pig_panel.fam) samples, $(wc -l < pig_panel.bim) biallelic SNPs"
"$PLINK2" --bfile pig_panel --freq --out pig_panel_freq
# joint PCA (all samples). #PCs must be < #samples and <= #SNPs; 10 is safe here.
"$PLINK2" --bfile pig_panel --pca 10 --threads 8 --out pig_panel_pca
[ -f pig_panel_pca.eigenvec ] || { echo "!! joint PCA failed"; exit 6; }
echo "  --- joint eigenvalues (rank-1 => broken; several real PCs => fixed) ---"; cat pig_panel_pca.eigenval

# modern-only PCA for the direct before/after vs Lily's degenerate modern panel.
# Filter the .fam (cols FID IID ...) by dropping the 6 archaic labels; pass as --keep.
grep -vF -e Altai -e Vindija33.19 -e Denisova3 -e Denisova25 -e Chagyrskaya8 -e Mez1 pig_panel.fam > modern.keep || true
if [ -s modern.keep ]; then
  "$PLINK2" --bfile pig_panel --keep modern.keep --pca 10 --threads 8 --out pig_panel_modern_pca 2>/dev/null \
    && { echo "  --- modern-only eigenvalues ($(wc -l < modern.keep) samples) ---"; cat pig_panel_modern_pca.eigenval; } \
    || echo "  (modern-only PCA skipped)"
fi

# ---- 6. place small outputs into the tracked repo -----------------------------------
echo "===== 6. copy coords -> output/pca_pig ====="
cp -f pig_panel_pca.eigenvec pig_panel_pca.eigenval pig_panel_freq.afreq "$OUT/" 2>/dev/null || true
[ -f pig_panel_modern_pca.eigenvec ] && cp -f pig_panel_modern_pca.eigenvec pig_panel_modern_pca.eigenval "$OUT/" 2>/dev/null || true
cp -f pos2rsid.tsv "$OUT/panel_pos2rsid.tsv" 2>/dev/null || true
# a small provenance note
{ echo "# B2a pigmentation-panel PCA — extracted at CORRECT hg19 coords from the panel map"; date;
  echo "samples=$NSAMP  panel_sites_merged=$NSITE  biallelic_snps_in_pca=$(wc -l < pig_panel.bim)";
  echo "archaics: $(cat arc.merge.list | wc -l) high-cov  |  moderns: $(cat sgdp.merge.list | wc -l) SGDP"; } > "$OUT/B2a_run.txt"
echo "  -> $OUT (eigenvec/eigenval/freq + provenance)"

# ---- 7. auto-sync results to GitHub so they survive an SSH-tunnel drop ---------------
# The SLURM job runs independently of the interactive SSH ControlMaster; this pushes the
# small result files straight to GitHub (BatchMode key, no agent) so they land in the repo
# even if Tina's tunnel is gone -- the figure + redeploy then happen laptop-side, no cluster.
echo "===== 7. auto-commit + push results to GitHub ====="
cd "$REPO" || exit 0
export GIT_SSH_COMMAND="ssh -o BatchMode=yes"
git pull --ff-only origin main >/dev/null 2>&1 || echo "  (pull skipped)"
if [ -n "$(git status --porcelain output/pca_pig 2>/dev/null)" ]; then
  git add output/pca_pig
  git -c user.name="SEPIA cluster" -c user.email="tlasisi@umich.edu" \
      commit -q -m "results(B2a): pigmentation-panel PCA output (auto-push from job ${SLURM_JOB_ID:-NA})"
  git push origin main 2>&1 | tail -2 && echo "  pushed to GitHub." || echo "  PUSH FAILED -- pull from cluster when reconnected."
else
  echo "  no output changes to push."
fi
echo "===== B2a done ====="; date
echo "NEXT (laptop, no cluster needed): git pull; Rscript plan/figures/make_pig_before_after.R; redeploy."
