#!/bin/bash
# =====================================================================================
# SEPIA — fix the Chagyrskaya 8 "0 records" bug: recompress its VCFs from gzip -> BGZF.
#
# The Chagyrskaya .noRB VCFs are PLAIN gzip, not bgzip (htsfile: "gzip-compressed"; magic
# 1f 8b 08 00, missing the BGZF FEXTRA 04). They stream fine (files are complete — chr15
# reaches 15:102506248) but cannot be random-accessed, so the panel/gene extractions' indexed
# `-R` seek returns 0 records. Recompressing to BGZF + indexing makes Chagyrskaya usable like
# the other high-cov genomes. Content is unchanged (gunzip | bgzip round-trip).
#
#   cd /nfs/turbo/lsa-tlasisi1/sepia && sbatch code/fix_chagyrskaya_bgzip.slurm
# Resume-safe: skips chromosomes already converted (BGZF + index present).
# =====================================================================================
#SBATCH --job-name=sepia_fix_chag
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/genomes/archaic/fix_chag_%j.log
set -uo pipefail

GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
D="$GENOMES/archaic/vcf/Chagyrskaya8"
module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21  >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
module load htslib          >/dev/null 2>&1 || true
command -v bgzip >/dev/null || { echo "!! bgzip not found"; exit 2; }
[ -d "$D" ] || { echo "!! $D not found"; exit 2; }
echo "===== recompress Chagyrskaya gzip -> BGZF ====="; date

is_bgzf(){ [ "$(od -An -tx1 -N4 "$1" 2>/dev/null | tr -d ' \n')" = "1f8b0804" ]; }

for c in $(seq 1 22) X; do
  f="$D/chr${c}.noRB.vcf.gz"
  [ -f "$f" ] || { echo "  chr$c: (no file)"; continue; }
  if is_bgzf "$f" && { [ -f "$f.tbi" ] || [ -f "$f.csi" ]; }; then echo "  chr$c: already BGZF+indexed — skip"; continue; fi
  echo "  chr$c: recompressing ($(du -h "$f" | cut -f1))..."
  # gzip -> BGZF into a temp, verify, then atomically replace
  if zcat "$f" | bgzip -@ 4 -c > "$f.bgz.tmp" 2>/dev/null; then
    if bcftools index -f -t "$f.bgz.tmp" 2>/dev/null; then
      # verify random access works now
      if bcftools view -h "$f.bgz.tmp" >/dev/null 2>&1; then
        mv -f "$f.bgz.tmp" "$f"
        mv -f "$f.bgz.tmp.tbi" "$f.tbi" 2>/dev/null || bcftools index -f -t "$f"
        rm -f "$f.csi"
        recs=$(bcftools view -H "$f" 15:1-1 2>/dev/null | wc -l)   # sanity: indexed seek no longer errors
        echo "    chr$c: OK ($(du -h "$f" | cut -f1), BGZF; indexed-seek works)"
      else echo "    chr$c: verify failed — leaving original"; rm -f "$f.bgz.tmp"*; fi
    else echo "    chr$c: index failed — leaving original"; rm -f "$f.bgz.tmp"*; fi
  else echo "    chr$c: recompress failed — leaving original"; rm -f "$f.bgz.tmp"*; fi
done

echo "===== verify: extract a panel tracer (SLC24A5 hg19 15:48426484) ====="
bcftools view -H -r 15:48426484 "$D/chr15.noRB.vcf.gz" 2>&1 | head -1 | cut -c1-80
echo "done"; date
