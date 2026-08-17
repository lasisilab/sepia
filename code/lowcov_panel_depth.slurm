#!/bin/bash
# =====================================================================================
# SEPIA — panel-locus read depth for the 3 LOW-coverage archaics (BAM-only samples).
#
# The high-cov archaics carry DP in their snpAD VCFs (already in output/panel_matrix/), but
# Denisova 11, Goyet Q56-1 and Les Cottés Z4-1514 are BAM-only. Here we stream samtools depth
# over each sample's run BAMs at the 222 panel hg19 positions and SUM the runs into one
# per-sample depth column, so these three extend the panel depth heatmap.
#
# Output: output/panel_matrix/lowcov_panel_depth.tsv  (chrom  pos  Denisova11  Goyet  LesCottes)
# Run from the cluster checkout:  sbatch code/lowcov_panel_depth.slurm
# =====================================================================================
#SBATCH --job-name=sepia_lowcov_depth
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_lowcov_depth_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
BAMROOT="$GENOMES/archaic/bam"
MAP="$REPO/data/panel_hg19_hg38_map.tsv"
OUT="$REPO/output/panel_matrix"; SCR="$REPO/scratch/lowcov_depth"
mkdir -p "$OUT" "$SCR"; cd "$SCR"

module load Bioinformatics >/dev/null 2>&1 || true
module load samtools       >/dev/null 2>&1 || true
command -v samtools >/dev/null || { echo "!! samtools not found"; exit 2; }
[ -f "$MAP" ] || { echo "!! map not found: $MAP"; exit 2; }

echo "===== SEPIA low-cov archaic panel depth ====="; date

# ---- 1. panel hg19 BED (0-based), both contig styles (map: col2=chrom, col4=pos_hg19) ----
awk -F'\t' 'NR>1 && $4!="" && $4!="NA"{print $2"\t"($4-1)"\t"$4}' "$MAP" | sort -k1,1 -k2,2n -u > panel.bare.bed
sed 's/^/chr/' panel.bare.bed > panel.chr.bed
echo "panel positions: $(wc -l < panel.bare.bed)"

# ---- 2. sample -> run BAMs ----
DEN="Denisova11/R5507_allruns.bam Denisova11/R5509_allruns.bam Denisova11/R5782.bam Denisova11/R9873_allruns.bam Denisova11/R9880_allruns.bam Denisova11/R9881_allruns.bam"
GOY="GoyetQ56-1/A9229_final_sorted.uniq.L35MQ0.bam"
LES="LesCottes_Z4-1514/A9230_final_sorted.uniq.L35MQ0.bam LesCottes_Z4-1514/A9291_final_sorted.uniq.L35MQ0.bam LesCottes_Z4-1514/A9350_final_sorted.uniq.L35MQ0.bam LesCottes_Z4-1514/A9394_final_sorted.uniq.L35MQ0.bam LesCottes_Z4-1514/A9420_final_sorted.uniq.L35MQ0.bam"

depth_for () {   # $1=label  $2..=relative BAM paths -> writes dep_<label>.txt (chrom pos summed_depth)
  local label="$1"; shift
  local bams=() first=""
  for b in "$@"; do
    local p="$BAMROOT/$b"
    if [ -f "$p" ]; then bams+=("$p"); [ -z "$first" ] && first="$p"; else echo "  MISSING $p"; fi
  done
  [ "${#bams[@]}" -gt 0 ] || { echo "  $label: no BAMs"; return; }
  # contig style from the first BAM's header
  local style bed
  style=$(samtools view -H "$first" | awk '/^@SQ/{print ($0 ~ /SN:chr/)?"chr":"bare"; exit}')
  bed=$([ "$style" = chr ] && echo panel.chr.bed || echo panel.bare.bed)
  echo "  $label: ${#bams[@]} BAM(s), contig style=$style"
  # -a = report every panel position (0 if uncovered); sum all run columns into one
  samtools depth -a -b "$bed" "${bams[@]}" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"}{s=0; for(i=3;i<=NF;i++) s+=$i; sub(/^chr/,"",$1); print $1,$2,s}' \
    | sort -k1,1 -k2,2n > "dep_$label.txt"
  echo "  $label: $(wc -l < dep_$label.txt) positions, mean depth $(awk '{n++;t+=$3}END{printf "%.2f", (n?t/n:0)}' dep_$label.txt)"
}

echo "===== 3. per-sample panel depth ====="
depth_for Denisova11 $DEN
depth_for Goyet      $GOY
depth_for LesCottes  $LES

# ---- 4. merge the three onto the panel loci (keyed by chrom:pos) ----
echo "===== 4. merge -> lowcov_panel_depth.tsv ====="
python3 - "$OUT" <<'PY'
import sys, os
OUT=sys.argv[1]; scr=os.getcwd()
def load(f):
    d={}
    if os.path.exists(f):
        for ln in open(f):
            c,p,v=ln.split()[:3]; d[(c,p)]=v
    return d
loci=[]
for ln in open("panel.bare.bed"):
    c,s,e=ln.split()[:3]; loci.append((c,e))
D1=load("dep_Denisova11.txt"); D2=load("dep_Goyet.txt"); D3=load("dep_LesCottes.txt")
with open(os.path.join(OUT,"lowcov_panel_depth.tsv"),"w") as o:
    o.write("chrom\tpos\tDenisova11\tGoyet\tLesCottes\n")
    for c,p in loci:
        o.write(f"{c}\t{p}\t{D1.get((c,p),'NA')}\t{D2.get((c,p),'NA')}\t{D3.get((c,p),'NA')}\n")
print("wrote", os.path.join(OUT,"lowcov_panel_depth.tsv"), "rows", len(loci))
PY
echo "done"; date
