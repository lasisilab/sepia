#!/bin/bash
# =====================================================================================
# SEPIA — panel-locus read depth for the 3 LOW-coverage archaics, by REMOTE extraction.
#
# The high-cov archaics carry DP in their snpAD VCFs (already in output/panel_matrix/), but
# Denisova 11, Goyet Q56-1 and Les Cottés Z4-1514 are BAM-only — and the local copies from
# acquire_archaic.slurm are TRUNCATED (incomplete downloads, missing BGZF EOF). Rather than
# re-download ~100 GB of full BAMs, we stream ONLY the reads overlapping the 222 panel
# positions straight from ENA's complete files (samtools 1.21 libcurl + .bai range requests),
# summing all of each individual's run libraries. Result reads from the authoritative files,
# so depth is correct even though our local BAMs are broken.
#
# Output: output/panel_matrix/lowcov_panel_depth.tsv  (chrom pos Denisova11 Goyet LesCottes)
# Run from the cluster checkout:  sbatch code/lowcov_panel_depth.slurm
# =====================================================================================
#SBATCH --job-name=sepia_lowcov_depth
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_lowcov_depth_%j.log
set -uo pipefail

REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
MAP="$REPO/data/panel_hg19_hg38_map.tsv"
OUT="$REPO/output/panel_matrix"; SCR="$REPO/scratch/lowcov_depth"
mkdir -p "$OUT" "$SCR/parts"; cd "$SCR"

module load Bioinformatics >/dev/null 2>&1 || true
module load samtools       >/dev/null 2>&1 || true
command -v samtools >/dev/null || { echo "!! samtools not found"; exit 2; }
[ -f "$MAP" ] || { echo "!! map not found: $MAP"; exit 2; }
echo "===== SEPIA low-cov archaic panel depth (remote) ====="; date

# ---- panel hg19 BED (0-based); ENA BAMs are hs37d5/bare contigs ----
awk -F'\t' 'NR>1 && $4!="" && $4!="NA"{print $2"\t"($4-1)"\t"$4}' "$MAP" | sort -k1,1 -k2,2n -u > panel.bed
echo "panel positions: $(wc -l < panel.bed)"

# sample -> ENA project (each project = one individual, verified via sample_title)
declare -A PROJ=( [Denisova11]=PRJEB24663 [LesCottes]=PRJEB21875 [Goyet]=PRJEB21870 )

bam_urls_for(){  # $1=accession -> https URLs of every run BAM (skip .bai)
  curl -sfL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=run_accession,submitted_ftp&format=tsv" 2>/dev/null \
    | awk -F'\t' 'NR>1{n=split($2,a,";"); for(i=1;i<=n;i++) if(a[i] ~ /\.bam$/) print "https://"a[i]}'
}

extract_run(){  # $1=out.bam  $2=url  — fetch only panel-region reads, retry 3x
  local out="$1" url="$2" t
  for t in 1 2 3; do
    if samtools view -b -M -L panel.bed "$url" -o "$out" 2>/dev/null && samtools quickcheck "$out" 2>/dev/null; then
      return 0
    fi
    echo "      retry $t: $(basename "$url")"; sleep $((5*t))
  done
  echo "      FAILED (skipped): $url"; rm -f "$out"; return 1
}

for G in Denisova11 LesCottes Goyet; do
  acc="${PROJ[$G]}"; echo "===== $G ($acc) ====="
  mapfile -t urls < <(bam_urls_for "$acc")
  echo "  runs at ENA: ${#urls[@]}"
  [ "${#urls[@]}" -gt 0 ] || { echo "  !! no BAM URLs resolved for $acc"; continue; }
  parts=(); i=0
  for url in "${urls[@]}"; do
    i=$((i+1)); ob="parts/${G}.$(printf '%02d' $i).bam"
    echo "    [$i/${#urls[@]}] $(basename "$url")"
    extract_run "$ob" "$url" && parts+=("$ob")
  done
  [ "${#parts[@]}" -gt 0 ] || { echo "  !! no panel reads extracted for $G"; continue; }
  # merge this individual's run parts, sort, depth over all panel positions (-a keeps zeros)
  samtools merge -f "$G.merged.bam" "${parts[@]}" 2>/dev/null
  samtools sort -@ 4 -o "$G.sorted.bam" "$G.merged.bam" 2>/dev/null
  samtools index "$G.sorted.bam" 2>/dev/null
  samtools depth -a -b panel.bed "$G.sorted.bam" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"}{sub(/^chr/,"",$1); print $1,$2,$3}' | sort -k1,1 -k2,2n > "dep_$G.txt"
  echo "  $G: $(wc -l < dep_$G.txt) positions, mean depth $(awk '{n++;t+=$3}END{printf "%.2f",(n?t/n:0)}' dep_$G.txt)"
done

# ---- merge the three onto the panel loci ----
echo "===== merge -> lowcov_panel_depth.tsv ====="
python3 - "$OUT" <<'PY'
import sys, os
OUT=sys.argv[1]
def load(f):
    d={}
    if os.path.exists(f):
        for ln in open(f):
            p=ln.split()
            if len(p)>=3: d[(p[0],p[1])]=p[2]
    return d
loci=[]
for ln in open("panel.bed"):
    c,s,e=ln.split()[:3]; loci.append((c,e))
D1=load("dep_Denisova11.txt"); D2=load("dep_Goyet.txt"); D3=load("dep_LesCottes.txt")
with open(os.path.join(OUT,"lowcov_panel_depth.tsv"),"w") as o:
    o.write("chrom\tpos\tDenisova11\tGoyet\tLesCottes\n")
    for c,p in loci:
        o.write(f"{c}\t{p}\t{D1.get((c,p),'NA')}\t{D2.get((c,p),'NA')}\t{D3.get((c,p),'NA')}\n")
print("wrote", os.path.join(OUT,"lowcov_panel_depth.tsv"), "rows", len(loci))
PY
echo "done"; date
