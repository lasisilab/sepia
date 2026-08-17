#!/bin/bash
# =====================================================================================
# SEPIA — re-acquire the LOW-coverage archaic BAMs in FULL (for Q7/Q8 aDNA QC).
#
# The copies from acquire_archaic.slurm are truncated (missing BGZF EOF) AND incomplete
# (several run libraries were never downloaded): its fetch() skipped any non-empty file, so
# it never resumed partial files, and resolve_bam only printed "BAD" without re-fetching.
# This re-pulls EVERY run library per individual from ENA with proper resume + integrity
# checks, so ANGSD (heterozygosity/contamination/sex) + mapDamage can run on complete BAMs.
#
# (The panel *depth* was already recovered by remote streaming — code/lowcov_panel_depth.slurm.
#  This job is only needed for the whole-genome QC steps, which need whole BAMs.)
#
# Resume-safe: re-running continues partial files and skips complete+valid ones.
#   cd /nfs/turbo/lsa-tlasisi1/sepia && sbatch code/reacquire_lowcov_bams.slurm
# =====================================================================================
#SBATCH --job-name=sepia_reacq_bam
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/genomes/archaic/reacquire_%j.log
set -uo pipefail

GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
A="$GENOMES/archaic"; BAMROOT="$A/bam"; MAN="$A/MANIFEST.tsv"
mkdir -p "$BAMROOT"
[ -f "$MAN" ] || printf "path\tsource_url\tmd5\tbuild\tdate\n" > "$MAN"
module load Bioinformatics >/dev/null 2>&1 || true
module load samtools       >/dev/null 2>&1 || true
command -v samtools >/dev/null || { echo "!! samtools not found"; exit 2; }
command -v curl     >/dev/null || { echo "!! curl not found"; exit 2; }
echo "===== SEPIA low-cov BAM re-acquisition ====="; date; df -h "$BAMROOT" | tail -1

# dir name (as already on disk) -> ENA project (each project = one individual, verified)
declare -A PROJ=( [Denisova11]=PRJEB24663 [GoyetQ56-1]=PRJEB21870 [LesCottes_Z4-1514]=PRJEB21875 )

get_bam(){  # $1=url  $2=dest  $3=expected_bytes  — resume+verify, up to 4 tries (fresh after 2 fails)
  local url="$1" dest="$2" want="$3" t
  if [ -f "$dest" ] && samtools quickcheck "$dest" 2>/dev/null; then
    echo "    OK (have) $(basename "$dest") [$(du -h "$dest"|cut -f1)]"; return 0
  fi
  for t in 1 2 3 4; do
    curl -C - -fSL --retry 5 --retry-delay 15 --connect-timeout 60 "$url" -o "$dest" 2>/dev/null || true
    if samtools quickcheck "$dest" 2>/dev/null; then
      local got; got=$(stat -c%s "$dest" 2>/dev/null || echo 0)
      [ -n "$want" ] && [ "$want" != "" ] && [ "$got" != "$want" ] && echo "    size ${got} != ENA ${want} (ok if ENA field blank)"
      echo "    fetched $(basename "$dest") [$(du -h "$dest"|cut -f1)] (try $t)"; return 0
    fi
    echo "    quickcheck failed (try $t)"; [ "$t" -ge 2 ] && { echo "    dropping partial → fresh"; rm -f "$dest"; }
    sleep $((15*t))
  done
  echo "    STILL BAD $(basename "$dest") — skipped"; return 1
}

for G in "${!PROJ[@]}"; do
  acc="${PROJ[$G]}"; dir="$BAMROOT/$G"; mkdir -p "$dir"
  echo "===== $G ($acc) ====="
  # resolve every run's BAM url + byte size from ENA
  curl -sfL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=run_accession,submitted_ftp,submitted_bytes&format=tsv" 2>/dev/null \
    | awk -F'\t' 'NR>1{ n=split($2,f,";"); split($3,b,";"); for(i=1;i<=n;i++) if(f[i] ~ /\.bam$/) print "https://"f[i]"\t"b[i] }' \
    > "$dir/.runs.tsv"
  echo "  runs at ENA: $(wc -l < "$dir/.runs.tsv")"
  while IFS=$'\t' read -r url bytes; do
    base=$(basename "$url"); dest="$dir/$base"
    if get_bam "$url" "$dest" "$bytes"; then
      # index: prefer ENA's .bai (small) else build
      if [ ! -f "$dest.bai" ]; then
        curl -fSL --retry 3 --connect-timeout 60 "$url.bai" -o "$dest.bai" 2>/dev/null || samtools index "$dest" 2>/dev/null || true
      fi
      grep -qF "$dest" "$MAN" 2>/dev/null || printf "%s\t%s\t%s\thg19\t%s\n" "$dest" "$url" "$(md5sum "$dest" 2>/dev/null|cut -d' ' -f1)" "$(date +%F)" >> "$MAN"
    fi
  done < "$dir/.runs.tsv"
  rm -f "$dir/.runs.tsv"
  echo "  $G now: $(ls "$dir"/*.bam 2>/dev/null | wc -l) BAMs, $(du -sh "$dir" 2>/dev/null | cut -f1)"
done

echo "===== re-acquisition summary ====="; date
echo "  total BAMs: $(find "$BAMROOT" -name '*.bam' 2>/dev/null | wc -l)   size: $(du -sh "$BAMROOT" 2>/dev/null | cut -f1)"
echo "  any still failing quickcheck:"; for b in "$BAMROOT"/*/*.bam; do samtools quickcheck "$b" 2>/dev/null || echo "    BAD $b"; done
echo "NEXT: ANGSD GLs + realSFS heterozygosity, contamination, X/Y sex (Q7); mapDamage (Q8)."
