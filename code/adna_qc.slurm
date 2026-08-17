#!/bin/bash
# =====================================================================================
# SEPIA — aDNA sample QC (Q7 + Q8) for the 3 low-coverage archaics (BAM samples).
#
#   Q7  heterozygosity (ANGSD -doSaf + realSFS, transversions-only) — the contamination
#       diagnostic: Neanderthals are expected LOW het; elevated het (esp. at low depth) flags
#       modern-human contamination. Plus sex (X/Y depth ratio) and mean depth.
#   Q8  mapDamage — post-mortem deamination authenticity: C→T at 5' ends / G→A at 3' ends.
#
# Runs on the FULL re-acquired run libraries per individual (ANGSD reads a bam-list, so no
# physical merge needed). High-cov archaics are VCF-only (no BAM) → their het comes from the
# VCFs separately; this job is the BAM-based QC the re-acquisition unblocks.
#
# Submit AFTER the BAM download, e.g.:
#   sbatch --dependency=afterok:<reacq_jobid> code/adna_qc.slurm
# =====================================================================================
#SBATCH --job-name=sepia_adna_qc
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/sepia/scratch/sepia_adna_qc_%j.log
set -uo pipefail

GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
REF="$GENOMES/reference/hs37d5.fa"
BAMROOT="$GENOMES/archaic/bam"
REPO="${SEPIA_REPO:-/nfs/turbo/lsa-tlasisi1/sepia}"
OUT="$REPO/output/adna_qc"; SCR="$REPO/scratch/adna_qc"
mkdir -p "$OUT" "$SCR"

# --- activate the QC env (angsd + mapdamage + samtools) ---
source ~/tlasisi/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate adna_qc 2>/dev/null || { echo "!! conda env 'adna_qc' not found — create it first"; exit 2; }
for t in angsd realSFS mapDamage samtools; do command -v "$t" >/dev/null || { echo "!! $t missing in adna_qc env"; exit 2; }; done
[ -f "$REF" ] || { echo "!! reference $REF missing"; exit 2; }
[ -f "$REF.fai" ] || samtools faidx "$REF"
echo "===== SEPIA aDNA QC (Q7/Q8) ====="; date; echo "angsd: $(command -v angsd)  mapDamage: $(command -v mapDamage)"

for c in $(seq 1 22); do echo "${c}:"; done > "$SCR/autosomes.rf"   # ANGSD region format

SUM="$OUT/adna_qc_summary.tsv"
printf "individual\tn_runs\tmean_depth_auto\tdepth_X\tdepth_Y\tsex_call\thet\thet_nsites\tdamage_5pCtoT_pos1\n" > "$SUM"

for G in Denisova11 GoyetQ56-1 LesCottes_Z4-1514; do
  dir="$BAMROOT/$G"; echo "===== $G ====="
  mapfile -t bams < <(for b in "$dir"/*.bam; do samtools quickcheck "$b" 2>/dev/null && echo "$b"; done)
  [ "${#bams[@]}" -gt 0 ] || { echo "  no valid BAMs"; printf "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" "$G" >> "$SUM"; continue; }
  printf '%s\n' "${bams[@]}" > "$SCR/$G.bamlist"; echo "  ${#bams[@]} run libraries"

  # ---- coverage + sex (stream all runs via samtools cat) ----
  samtools cat "${bams[@]}" 2>/dev/null | samtools coverage - 2>/dev/null > "$SCR/$G.coverage" || true
  read mdep xdep ydep sex < <(awk 'NR>1{
      if($1~/^[0-9]+$/){a+=$7*$3; al+=$3} if($1=="X"){x=$7} if($1=="Y"){y=$7}}
      END{am=(al>0)?a/al:0; xr=(am>0)?x/am:0; yr=(am>0)?y/am:0;
          s=(xr<0.75 && yr>0.05)?"XY_male":((xr>=0.75 && yr<0.05)?"XX_female":"uncertain");
          printf "%.3f %.3f %.3f %s", am, x+0, y+0, s}' "$SCR/$G.coverage")
  echo "  depth auto=$mdep X=$xdep Y=$ydep  -> $sex"

  # ---- heterozygosity (ANGSD SAF -> realSFS, transversions only) ----
  het="NA"; nsites="NA"
  angsd -bam "$SCR/$G.bamlist" -ref "$REF" -anc "$REF" -out "$SCR/$G.het" \
    -GL 1 -doSaf 1 -minMapQ 30 -minQ 30 -noTrans 1 -C 50 -baq 1 -rf "$SCR/autosomes.rf" -nThreads 8 \
    > "$SCR/$G.angsd.log" 2>&1 || echo "  angsd failed (see $G.angsd.log)"
  if [ -f "$SCR/$G.het.saf.idx" ]; then
    realSFS "$SCR/$G.het.saf.idx" -fold 1 > "$SCR/$G.sfs" 2> "$SCR/$G.realsfs.log" || true
    read het nsites < <(awk '{t=0; for(i=1;i<=NF;i++)t+=$i; if(t>0 && NF>=2) printf "%.6f %d", $2/t, t; else printf "NA 0"}' "$SCR/$G.sfs")
  fi
  echo "  het=$het (nsites=$nsites)"

  # ---- mapDamage authenticity (subsample the largest run for speed) ----
  dmg="NA"; big=$(ls -S "${bams[@]}" 2>/dev/null | head -1)
  if [ -n "$big" ]; then
    samtools view -s 0.02 -b "$big" -o "$SCR/$G.sub.bam" 2>/dev/null && samtools index "$SCR/$G.sub.bam" 2>/dev/null
    mapDamage -i "$SCR/$G.sub.bam" -r "$REF" -d "$OUT/$G.mapdamage" --no-stats -y 0.3 > "$SCR/$G.mapdamage.log" 2>&1 || echo "  mapDamage failed (see $G.mapdamage.log)"
    [ -f "$OUT/$G.mapdamage/5pCtoT_freq.txt" ] && dmg=$(awk 'NR==2{print $2}' "$OUT/$G.mapdamage/5pCtoT_freq.txt")
    rm -f "$SCR/$G.sub.bam"*
  fi
  echo "  5' C->T (pos1) = $dmg"

  printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$G" "${#bams[@]}" "$mdep" "$xdep" "$ydep" "$sex" "$het" "$nsites" "$dmg" >> "$SUM"
done

echo "===== QC SUMMARY ====="; column -t "$SUM" 2>/dev/null || cat "$SUM"
echo "Interpret: Neanderthals expect LOW het (elevated => contamination); 5' C->T should be"
echo "clearly elevated (authentic aDNA deamination). mapDamage plots in output/adna_qc/<ind>.mapdamage/."
date
