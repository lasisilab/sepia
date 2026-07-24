#!/bin/bash
# =====================================================================================
# SEPIA — acquire archaic-genome data into the SHARED lab resource
#   $SEPIA_GENOMES/archaic/{vcf,mask,bam}/<Genome>/   (default /nfs/turbo/lsa-tlasisi1/genomes)
#
#   HIGH-coverage genomes -> published damage-aware snpAD genotype VCFs + FilterBed masks
#                            (no BAM needed; genotypes + heterozygosity come from the VCF).
#   LOW-coverage genomes   -> full BAMs (for ANGSD genotype-likelihoods + mapDamage).
#
# Resume-safe (curl -C -), integrity-checked (rebuild index / re-fetch truncated), md5 +
# provenance recorded to MANIFEST.tsv (the archaic half of Q1). Everything is hg19/hs37d5;
# the project reports on hg38 later via the rsID map (item B2 — NOT this job).
# Included subset per plan/papers/SAMPLE_INCLUSION.md.
#
# Run from the cluster checkout:
#   cd /nfs/turbo/lsa-tlasisi1/sepia && sbatch code/acquire_archaic.slurm
# Re-running is safe: it resumes partial files and skips complete ones.
# =====================================================================================
#SBATCH --job-name=sepia_acquire_archaic
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/genomes/archaic/acquire_%j.log
set -uo pipefail

GENOMES="${SEPIA_GENOMES:-/nfs/turbo/lsa-tlasisi1/genomes}"
SRC=/nfs/turbo/lsa-tlasisi1/lheald_thesis/aDNA_data          # for the hs37d5 reference copy
A="$GENOMES/archaic"; MAN="$A/MANIFEST.tsv"
mkdir -p "$A"/vcf "$A"/mask "$A"/bam "$GENOMES/reference"
[ -f "$MAN" ] || printf "path\tsource_url\tmd5\tbuild\tdate\n" > "$MAN"

module load Bioinformatics >/dev/null 2>&1 || true
module load bcftools/1.21  >/dev/null 2>&1 || module load bcftools >/dev/null 2>&1 || true
module load samtools       >/dev/null 2>&1 || true
echo "===== SEPIA archaic acquisition ====="; date; echo "dest: $A"
for t in curl bcftools samtools md5sum; do printf "%s: " "$t"; command -v "$t" || echo MISSING; done

fetch(){  # $1=url  $2=dest  — resume-safe download + manifest row (skip if already complete)
  [ -s "$2" ] && return 0
  curl -C - -fSL --retry 3 --retry-delay 5 -s "$1" -o "$2" || { echo "  MISS $1"; return 1; }
  printf "%s\t%s\t%s\thg19\t%s\n" "$2" "$1" "$(md5sum "$2" 2>/dev/null | cut -d' ' -f1)" "$(date +%F)" >> "$MAN"
}
vcf_ok(){  # ensure an index exists, but ONLY build it if missing. The old `-f` rebuilt the
           # index on EVERY re-run => re-read all ~257 GB of VCFs off NFS each time (hours),
           # never reaching the BAMs. An index already present == the file was validated on a
           # prior run, so skip it. New/truncated downloads have no index yet and get built once.
  { [ -f "$1.tbi" ] || [ -f "$1.csi" ]; } && return 0
  bcftools index "$1" >/dev/null 2>&1
}

# ---- 0. reference (hs37d5, hg19) ----------------------------------------------------
if [ ! -s "$GENOMES/reference/hs37d5.fa" ]; then
  echo "===== 0. copy reference hs37d5.fa ====="
  cp -n "$SRC/reference/hs37d5.fa"     "$GENOMES/reference/" 2>/dev/null || echo "  (reference copy failed — fetch hs37d5 manually)"
  cp -n "$SRC/reference/hs37d5.fa.fai" "$GENOMES/reference/" 2>/dev/null
fi

CHRS="$(seq 1 22) X"

# ---- 1. HIGH-COVERAGE snpAD genotype VCFs (+ .tbi) ----------------------------------
echo "===== 1. high-coverage snpAD VCFs ====="
VB=http://cdna.eva.mpg.de/neandertal/Vindija/VCF          # Altai, Vindija33.19, Denisova(=Den3), Mez1
for G in Altai Vindija33.19 Denisova Mez1; do
  mkdir -p "$A/vcf/$G"
  for c in $CHRS; do
    f="chr${c}_mq25_mapab100.vcf.gz"
    fetch "$VB/$G/$f" "$A/vcf/$G/$f" && fetch "$VB/$G/$f.tbi" "$A/vcf/$G/$f.tbi"
    vcf_ok "$A/vcf/$G/$f" || { echo "  RE-FETCH truncated $G/$f"; rm -f "$A/vcf/$G/$f" "$A/vcf/$G/$f".{tbi,csi}; fetch "$VB/$G/$f" "$A/vcf/$G/$f" && vcf_ok "$A/vcf/$G/$f" || echo "  BAD $G/$f"; }
  done
  echo "  $G: $(ls "$A/vcf/$G"/*.vcf.gz 2>/dev/null | wc -l) chr VCFs"
done
# Chagyrskaya 8 — separate release; use http (ftp egress is blocked; http 302-redirects, -L follows)
mkdir -p "$A/vcf/Chagyrskaya8"
CH=http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF
for c in $CHRS; do f="chr${c}.noRB.vcf.gz"; fetch "$CH/$f" "$A/vcf/Chagyrskaya8/$f" && fetch "$CH/$f.tbi" "$A/vcf/Chagyrskaya8/$f.tbi"; vcf_ok "$A/vcf/Chagyrskaya8/$f"; done
echo "  Chagyrskaya8: $(ls "$A/vcf/Chagyrskaya8"/*.vcf.gz 2>/dev/null | wc -l) chr VCFs (0 => ftp egress blocked; use xfer node)"
# Denisova 25 — PROVISIONAL (preprint; may be gated / not yet public)
mkdir -p "$A/vcf/Denisova25"
D25=http://cdna.eva.mpg.de/denisova/Den25/VCF
for c in $CHRS Y; do f="chr${c}.Den25.L35MQ25.B30.map35_100.vcf.gz"; fetch "$D25/$f" "$A/vcf/Denisova25/$f" && vcf_ok "$A/vcf/Denisova25/$f"; done
echo "  Denisova25: $(ls "$A/vcf/Denisova25"/*.vcf.gz 2>/dev/null | wc -l) chr VCFs (0 => gated/not-public — defer, per plan)"

# ---- 2. FilterBed masks (best-effort: probe each release's dir listing) --------------
echo "===== 2. FilterBed masks ====="
probe_masks(){ mkdir -p "$2"; for u in $(curl -sfL "$1/" 2>/dev/null | grep -oiE 'href="[^"?/][^"]*\.bed(\.gz)?"' | sed -E 's/href="//;s/"//'); do fetch "$1/$u" "$2/$(basename "$u")"; done; }
# Vindija-tree FilterBed is per-genome subdirs (FilterBed/Altai/, /Denisova/, /Vindija33.19/)
for G in Altai Denisova Vindija33.19 Mez1; do probe_masks "http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/$G" "$A/mask/$G"; done
probe_masks "http://ftp.eva.mpg.de/neandertal/Chagyrskaya/FilterBed"  "$A/mask/Chagyrskaya8"
probe_masks "http://cdna.eva.mpg.de/denisova/Den25/FilterBed"         "$A/mask/Denisova25"
echo "  masks downloaded: $(find "$A/mask" -type f -name '*.bed*' 2>/dev/null | wc -l) (0 => confirm mask filenames from the release README)"

# ---- 3. LOW-COVERAGE BAMs (resolve via ENA filereport; verify + index) ---------------
echo "===== 3. low-coverage BAMs ====="
resolve_bam(){  # $1=genome  $2=ENA-accession
  local G="$1" acc="$2"; mkdir -p "$A/bam/$G"
  echo "  $G ($acc): resolving BAM URL via ENA filereport..."
  local urls; urls=$(curl -sfL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=submitted_ftp&format=tsv" 2>/dev/null | tr ';' '\n' | grep -iE '\.bam$' | sed -E 's#^(ftp://)?#https://#')
  [ -z "$urls" ] && { echo "    no submitted BAM at ENA (likely FASTQ, or EVA-only) -> resolve from EVA mirror manually"; return; }
  for u in $urls; do b=$(basename "$u"); fetch "$u" "$A/bam/$G/$b" && samtools quickcheck "$A/bam/$G/$b" 2>/dev/null && samtools index "$A/bam/$G/$b" 2>/dev/null && echo "    OK $b" || echo "    BAD/again $b"; done
}
resolve_bam LesCottes_Z4-1514 PRJEB21875
resolve_bam Denisova11        PRJEB24663
resolve_bam GoyetQ56-1        PRJEB21870

echo "===== acquisition summary ====="; date
echo "  VCF genomes: $(ls -d "$A"/vcf/*/ 2>/dev/null | wc -l)   BAM genomes with data: $(find "$A/bam" -name '*.bam' 2>/dev/null | wc -l)   masks: $(find "$A/mask" -type f 2>/dev/null | wc -l)"
echo "  manifest rows: $(($(wc -l < "$MAN") - 1))"; du -sh "$A" 2>/dev/null
echo "NEXT: apply FilterBed + intersect to a common mask; extract pigmentation loci on hg38 (B2a);"
echo "      ANGSD GLs + mapDamage for the low-cov BAMs. Any '0'/'MISS' lines above => confirm URLs and re-run (resume-safe)."
