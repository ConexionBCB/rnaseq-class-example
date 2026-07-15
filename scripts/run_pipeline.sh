#!/usr/bin/env bash
# End-to-end RNA-seq pipeline:
#   FastQC -> Trimmomatic -> STAR -> featureCounts -> (edgeR + plots in R)
#
# Expects paired-end reads in data/ named:  <sample>_R1.fastq.gz / <sample>_R2.fastq.gz
# First run:  scripts/00_get_data.sh  (reference + index)  and
#             scripts/get_reads.sh    (download the GSE135568 samples)
set -euo pipefail

# --- make sure the pipeline's tools are on PATH --------------------------
# In a fresh Codespaces/Binder terminal the 'rnaseq' conda env may not be
# active yet — that is what causes "STAR: not found". Activate it if needed.
if [ "${CONDA_DEFAULT_ENV:-}" != "rnaseq" ]; then
  __base="$(conda info --base 2>/dev/null || echo /opt/conda)"
  # shellcheck disable=SC1091
  source "${__base}/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate rnaseq 2>/dev/null || true
fi

THREADS="${THREADS:-4}"
IDX=ref/star_index
GTF=ref/genes.gtf

# Trimmomatic's bundled adapters live under share/. bioconda ships them in a
# version-stamped directory (e.g. share/trimmomatic-0.39-2/), so resolve the
# path with a glob and fall back gracefully instead of hard-coding a version.
BIN_DIR="$(dirname "$(command -v trimmomatic)")"
ADAPTERS="${BIN_DIR}/../share/trimmomatic/adapters/TruSeq3-PE.fa"
if [ ! -f "$ADAPTERS" ]; then
  ADAPTERS="$(ls "${BIN_DIR}/../share/"trimmomatic*/adapters/TruSeq3-PE.fa 2>/dev/null | head -n1 || true)"
fi
if [ -z "${ADAPTERS:-}" ] || [ ! -f "$ADAPTERS" ]; then
  echo "Could not locate Trimmomatic's TruSeq3-PE.fa adapters. Is the 'rnaseq' env active?" >&2
  exit 1
fi

mkdir -p qc trimmed aligned counts

shopt -s nullglob
reads=(data/*_R1.fastq.gz)
if [ ${#reads[@]} -eq 0 ]; then
  echo "No reads found. Run scripts/get_reads.sh, or put <sample>_R1/_R2.fastq.gz in data/." >&2
  exit 1
fi

for R1 in "${reads[@]}"; do
  SAMPLE="$(basename "$R1" _R1.fastq.gz)"
  R2="data/${SAMPLE}_R2.fastq.gz"
  echo "==== ${SAMPLE} ===="

  echo ">>> FastQC"
  fastqc -q -t "$THREADS" -o qc "$R1" "$R2"

  echo ">>> Trimmomatic"
  trimmomatic PE -threads "$THREADS" "$R1" "$R2" \
    "trimmed/${SAMPLE}_R1.fq.gz" "trimmed/${SAMPLE}_R1.unpaired.fq.gz" \
    "trimmed/${SAMPLE}_R2.fq.gz" "trimmed/${SAMPLE}_R2.unpaired.fq.gz" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

  echo ">>> STAR"
  STAR --runThreadN "$THREADS" --genomeDir "$IDX" \
    --readFilesIn "trimmed/${SAMPLE}_R1.fq.gz" "trimmed/${SAMPLE}_R2.fq.gz" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "aligned/${SAMPLE}_"
done

echo ">>> featureCounts (all samples into one matrix)"
featureCounts -T "$THREADS" -p --countReadPairs -a "$GTF" \
  -o counts/counts.txt aligned/*_Aligned.sortedByCoord.out.bam

echo ">>> edgeR + plots (volcano, heatmap)"
Rscript scripts/de_analysis.R

echo ">>> pipeline complete. See counts/de_results.csv and figures/*.png"
