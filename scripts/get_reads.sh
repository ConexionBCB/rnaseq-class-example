#!/usr/bin/env bash
# Download the teaching dataset and subsample it for class use.
#
# Dataset: GEO GSE135568 (S. cerevisiae).
#   controls = heat-shocked wild-type
#   cases    = mip6-delta  (knockout of Mip6, an RNA-binding protein)
#
# Reads samples.tsv  (sample <TAB> gsm <TAB> group)  and, for each row,
# resolves GSM -> SRR, downloads it, and keeps N random read pairs so the
# whole pipeline finishes in minutes on a free machine.
#
# Requires: pysradb, sra-tools, seqtk  (all in environment.yml).
set -euo pipefail

N="${N:-1000000}"          # read pairs to keep per sample
mkdir -p data raw

tail -n +2 samples.tsv | while IFS=$'\t' read -r SAMPLE GSM GROUP; do
  [ -z "${SAMPLE:-}" ] && continue
  echo "==== ${SAMPLE}  (${GSM}, ${GROUP}) ===="

  SRR=$(pysradb gsm-to-srr "$GSM" | awk 'NR==2{print $2}')
  if [ -z "$SRR" ]; then echo "could not resolve $GSM" >&2; exit 1; fi
  echo ">>> ${GSM} -> ${SRR}"

  prefetch -O raw "$SRR"
  fasterq-dump --split-files -O raw "raw/${SRR}/${SRR}.sra"

  if [ ! -f "raw/${SRR}_2.fastq" ]; then
    echo "WARNING: ${SRR} is single-end; this pipeline is paired-end." >&2
    echo "         Use the single-end branch noted in README.md." >&2
    seqtk sample -s100 "raw/${SRR}.fastq" "$N" | gzip > "data/${SAMPLE}_R1.fastq.gz"
    continue
  fi

  # same seed for both mates keeps pairs in sync
  seqtk sample -s100 "raw/${SRR}_1.fastq" "$N" | gzip > "data/${SAMPLE}_R1.fastq.gz"
  seqtk sample -s100 "raw/${SRR}_2.fastq" "$N" | gzip > "data/${SAMPLE}_R2.fastq.gz"
done

echo ">>> reads ready in data/  (run scripts/run_pipeline.sh next)"
