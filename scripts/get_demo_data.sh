#!/usr/bin/env bash
# Generate a TINY, self-contained demo dataset so the whole pipeline runs in
# ~1 minute with NO large downloads — ideal for a live Codespaces/Binder class.
#
# It simulates paired-end reads from the yeast reference (downloaded by
# 00_get_data.sh) with a built-in expression difference between the two groups,
# so STAR -> featureCounts -> edgeR produce a real volcano plot and heatmap.
#
# For the REAL GSE135568 data instead, use scripts/get_reads.sh.
set -euo pipefail

# --- make sure the pipeline's tools are on PATH --------------------------
# In a fresh Codespaces/Binder terminal the 'rnaseq' conda env may not be
# active yet. Activate it if needed so python is the env's python.
if [ "${CONDA_DEFAULT_ENV:-}" != "rnaseq" ]; then
  __base="$(conda info --base 2>/dev/null || echo /opt/conda)"
  # shellcheck disable=SC1091
  source "${__base}/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate rnaseq 2>/dev/null || true
fi

# The simulator reads the reference produced by 00_get_data.sh. Build it first
# if it isn't there yet (small, reliable Ensembl download — unlike the SRA runs).
if [ ! -f ref/genome.fa ] || [ ! -f ref/genes.gtf ]; then
  echo ">>> reference not found — running scripts/00_get_data.sh first…"
  bash scripts/00_get_data.sh
fi

mkdir -p data
python scripts/simulate_reads.py --samples samples.tsv --outdir data

echo ">>> demo reads ready in data/  (run scripts/run_pipeline.sh next)"
