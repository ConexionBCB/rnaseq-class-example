#!/usr/bin/env bash
# Download a SMALL reference (S. cerevisiae, ~12 Mb) and build a STAR index.
# Yeast is deliberate: the genome + index fit inside free-tier / Binder RAM,
# so the whole pipeline runs in a classroom without a paid machine.
# Swap the URLs for your organism of choice once students have the idea.
set -euo pipefail

# --- make sure the pipeline's tools are on PATH --------------------------
# In Codespaces the tools live in the 'rnaseq' conda env; on MyBinder
# repo2docker installs them into the base env (already on PATH). Activate
# 'rnaseq' only if a tool is missing — a clean no-op on Binder.
if ! command -v STAR >/dev/null 2>&1; then
  __base="$(conda info --base 2>/dev/null || echo /opt/conda)"
  # shellcheck disable=SC1091
  source "${__base}/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate rnaseq 2>/dev/null || true
fi

mkdir -p ref
cd ref

REL=110
DNA="https://ftp.ensembl.org/pub/release-${REL}/fasta/saccharomyces_cerevisiae/dna"
GTF="https://ftp.ensembl.org/pub/release-${REL}/gtf/saccharomyces_cerevisiae"

echo ">>> genome FASTA"
curl -sL "${DNA}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" -o genome.fa.gz
gunzip -f genome.fa.gz

echo ">>> gene annotation GTF"
curl -sL "${GTF}/Saccharomyces_cerevisiae.R64-1-1.${REL}.gtf.gz" -o genes.gtf.gz
gunzip -f genes.gtf.gz

echo ">>> STAR index"
# genomeSAindexNbases must be reduced for small genomes:
#   min(14, log2(genomeLength)/2 - 1)  ->  ~10 for the 12 Mb yeast genome.
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile genes.gtf \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 10 \
  --runThreadN "${THREADS:-4}"

echo ">>> done. Reference in ref/, index in ref/star_index/"
