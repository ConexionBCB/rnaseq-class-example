#!/usr/bin/env bash
# Runs ONCE, when the Codespace / Dev Container is first created.
#   1) build the bioconda 'rnaseq' environment (all pipeline tools + R)
#   2) make that environment activate automatically in every new terminal,
#      so students never hit "STAR: not found" / "cannot find R".
set -euo pipefail

echo ">>> creating the 'rnaseq' conda environment (first build ~3-5 min)…"
conda env create -f environment.yml

# Make `conda activate` work in interactive shells (idempotent), then
# auto-activate rnaseq. Guarded so re-runs don't append duplicate lines.
conda init bash
if ! grep -q 'conda activate rnaseq' ~/.bashrc; then
  echo 'conda activate rnaseq' >> ~/.bashrc
fi

echo ">>> done. Open a NEW terminal — the prompt should read (rnaseq)."
