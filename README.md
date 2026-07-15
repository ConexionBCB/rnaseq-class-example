# rnaseq-class — a teaching RNA-seq pipeline in the cloud

A minimal, **runnable** RNA-seq pipeline for teaching, designed to fit inside
free cloud tiers (GitHub Codespaces and MyBinder), with a matching Google Colab
notebook.

```
FastQC → Trimmomatic → STAR → featureCounts → edgeR → volcano + heatmap
```

The reference is *Saccharomyces cerevisiae* (yeast, ~12 Mb) on purpose: the
genome and STAR index fit in ~2 GB of RAM, so the whole thing runs on a free
machine.

## The teaching dataset

**GEO [GSE135568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135568)** — yeast.

| sample | GSM | group |
| --- | --- | --- |
| ctrl_1 | GSM4017984 | control — heat-shocked wild-type |
| ctrl_2 | GSM4017991 | control |
| ctrl_3 | GSM4018003 | control |
| mip6_1 | GSM4017990 | case — *mip6Δ* (Mip6 RNA-binding protein KO) |
| mip6_2 | GSM4017999 | case |
| mip6_3 | GSM4018006 | case |

`samples.tsv` holds this mapping. `scripts/get_reads.sh` resolves each GSM to
its SRR run (via `pysradb`), downloads it, and keeps 1,000,000 random read
pairs per sample so the pipeline finishes in minutes. Change `N=` to keep more.

The contrast in `scripts/de_analysis.R` is **mip6Δ vs control**; edit the
`group` vector there if you change the samples.

> If the GEO runs turn out to be single-end, `get_reads.sh` warns you. The
> pipeline as written is paired-end — for single-end, drop the `R2`/PE
> arguments from the Trimmomatic and STAR calls in `run_pipeline.sh`.

## What's here

| File | Purpose |
| --- | --- |
| `.devcontainer/devcontainer.json` | Codespaces / Dev Containers definition |
| `environment.yml` | conda (bioconda) environment — all the tools |
| `postBuild` | MyBinder build hook (builds the index into the image) |
| `samples.tsv` | GSM → sample → group mapping |
| `scripts/00_get_data.sh` | download yeast reference + build STAR index |
| `scripts/get_reads.sh` | download & subsample the GSE135568 samples |
| `scripts/run_pipeline.sh` | FastQC → Trimmomatic → STAR → featureCounts → R |
| `scripts/de_analysis.R` | edgeR + volcano plot + heatmap |
| `rnaseq-class.ipynb` | the same pipeline as a Google Colab notebook |
| `data/` | subsampled FASTQ land here |

## Run it in GitHub Codespaces

1. Push this repo to GitHub.
2. **Code ▸ Codespaces ▸ Create codespace on main.** The container builds and
   the `rnaseq` conda environment is created automatically (first build ≈ 3–5 min).
3. **Open a _new_ terminal** once the build finishes — the prompt should read
   `(rnaseq)`, meaning the environment is active. Then run (use `bash`, not `sh`):
   ```bash
   bash scripts/00_get_data.sh      # reference + STAR index
   bash scripts/get_reads.sh        # download the GSE135568 samples
   bash scripts/run_pipeline.sh     # the full pipeline + figures
   ```
   Results: `counts/de_results.csv`, `figures/volcano.png`, `figures/heatmap.png`.

   The scripts also re-activate the `rnaseq` env themselves, so they still work
   in a terminal that was opened before the build finished.

### Troubleshooting

- **`STAR: not found` (or any tool not found)** — the `rnaseq` environment
  isn't active. Open a fresh terminal (prompt shows `(rnaseq)`), or run
  `conda activate rnaseq`. Always launch scripts with `bash scripts/…`, never
  `sh scripts/…` — `sh` ignores the env and the script's `bash` shebang.
- **"Cannot find R… change `r.rpath.linux`"** — this is fixed by the R-path
  setting in `.devcontainer/devcontainer.json`. If you see it once on the very
  first build (before the env exists), reload the window (Command Palette ▸
  *Developer: Reload Window*) after the build completes.

## Run it on MyBinder

Make the repo public, then open
`https://mybinder.org/v2/gh/<you>/rnaseq-class/HEAD`.
`postBuild` builds the index at image-build time. In a terminal run
`bash scripts/get_reads.sh && bash scripts/run_pipeline.sh`. Download results
before you leave — Binder does not persist. (Downloading six SRA runs may be
tight in 2 GB; use fewer samples or a smaller `N` on Binder.)

## Run it in Google Colab

Open `rnaseq-class.ipynb` in Colab, set the repo URL in the first code cell,
and **Runtime ▸ Run all**. Cell 1 installs conda (the runtime restarts once).
