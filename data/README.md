# Sequencing reads land here

Two scripts can fill this folder, both producing the same paired-end,
gzip-compressed FASTQ layout described below:

- `scripts/get_demo_data.sh` — a **tiny simulated** dataset (fast, no big
  downloads); the default for a live class.
- `scripts/get_reads.sh` — the **real** GEO **GSE135568** samples (see the
  table in the top-level README); slower and network-dependent.

Both write files named:

    <sample>_R1.fastq.gz
    <sample>_R2.fastq.gz

e.g. `ctrl_1_R1.fastq.gz`, `mip6_2_R2.fastq.gz`. Each is subsampled to 1,000,000
read pairs so the pipeline runs in minutes on a free machine.

To use your own data instead, drop paired FASTQ here with the same naming and
update `samples.tsv` + the `group` vector in `scripts/de_analysis.R`.
