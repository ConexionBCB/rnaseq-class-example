# Sequencing reads land here

`scripts/get_reads.sh` fills this folder from GEO **GSE135568** (see the table
in the top-level README), producing paired-end, gzip-compressed FASTQ named:

    <sample>_R1.fastq.gz
    <sample>_R2.fastq.gz

e.g. `ctrl_1_R1.fastq.gz`, `mip6_2_R2.fastq.gz`. Each is subsampled to 1,000,000
read pairs so the pipeline runs in minutes on a free machine.

To use your own data instead, drop paired FASTQ here with the same naming and
update `samples.tsv` + the `group` vector in `scripts/de_analysis.R`.
