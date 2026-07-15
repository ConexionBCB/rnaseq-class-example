#!/usr/bin/env python3
"""Simulate a tiny paired-end RNA-seq dataset for the teaching pipeline.

Reads are drawn from the yeast reference downloaded by scripts/00_get_data.sh
(ref/genome.fa + ref/genes.gtf). A built-in expression difference between the
two groups (control vs mip6) is baked in, so the downstream edgeR step finds
real differentially-expressed genes and the volcano / heatmap are meaningful.

This is deliberately simple and fast (no sequencing errors, fixed quality) so
the whole pipeline finishes in ~1 minute on a free Codespace or Binder.
For the REAL GSE135568 data instead, use scripts/get_reads.sh.
"""
import argparse
import gzip
import os
import random
import re

READLEN = 90        # length of each mate
FRAGLEN = 180       # insert size (R1 and R2 sit end-to-end, no overlap)
N_GENES = 250       # genes to simulate reads from
N_UP = 30           # genes up-regulated in mip6 vs control
N_DOWN = 30         # genes down-regulated in mip6 vs control
DE_FOLD = 6         # fold-change applied to the DE genes
QUAL = "I" * READLEN  # Phred 40 for every base

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def read_fasta(path):
    """Return {seqname: sequence} from a (plain) FASTA file."""
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_exons(path):
    """Longest exon per gene_id: {gene_id: (chrom, start, end)} (1-based)."""
    gid_re = re.compile(r'gene_id "([^"]+)"')
    best = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = gid_re.search(f[8])
            if not m:
                continue
            gid = m.group(1)
            chrom, start, end = f[0], int(f[3]), int(f[4])
            length = end - start + 1
            if gid not in best or length > (best[gid][2] - best[gid][1] + 1):
                best[gid] = (chrom, start, end)
    return best


def read_samples(path):
    """Return [(sample, group), ...] from samples.tsv (skips the header)."""
    rows = []
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0]:
                rows.append((parts[0], parts[2]))
    return rows


def pick_genes(exons, genome):
    """Deterministically choose N_GENES whose longest exon can hold a fragment."""
    ok = []
    for gid in sorted(exons):
        chrom, start, end = exons[gid]
        if chrom not in genome:
            continue
        if (end - start + 1) < FRAGLEN + 10:
            continue
        ok.append(gid)
        if len(ok) >= N_GENES:
            break
    return ok


def gene_counts(idx, group):
    """Expected read-pair count for gene #idx in the given group."""
    base = 50 + (idx % 11) * 20            # 50..250, varies across genes
    if idx < N_UP:                         # up in mip6
        return base * DE_FOLD if group == "mip6" else base
    if idx < N_UP + N_DOWN:                # down in mip6
        return max(5, base // DE_FOLD) if group == "mip6" else base
    return base                            # not differentially expressed


def sample_fragment(seq, start, end, rng):
    """A FRAGLEN window inside the 1-based [start, end] exon, avoiding Ns."""
    lo, hi = start - 1, end - FRAGLEN      # 0-based inclusive start range
    for _ in range(6):
        s = rng.randint(lo, hi)
        frag = seq[s:s + FRAGLEN]
        if "N" not in frag.upper():
            return frag
    return frag  # give up avoiding Ns; STAR will soft-clip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="samples.tsv")
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--genome", default="ref/genome.fa")
    ap.add_argument("--gtf", default="ref/genes.gtf")
    args = ap.parse_args()

    genome = read_fasta(args.genome)
    exons = read_exons(args.gtf)
    genes = pick_genes(exons, genome)
    if len(genes) < 50:
        raise SystemExit(
            f"Only found {len(genes)} usable genes in {args.gtf}; expected the "
            "yeast reference from scripts/00_get_data.sh."
        )

    samples = read_samples(args.samples)
    os.makedirs(args.outdir, exist_ok=True)
    print(f">>> simulating {len(genes)} genes across {len(samples)} samples")

    for si, (sample, group) in enumerate(samples):
        rng = random.Random(1000 + si)     # reproducible, independent per sample
        r1_path = os.path.join(args.outdir, f"{sample}_R1.fastq.gz")
        r2_path = os.path.join(args.outdir, f"{sample}_R2.fastq.gz")
        total = 0
        with gzip.open(r1_path, "wt") as o1, gzip.open(r2_path, "wt") as o2:
            for gi, gid in enumerate(genes):
                chrom, start, end = exons[gid]
                seq = genome[chrom]
                # small per-replicate noise so dispersion is estimable
                expected = gene_counts(gi, group) * rng.uniform(0.85, 1.15)
                n = max(0, int(round(expected)))
                for k in range(n):
                    frag = sample_fragment(seq, start, end, rng)
                    r1 = frag[:READLEN]
                    r2 = revcomp(frag[-READLEN:])
                    name = f"{sample}:{gid}:{total + k}"
                    o1.write(f"@{name}/1\n{r1}\n+\n{QUAL}\n")
                    o2.write(f"@{name}/2\n{r2}\n+\n{QUAL}\n")
                total += n
        print(f"    {sample:8s} ({group:7s}) -> {total:6d} read pairs")

    print(f">>> demo reads written to {args.outdir}/")


if __name__ == "__main__":
    main()
