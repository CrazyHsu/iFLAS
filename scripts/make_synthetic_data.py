"""Generate a tiny synthetic iFLAS dataset with KNOWN alternative-splicing events.

Ground truth (0-based, half-open exons):
  geneA  control, single isoform
  geneB  exon skipping (SE)
  geneC  alternative 5' splice site (A5SS, alt donor)
  geneE  alternative 3' splice site (A3SS, alt acceptor)
  geneD  intron retention (IR)
  genePU01..genePU20  PU-filter background pairs (canonical iso1 + novel SE variant)

Outputs (into <outdir>):
  genome.fa, genes.gtf, flnc.fa (full-length reads), ngs_R1/R2.fastq.gz, primer.fa

Usage:  python scripts/make_synthetic_data.py test_data/synthetic
"""

from __future__ import annotations

import gzip
import os
import random
import sys

CHROM_LEN = {"chr1": 6000, "chr2": 3000, "chr3": 22000}
PU_BACKGROUND_GENE_COUNT = 20


def _make_pu_background_genes(n: int = PU_BACKGROUND_GENE_COUNT) -> dict:
    """Extra isoform pairs used to keep PU-filter synthetic tests above OOB limits."""
    genes = {}
    for i in range(n):
        gene = f"genePU{i + 1:02d}"
        base = 500 + i * 1000
        exons = [(base, base + 120), (base + 300, base + 420), (base + 650, base + 780)]
        genes[gene] = [
            (f"{gene}.iso1", "chr3", "+", exons),
            (f"{gene}.iso2_SE", "chr3", "+", [exons[0], exons[2]]),
        ]
    return genes

# gene -> [(isoform_name, chrom, strand, [(exon_start, exon_end), ...]), ...]
CORE_GENES = {
    "geneA": [("geneA.iso1", "chr1", "+", [(100, 300), (800, 1000), (1500, 1700)])],
    "geneB": [("geneB.iso1", "chr1", "+", [(2000, 2200), (2500, 2650), (3000, 3200)]),
              ("geneB.iso2_SE", "chr1", "+", [(2000, 2200), (3000, 3200)])],
    "geneC": [("geneC.iso1", "chr1", "+", [(3700, 4000), (4500, 4700)]),
              ("geneC.iso2_A5SS", "chr1", "+", [(3700, 4050), (4500, 4700)])],
    "geneE": [("geneE.iso1", "chr1", "+", [(4300, 4600), (5000, 5200)]),
              ("geneE.iso2_A3SS", "chr1", "+", [(4300, 4600), (4950, 5200)])],
    "geneD": [("geneD.iso1", "chr2", "+", [(100, 300), (700, 900)]),
              ("geneD.iso2_IR", "chr2", "+", [(100, 900)])],
}
GENES = {**CORE_GENES, **_make_pu_background_genes()}

_RC = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _tx_seq(genome: dict, chrom: str, exons) -> str:
    return "".join(genome[chrom][st:en] for (st, en) in exons)


def generate(outdir: str, tgs_copies: int = 5, ngs_pairs: int = 15,
             read_len: int = 100, seed: int = 42) -> str:
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    genome = {c: _rand_seq(L, rng) for c, L in CHROM_LEN.items()}

    # make introns canonical (GT..AG) so splice-aware alignment places junctions exactly
    for chrom, introns in expected_introns().items():
        lst = list(genome[chrom])
        for (s, e) in introns:
            lst[s], lst[s + 1] = "G", "T"     # donor
            lst[e - 2], lst[e - 1] = "A", "G"  # acceptor
        genome[chrom] = "".join(lst)

    # genome.fa (60-col wrapped)
    with open(os.path.join(outdir, "genome.fa"), "w") as f:
        for c, s in genome.items():
            f.write(f">{c}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")

    # genes.gtf (1-based inclusive)
    with open(os.path.join(outdir, "genes.gtf"), "w") as f:
        for gene, isos in GENES.items():
            for name, chrom, strand, exons in isos:
                for (st, en) in exons:
                    f.write(f'{chrom}\tsyn\texon\t{st + 1}\t{en}\t.\t{strand}\t.\t'
                            f'gene_id "{gene}"; transcript_id "{name}";\n')

    # flnc.fa : full-length transcript reads (tgs_copies per isoform)
    with open(os.path.join(outdir, "flnc.fa"), "w") as f:
        for gene, isos in GENES.items():
            for name, chrom, strand, exons in isos:
                seq = _tx_seq(genome, chrom, exons)
                for k in range(tgs_copies):
                    f.write(f">{name}/ccs/{k}\n{seq}\n")

    # ngs_R1/R2.fastq.gz : paired short reads from transcripts
    r1p = os.path.join(outdir, "ngs_R1.fastq.gz")
    r2p = os.path.join(outdir, "ngs_R2.fastq.gz")
    qual = "I" * read_len
    with gzip.open(r1p, "wt") as o1, gzip.open(r2p, "wt") as o2:
        for gene, isos in GENES.items():
            for name, chrom, strand, exons in isos:
                seq = _tx_seq(genome, chrom, exons)
                if len(seq) < 2 * read_len:
                    continue
                for p in range(ngs_pairs):
                    start = rng.randint(0, len(seq) - 2 * read_len)
                    frag_len = min(len(seq) - start, rng.randint(2 * read_len, 400))
                    frag = seq[start:start + frag_len]
                    r1 = frag[:read_len]
                    r2 = revcomp(frag[-read_len:])
                    rid = f"{name}_{p}"
                    o1.write(f"@{rid}/1\n{r1}\n+\n{qual}\n")
                    o2.write(f"@{rid}/2\n{r2}\n+\n{qual}\n")

    # primer.fa (dummy IsoSeq primers)
    with open(os.path.join(outdir, "primer.fa"), "w") as f:
        f.write(">5p\nAAGCAGTGGTATCAACGCAGAGTAC\n>3p\nGTACTCTGCGTTGATACCACTGCTT\n")

    return outdir


# expected introns (0-based start=end-of-exon_i, end=start-of-exon_{i+1}) per chrom, for tests
def expected_introns() -> dict:
    introns: dict[str, set] = {}
    for gene, isos in GENES.items():
        for name, chrom, strand, exons in isos:
            for i in range(len(exons) - 1):
                introns.setdefault(chrom, set()).add((exons[i][1], exons[i + 1][0]))
    return introns


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "test_data/synthetic"
    generate(out)
    print(f"synthetic dataset written to {out}")
