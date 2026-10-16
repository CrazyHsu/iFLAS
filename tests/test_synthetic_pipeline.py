"""Integration smoke test on synthetic data: generate -> init-ref -> splice-align -> junctions.

Validates the synthetic dataset is biologically coherent (designed AS events recovered)
and that the runnable pipeline parts (init-ref, alignment) work end-to-end with env tools.
"""

import os
import subprocess
import sys

import pysam

from iflas.cli import main

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402


def _introns(bam_path):
    out = set()
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam:
        if r.is_unmapped:
            continue
        blocks = r.get_blocks()
        for i in range(len(blocks) - 1):
            out.add((r.reference_name, blocks[i][1], blocks[i + 1][0]))
    return out


def test_synthetic_generate_initref_align(tmp_path):
    d = str(tmp_path / "syn")
    syn.generate(d)
    for fn in ["genome.fa", "genes.gtf", "flnc.fa", "ngs_R1.fastq.gz", "ngs_R2.fastq.gz", "primer.fa"]:
        assert os.path.getsize(os.path.join(d, fn)) > 0

    # init-ref derives aux files
    assert main(["init-ref", "-f", f"{d}/genome.fa", "--gtf", f"{d}/genes.gtf", "-o", d]) == 0
    assert os.path.isfile(f"{d}/genome.mm2") and os.path.isfile(f"{d}/genome.len.txt")

    # splice-align the full-length reads
    bam = f"{d}/flnc.bam"
    subprocess.run(f"minimap2 -ax splice -uf {d}/genome.fa {d}/flnc.fa | samtools sort -o {bam} -",
                   shell=True, check=True, capture_output=True)
    pysam.index(bam)
    introns = _introns(bam)

    # every designed (spliced) intron is recovered exactly
    for chrom, iset in syn.expected_introns().items():
        for s, e in iset:
            assert (chrom, s, e) in introns, f"missing designed intron {chrom}:{s}-{e}"

    # AS-diagnostic junctions are present
    assert ("chr1", 2200, 3000) in introns   # SE (exon-skipping) junction
    assert ("chr1", 4050, 4500) in introns    # A5SS alternative donor
    assert ("chr1", 4600, 4950) in introns    # A3SS alternative acceptor
