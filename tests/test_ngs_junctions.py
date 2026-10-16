"""NGS junction generation on synthetic data (srMapping + getJuncFromRegtools port).

Chain: hisat2 (NGS short reads) -> per-repeat stringtie -> stringtie --merge +
samtools cat -> tmp.bam -> regtools junctions -> junctions.bed. The synthetic NGS
reads are sampled from spliced transcripts, so regtools should recover the designed
introns. tmp.bam + junctions.bed are what pu_filter / find_as NGS filtering consume.

Requires hisat2 + stringtie + regtools + samtools + bedtools (core iflas env).
"""

import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_synthetic_data as syn  # noqa: E402

from iflas.io.intervals import Junction
from iflas.stages.align import ngs_reassembly_junctions

_TOOLS = ["hisat2", "hisat2-build", "stringtie", "regtools", "samtools", "bedtools"]
pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in _TOOLS),
                                reason="needs hisat2/stringtie/regtools/samtools/bedtools")


def _recovered_introns(junctions_bed):
    introns = set()
    with open(junctions_bed) as f:
        for line in f:
            if not line.strip():
                continue
            j = Junction(line)
            introns.add((j.chrom, j.donor, j.acceptor))
    return introns


def test_ngs_reassembly_junctions_recovers_designed_introns(tmp_path):
    d = str(tmp_path)
    syn.generate(d)
    idx = os.path.join(d, "hisat2_idx")
    subprocess.run(["hisat2-build", "-q", f"{d}/genome.fa", idx], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    out = ngs_reassembly_junctions(
        idx, os.path.join(d, "reassembly"),
        repeats=[{"left": f"{d}/ngs_R1.fastq.gz", "right": f"{d}/ngs_R2.fastq.gz"}],
        threads=2, filter_count=0)

    assert os.path.isfile(out["tmp_bam"]) and os.path.getsize(out["tmp_bam"]) > 0
    assert os.path.isfile(out["tmp_bam"] + ".bai")
    assert os.path.isfile(out["junctions_bed"])

    recovered = _recovered_introns(out["junctions_bed"])
    expected = {(c, s, e) for c, introns in syn.expected_introns().items() for (s, e) in introns}
    # NGS reads are sampled from spliced transcripts -> at least some designed
    # introns must be recovered, and every recovered junction must be a real one.
    assert recovered & expected, f"no designed introns recovered; got {recovered}"
    assert recovered <= expected, f"unexpected (false) junctions: {recovered - expected}"
