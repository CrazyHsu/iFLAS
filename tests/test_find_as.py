"""Unit tests for the find_as AS-event detectors on synthetic ground truth.

These drive the four detectors (findIR/findSE/findA5SS/findA3SS) directly on
hand-built Gene2Reads objects whose isoform structures match the designed
synthetic dataset (scripts/make_synthetic_data.py). No external tools needed —
the detector logic is pure and deterministic, so this locks the ported
algorithm's behavior precisely.

Ground truth (0-based, half-open exons; all + strand):
  geneB SE   : iso1 [(2000,2200),(2500,2650),(3000,3200)] vs iso2 [(2000,2200),(3000,3200)]
  geneC A5SS : iso1 [(3700,4000),(4500,4700)]   vs iso2 [(3700,4050),(4500,4700)]   (alt donor 4000/4050)
  geneE A3SS : iso1 [(4300,4600),(5000,5200)]   vs iso2 [(4300,4600),(4950,5200)]   (alt acceptor 5000/4950)
  geneD IR   : iso1 [(100,300),(700,900)]       vs iso2 [(100,900)]                 (intron 300-700 retained)
"""

from io import StringIO

import pytest

from iflas.io.intervals import Bed12
from iflas.io.records import Gene2Reads
from iflas.stages import find_as


def _bed12(name, chrom, strand, exons):
    chrom_start = exons[0][0]
    chrom_end = exons[-1][1]
    sizes = ",".join(str(e - s) for s, e in exons)
    rel_starts = ",".join(str(s - chrom_start) for s, e in exons)
    fields = [chrom, chrom_start, chrom_end, name, 0, strand,
              chrom_start, chrom_end, "0,0,0", len(exons), sizes, rel_starts]
    return Bed12("\t".join(map(str, fields)))


def _gene2reads(gname, chrom, strand, isoforms):
    g = Gene2Reads(gname)
    for name, exons in isoforms:
        g.update(_bed12(name, chrom, strand, exons))
    return g


def _run(detector, as_type, gene2reads):
    out = StringIO()
    detector({gene2reads.geneName: gene2reads}, as_type, out, anno=False, offset=5)
    lines = [l.split("\t") for l in out.getvalue().splitlines()]
    return lines


def test_find_se_detects_skipped_exon():
    g = _gene2reads("geneB", "chr1", "+",
                    [("geneB.iso1", [(2000, 2200), (2500, 2650), (3000, 3200)]),
                     ("geneB.iso2_SE", [(2000, 2200), (3000, 3200)])])
    rows = _run(find_as.findSE, "SE", g)
    assert len(rows) == 1
    row = rows[0]
    assert row[0] == "chr1" and row[1] == "2500" and row[2] == "2650"
    assert row[3] == "geneB:2200@2500-2650@3000"
    assert row[4] == "500.0"
    assert row[5] == "+"
    assert row[12] == "geneB"
    assert row[13] == "2200" and row[14] == "3000"   # flanking splice sites
    assert row[15] == "geneB.iso1" and row[16] == "1"   # keep (inclusion) reads
    assert row[17] == "geneB.iso2_SE" and row[18] == "1"  # skip (exclusion) reads


def test_find_ir_detects_retained_intron():
    g = _gene2reads("geneD", "chr2", "+",
                    [("geneD.iso1", [(100, 300), (700, 900)]),
                     ("geneD.iso2_IR", [(100, 900)])])
    rows = _run(find_as.findIR, "IR", g)
    assert len(rows) == 1
    row = rows[0]
    assert row[0] == "chr2" and row[1] == "300" and row[2] == "700"
    assert row[3] == "geneD:300-700"
    assert row[4] == "500.0"
    assert row[5] == "+"
    assert row[6] == "geneD"
    assert row[7] == "geneD.iso2_IR" and row[8] == "1"   # retention reads
    assert row[9] == "geneD.iso1" and row[10] == "1"     # spliced reads


def test_find_a5ss_detects_alt_donor():
    g = _gene2reads("geneC", "chr1", "+",
                    [("geneC.iso1", [(3700, 4000), (4500, 4700)]),
                     ("geneC.iso2_A5SS", [(3700, 4050), (4500, 4700)])])
    rows = _run(find_as.findA5SS, "A5SS", g)
    assert len(rows) == 1
    row = rows[0]
    assert row[0] == "chr1" and row[1] == "4000" and row[2] == "4050"
    assert row[3] == "geneC:4000-4050"
    assert row[4] == "500.0"
    assert row[5] == "+"
    assert row[6] == "geneC"
    assert row[7] == "geneC.iso2_A5SS" and row[8] == "1"   # inclusion (longer donor)
    assert row[9] == "geneC.iso1" and row[10] == "1"       # exclusion
    assert row[11] == "inc:4050-4500"
    assert row[12] == "exc:4000-4500"


def test_find_a3ss_detects_alt_acceptor():
    g = _gene2reads("geneE", "chr1", "+",
                    [("geneE.iso1", [(4300, 4600), (5000, 5200)]),
                     ("geneE.iso2_A3SS", [(4300, 4600), (4950, 5200)])])
    rows = _run(find_as.findA3SS, "A3SS", g)
    assert len(rows) == 1
    row = rows[0]
    assert row[0] == "chr1" and row[1] == "4950" and row[2] == "5000"
    assert row[3] == "geneE:4950-5000"
    assert row[4] == "500.0"
    assert row[5] == "+"
    assert row[6] == "geneE"
    assert row[7] == "geneE.iso2_A3SS" and row[8] == "1"   # inclusion (earlier acceptor)
    assert row[9] == "geneE.iso1" and row[10] == "1"       # exclusion
    assert row[11] == "inc:4600-4950"
    assert row[12] == "exc:4600-5000"


def test_no_event_for_single_isoform_gene():
    g = _gene2reads("geneA", "chr1", "+",
                    [("geneA.iso1", [(100, 300), (800, 1000), (1500, 1700)])])
    assert _run(find_as.findSE, "SE", g) == []
    assert _run(find_as.findIR, "IR", g) == []
    assert _run(find_as.findA5SS, "A5SS", g) == []
    assert _run(find_as.findA3SS, "A3SS", g) == []


@pytest.mark.parametrize("offset,expected", [(0, 12)])
def test_consensus_intron_n_identical_structures(offset, expected):
    starts = [100, 800, 1500]
    ends = [300, 1000, 1700]
    # identical structures: 2 introns matched under IR+SE, 2 sites x2 under A3SS+A5SS
    n = find_as.getConsensusIntronNfromAS(starts, ends, starts, ends, offset)
    assert n == expected


def test_consensus_intron_n_disjoint_structures():
    starts1, ends1 = [100, 800, 1500], [300, 1000, 1700]
    starts2, ends2 = [5000], [5100]   # single exon, no introns
    assert find_as.getConsensusIntronNfromAS(starts1, ends1, starts2, ends2, 0) == 0
