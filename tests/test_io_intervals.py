"""Unit tests for iflas.io.intervals (ported BED parsers).

Fixture: a 2-exon transcript on chr1, + strand.
  exon1 100-150, intron 150-260, exon2 260-300; thick(CDS) 120-280.
Expected values are hand-computed to lock the ported behavior.
"""

from iflas.io.intervals import Bed6, Bed12, Bed12Plus, Junction

BED12 = "chr1\t100\t300\ttx1\t0\t+\t120\t280\t0,0,0\t2\t50,40\t0,160"


def test_bed6_parsing():
    b = Bed6("chr2\t10\t20\tfeat\t30\t-")
    assert (b.chrom, b.chromStart, b.chromEnd, b.name, b.score, b.strand) == \
        ("chr2", 10, 20, "feat", 30, "-")
    assert b.pos == "chr2:10-20"


def test_bed12_structure():
    b = Bed12(BED12)
    assert b.exonStarts == [100, 260]
    assert b.exonEnds == [150, 300]
    assert b.introns == [(150, 260)]
    assert b.exons == [(100, 150), (260, 300)]
    assert b.exonChain == "100-150;260-300"
    assert b.juncChain == "150-260"


def test_bed12_lengths():
    b = Bed12(BED12)
    assert len(b) == 90                 # 50 + 40
    assert b.get_all_exons_len() == 90
    assert b.cds_len() == 50            # CDS 120-280 intersected with exons
    assert b.utr_5_len() == 20          # + strand: 100..120
    assert b.utr_3_len() == 20          # + strand: 280..300


def test_bed12_repr_roundtrip():
    b = Bed12(BED12)
    b2 = Bed12(repr(b))
    assert b2.exonChain == b.exonChain
    assert b2.juncChain == b.juncChain
    assert len(b2) == len(b)
    assert (b2.chromStart, b2.chromEnd, b2.strand) == (100, 300, "+")


def test_bed12plus_other_columns():
    b = Bed12Plus(BED12 + "\textraA\textraB")
    assert b.otherList == ["extraA", "extraB"]
    assert b.exonChain == "100-150;260-300"


def test_junction_donor_acceptor():
    j = Junction(BED12)
    assert j.donor == 150
    assert j.acceptor == 260
    assert j.jPos == "chr1:150-260"
    assert (j.dn, j.ac) == (150, 260)   # + strand
    assert j.dnPos == "chr1:150:+"
