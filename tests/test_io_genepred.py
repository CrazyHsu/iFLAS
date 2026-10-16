"""Unit tests for iflas.io.genepred (ported GPE parsers) + intervals<->genepred round-trip.

GPE fixture mirrors the BED12 fixture: chr1 +, tx 100-300, CDS 120-280, exons 100-150 & 260-300.
"""

import textwrap

from iflas.io.genepred import GenePredExtLine, GenePredObj
from iflas.io.intervals import Bed12

# gpeExt without bin column (15 fields)
GPE = "tx1\tchr1\t+\t100\t300\t120\t280\t2\t100,260,\t150,300,\t0\tgene1\tcmpl\tcmpl\t0,0,"
BED12 = "chr1\t100\t300\ttx1\t0\t+\t120\t280\t0,0,0\t2\t50,40\t0,160"


def test_gpe_parsing():
    g = GenePredExtLine(GPE, bincolumn=False)
    assert g.transName == "tx1"
    assert g.chrom == "chr1" and g.strand == "+"
    assert g.txStart == 100 and g.txEnd == 300
    assert g.cdsStart == 120 and g.cdsEnd == 280
    assert g.exonStarts == [100, 260] and g.exonEnds == [150, 300]
    assert g.blockSizes == [50, 40] and g.blockStarts == [0, 160]
    assert g.introns == [(150, 260)]
    assert g.geneName == "gene1"


def test_gpe_lengths_and_flags():
    g = GenePredExtLine(GPE, bincolumn=False)
    assert len(g) == 90
    assert g.cds_len() == 50
    assert g.utr_5_len() == 20
    assert g.utr_3_len() == 20
    assert g.is_standard() is True


def test_gpe_bincolumn():
    g = GenePredExtLine("585\t" + GPE, bincolumn=True)
    assert g.bin == "585"
    assert g.transName == "tx1" and g.exonStarts == [100, 260]


def test_gpe_repr_roundtrip():
    g = GenePredExtLine(GPE, bincolumn=False)
    g2 = GenePredExtLine(repr(g), bincolumn=False)
    assert (g2.txStart, g2.txEnd, g2.exonStarts, g2.exonEnds) == \
        (g.txStart, g.txEnd, g.exonStarts, g.exonEnds)


def test_gpe_to_bed_structure():
    g = GenePredExtLine(GPE, bincolumn=False)
    bed = g.go_to_bed()
    assert isinstance(bed, Bed12)
    assert bed.chromStart == 100 and bed.chromEnd == 300
    assert bed.thickStart == 120 and bed.thickEnd == 280
    assert bed.blockSizes == [50, 40] and bed.blockStarts == [0, 160]
    assert repr(bed).split("\t")[:3] == ["chr1", "100", "300"]



def test_genepredobj_from_file(tmp_path):
    gpe_file = tmp_path / "sample.gpe"
    gpe_file.write_text(textwrap.dedent(GPE) + "\n")
    obj = GenePredObj(str(gpe_file), bincolumn=False)
    assert "gene1" in obj.geneName2gpeObj
    assert len(obj.geneName2gpeObj["gene1"]) == 1
    assert "chr1" in obj.genePredDict and "+" in obj.genePredDict["chr1"]
    beds = obj.toBed()
    assert len(beds) == 1 and beds[0].startswith("chr1\t100\t300")
