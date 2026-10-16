"""Unit tests for iflas.io.records (BedFile, Gene2Reads, ReadLineStruc).

Fixture: a BED12+ read line with geneName in column 13 (chr1 +, tx 100-300,
CDS 120-280, exons 100-150 & 260-300, gene1).
"""

from iflas.io.records import BedFile, Gene2Reads, ReadLineStruc

BED12P = "chr1\t100\t300\tread1\t0\t+\t120\t280\t0,0,0\t2\t50,40\t0,160\tgene1"


def test_readlinestruc_parsing():
    r = ReadLineStruc(BED12P)
    assert r.name == "read1"
    assert r.geneName == "gene1"
    assert r.exonLen == 90
    assert r.readLen == 200
    assert r.introns == [(150, 260)]



def test_gene2reads_aggregation():
    g = Gene2Reads("gene1")
    g.update(ReadLineStruc(BED12P))
    assert g.chrom == "chr1" and g.strand == "+"
    assert g.minpos == 100 and g.maxpos == 300
    assert g.readNames == ["read1"] and len(g.reads) == 1
    assert str(g) == "chr1:100-300 (gene1) (+) (1 reads)"


def test_gene2reads_asdict_types_do_not_share_mutable_defaults():
    g = Gene2Reads("gene1")
    g.asDict["IR"]["event1"] = "read1"

    assert g.asDict["IR"] == {"event1": "read1"}
    assert g.asDict["SE"] == {}
    assert g.asDict["A3SS"] == {}
    assert g.asDict["A5SS"] == {}


def test_bedfile_views(tmp_path):
    bed = tmp_path / "reads.bed12+"
    bed.write_text(BED12P + "\n")
    bf = BedFile(str(bed), type="bed12+")
    assert set(bf.reads) == {"read1"}
    assert bf.reads["read1"].otherList == ["gene1"]

    assert set(bf.getAllJuncDict(onebase=True)) == {"chr1:151-260"}
    assert set(bf.getAllExonDict()) == {"chr1:100-150", "chr1:260-300"}
    assert bf.getJunc2trans(onebase=True) == {"chr1:151-260": ["read1"]}
    assert bf.getGenePos(bedType="bed12+", geneCol=13) == {"gene1": ["chr1", 120, 280]}
