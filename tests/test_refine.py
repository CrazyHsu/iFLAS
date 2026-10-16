"""Unit tests for refine's pure splice-site / junction-validity helpers.

These are faithful ports of the published refine.py logic and need no external
tools, so they lock the canonical-motif and exon-validity rules precisely.
(strandAdjust / readsAssign / identifyNovelIsoformsByJunctions are exercised
end-to-end in tests/test_refine_driver.py against a real reference.)
"""

import pytest

from iflas.stages.refine import (
    identifyNovelIsoformsByJunctions,
    isCanonicalSite,
    validJunction,
)


@pytest.mark.parametrize("dinuc,expected", [
    ("GT", 1), ("GC", 1), ("AT", 1),   # + donor: GT/GC/AT canonical
    ("AA", 0), ("CT", 0),
])
def test_canonical_plus_donor(dinuc, expected):
    assert isCanonicalSite("+", "donor", dinuc) == expected


@pytest.mark.parametrize("dinuc,expected", [
    ("AG", 1), ("AC", 1),              # + acceptor: AG/AC canonical
    ("GG", 0), ("TT", 0),
])
def test_canonical_plus_acceptor(dinuc, expected):
    # any non-"donor" type takes the acceptor branch (legacy passes "accept")
    assert isCanonicalSite("+", "acceptor", dinuc) == expected
    assert isCanonicalSite("+", "accept", dinuc) == expected


@pytest.mark.parametrize("dinuc,expected", [
    ("AC", 1), ("GC", 1), ("AT", 1),   # - donor: AC/GC/AT canonical
    ("GT", 0), ("AA", 0),
])
def test_canonical_minus_donor(dinuc, expected):
    assert isCanonicalSite("-", "donor", dinuc) == expected


@pytest.mark.parametrize("dinuc,expected", [
    ("CT", 1), ("GT", 1),              # - acceptor: CT/GT canonical
    ("AG", 0), ("AA", 0),
])
def test_canonical_minus_acceptor(dinuc, expected):
    assert isCanonicalSite("-", "acceptor", dinuc) == expected


def test_valid_junction_all_positive_exons():
    # start=100,end=900, intron 300-700 -> exons [200, 200], all > 0
    assert validJunction(100, 900, "300-700") is True
    # two introns -> three exons, all positive
    assert validJunction(2000, 3200, "2200-2500;2650-3000") is True


def test_valid_junction_zero_or_negative_exon():
    # start coincides with first donor -> first exon length 0 -> invalid
    assert validJunction(300, 900, "300-700") is False
    # start past first donor -> negative exon -> invalid
    assert validJunction(400, 900, "300-700") is False
    # middle exon collapses (intron2 start <= intron1 end) -> invalid
    assert validJunction(2000, 3200, "2200-2700;2650-3000") is False


# --- BUG-115: identifyNovelIsoformsByJunctions must not crash when the
# reference has no mono-exon transcript but a read is mono-exon. -------------

# genePredExt (bincolumn=False): one MULTI-exon transcript only -> the legacy
# loop never creates the "monoExon" key. Fields:
# transName chrom strand txStart txEnd cdsStart cdsEnd exonCount
#   exonStarts exonEnds score geneName cdsStartStat cdsEndStat exonFrames
_MULTI_EXON_GPE = "tx1\tchr1\t+\t1000\t2000\t1000\t2000\t2\t1000,1500,\t1200,2000,\t0\tgeneX\tcmpl\tcmpl\t0,0,\n"

# bed12+: a mono-exon read (blockCount=1) + a multi-exon read whose intron
# chain (chr1:1200-1500) matches tx1 above.
_MONO_READ = "chr1\t3000\t3500\tread_mono\t0\t+\t3000\t3500\t0,0,0\t1\t500,\t0,\tgeneX\n"
_MATCH_READ = "chr1\t1000\t2000\tread_anno\t0\t+\t1000\t2000\t0,0,0\t2\t200,500,\t0,500,\tgeneX\n"


def test_identify_novel_no_mono_exon_reference(tmp_path):
    """Reference with no mono-exon transcript + a mono-exon read.

    Legacy code does ``knownJuncDict["monoExon"]`` directly, which KeyErrors
    here. After the BUG-115 fix it must run and classify the mono-exon read as
    novel (no known mono-exon transcript to match it), while the multi-exon
    read whose intron chain is known is annotated.
    """
    gpe = tmp_path / "ref.gpe"
    bed = tmp_path / "reads.bed12+"
    gpe.write_text(_MULTI_EXON_GPE)
    bed.write_text(_MONO_READ + _MATCH_READ)
    anno = tmp_path / "anno.bed"
    novel = tmp_path / "novel.bed"

    identifyNovelIsoformsByJunctions(str(gpe), str(bed), anno=str(anno), novel=str(novel))

    anno_names = {l.split("\t")[3] for l in anno.read_text().splitlines() if l.strip()}
    novel_names = {l.split("\t")[3] for l in novel.read_text().splitlines() if l.strip()}
    assert novel_names == {"read_mono"}
    assert anno_names == {"read_anno"}
