"""Unit tests for iflas.intervalops (pure exon/intron helpers; behavior lock for find_as/refine)."""

from iflas.intervalops import (
    getBlockLength, getConsensusIntronN, getIntrons, getOverlapOfTuple, getRelStarts,
    getSizes, isOverlap, listStr2Int, sortTupleList,
)


def test_sort_and_merge():
    assert sortTupleList([(3, 4), (1, 2)]) == [(1, 2), (3, 4)]


def test_overlap():
    assert isOverlap([1, 5], [3, 8]) is True
    assert isOverlap([1, 2], [5, 8]) is False
    assert getOverlapOfTuple([(1, 5), (10, 15)], [(3, 12)]) == [[3, 5], [10, 12]]
    assert getOverlapOfTuple([(1, 2)], [(5, 6)]) == []


def test_introns_and_sizes():
    assert getIntrons([100, 260], [150, 300]) == ([150], [260])
    assert getSizes([100, 260], [150, 300]) == [50, 40]
    assert getBlockLength([(100, 150), (260, 300)]) == 90
    assert getRelStarts([10, 170]) == [0, 160]


def test_consensus_introns():
    s, e = [100, 260, 400], [150, 300, 450]
    assert getConsensusIntronN(s, e, s, e, 0) == 2            # identical -> all introns match
    # shift the 2nd structure DOWN by 1bp (its intron starts stay <= ref, so the inner
    # scan doesn't break early): offset 0 -> no match, offset 2 -> both match
    s2 = [99, 259, 399]
    e2 = [149, 299, 449]
    assert getConsensusIntronN(s, e, s2, e2, 0) == 0
    assert getConsensusIntronN(s, e, s2, e2, 2) == 2


def test_list_conversions():
    assert listStr2Int(["1", "2", "3"]) == [1, 2, 3]
