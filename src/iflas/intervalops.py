"""Pure interval / exon-intron helpers — faithfully ported from legacy commonFuncs.py.

These underpin AS detection (find_as) and junction refinement (refine). Pure functions,
no I/O — unit-tested here so their behavior is locked before the stages are migrated.
"""

from __future__ import annotations


def sortTupleList(tupleList):
    return sorted(tupleList)



def isOverlap(listA, listB):
    sortedList = sorted([listA, listB])
    return True if sortedList[0][1] >= sortedList[1][0] else False


def getOverlapOfTuple(tupleListA, tupleListB):
    res = []
    i = j = 0
    while i < len(tupleListA) and j < len(tupleListB):
        lo = max(tupleListA[i][0], tupleListB[j][0])
        hi = min(tupleListA[i][1], tupleListB[j][1])
        if lo <= hi:
            res.append([lo, hi])
        if tupleListA[i][1] < tupleListB[j][1]:
            i += 1
        else:
            j += 1
    return res


def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]


def getBlockLength(blockList):
    return sum([int(x[1]) - int(x[0]) for x in blockList])


def getConsensusIntronN(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    """Count introns of structure 1 that match structure 2 within +/- offset."""
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    j, consensusN = 0, 0
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronStarts1[i] + offset:
                break
            if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronStarts1[i] + offset and \
                    intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                consensusN += 1
                j += 1
    return consensusN


def getSizes(starts, ends):
    return [int(ends[x]) - int(starts[x]) for x in range(len(starts))]


def getRelStarts(blockStarts):
    return [int(blockStarts[x]) - int(blockStarts[0]) for x in range(len(blockStarts))]


def listStr2Int(myList):
    return [int(i) for i in myList]


