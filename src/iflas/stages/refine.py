"""refine stage — strand validation + read/isoform-to-gene assignment (ported from refine.py).

The behavioral functions (isCanonicalSite/validJunction/identifyNovelIsoformsByJunctions/
strandAdjust/readsAssign/filterJunctionFile/refineWithJunc) and the splice-motif helper
keep the published algorithm apart from the fixes listed in CHANGELOG.md; their results are pinned by
tests/test_refine_characterization.py.

The driver ``refine_isoforms`` reproduces the LR core of legacy ``refineJunc`` with
explicit paths (analog to find_as_lr). The legacy perl shims are replaced with the
Python equivalents the project already standardized on (gpe2bed.pl -> GenePredObj.toBed).

Outputs the REAL ``refine/isoformGrouped.bed12+`` (bed12 + gene) and
``tofu.collapsed.assigned.unambi.bed12+`` — the latter feeds pu_filter (iso_pu.getInputBed).
NOTE: ``hqIsoforms/hq.collapsed.bed12+`` (read by find_as) is produced LATER by pu_filter,
not here. The NGS-junction refinement (filterJunctionFile/refineWithJunc) is ported but
runs only when NGS junctions are supplied (the synthetic/TGS-only path skips it).
"""

from __future__ import annotations

import copy
import os
import re

import pybedtools

from ..intervalops import (getBlockLength, getConsensusIntronN, getIntrons,
                           getOverlapOfTuple, isOverlap)
from ..io.files import makeLink, removeFiles, resolveDir
from ..io.genepred import GenePredObj
from ..io.intervals import Bed12, Bed12Plus
from ..io.records import BedFile
from ..pu.cpc2 import _find_utils_dir
from ..runtime import ToolRunner, getCurrentTime, quote as q


def getSpliceMotifFromJuncDict(juncDict, genomeFasta, withStrand=False):
    dinucleotideBedList = []
    for juncName in juncDict:
        juncBed = juncDict[juncName]
        chrom = juncBed.split("\t")[0]
        juncStart = int(juncBed.split("\t")[1])
        juncEnd = int(juncBed.split("\t")[2])
        strand = juncBed.split("\t")[-1]
        leftDinucleotidePos = "\t".join(map(str, [chrom, juncStart, juncStart + 2, ":".join([juncName, "left"]), ".", strand]))
        rightDinucleotidePos = "\t".join(map(str, [chrom, juncEnd - 2, juncEnd, ":".join([juncName, "right"]), ".", strand]))
        dinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])
    dinucleotideBedObj = pybedtools.BedTool("\n".join(dinucleotideBedList), from_string=True)
    dinucleotideBedRes = dinucleotideBedObj.sequence(genomeFasta, name=True, tab=True, s=True)

    juncMotifDict = {}
    for i in str(open(dinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
        infoList = str(i).strip("\n").split("\t")
        strand = infoList[0][infoList[0].find("(") + 1:infoList[0].find(")")]
        if not withStrand:
            juncName = "{}:{}".format(":".join(infoList[0].split(":")[:2]), strand)
            dinucleotideType = infoList[0].split(":")[2]
        else:
            juncName = ":".join(infoList[0].split(":")[:3])
            dinucleotideType = infoList[0].split(":")[3]

        if juncName not in juncMotifDict:
            juncMotifDict.update({juncName: {dinucleotideType: infoList[1]}})
        else:
            juncMotifDict[juncName].update({dinucleotideType: infoList[1]})

    junc2motif = {}
    for juncName in juncMotifDict:
        strand = juncName.split(":")[-1]
        if strand == "+":
            spliceMotif = "{}-{}".format(juncMotifDict[juncName]["left"], juncMotifDict[juncName]["right"])
        else:
            spliceMotif = "{}-{}".format(juncMotifDict[juncName]["right"], juncMotifDict[juncName]["left"])
        if juncName not in junc2motif:
            junc2motif[juncName] = spliceMotif

    return junc2motif


def isCanonicalSite(strand, dinucleotideType, dinucleotide):
    dinucleotide = dinucleotide.upper()                 # soft-masked genomes hold repeats in lowercase
    if strand == "+":
        if dinucleotideType == "donor":
            if re.match("(G[TC])|(AT)", dinucleotide):
                return 1
        else:
            if re.match("A[GC]", dinucleotide):
                return 1
    else:
        if dinucleotideType == "donor":
            if re.match("([AG]C)|(AT)", dinucleotide):
                return 1
        else:
            if re.match("[CG]T", dinucleotide):
                return 1
    return 0


def validJunction(start, end, juncChain):
    exonLen = []
    juncStarts = [int(i.split("-")[0]) for i in juncChain.split(";")]
    juncEnds = [int(i.split("-")[1]) for i in juncChain.split(";")]
    exonLen.append(juncStarts[0] - int(start))
    for i in range(0, len(juncStarts)-1):
        exonLen.append(juncStarts[i+1] - juncEnds[i])
    exonLen.append(int(end) - juncEnds[-1])
    return False if sum([i <= 0 for i in exonLen]) else True


def _intron_chain(chrom, strand, introns):
    return chrom + strand + ":" + ";".join(f"{st}-{ed}" for st, ed in introns)


def _known_structures(gpeFile):
    """Intron chains of the reference transcripts and, per chromosome, the spans of its single-exon transcripts."""
    chains, monoExon = set(), {}
    for transList in GenePredObj(gpeFile, bincolumn=False).geneName2gpeObj.values():
        for trans in transList:
            if trans.introns:
                chains.add(_intron_chain(trans.chrom, trans.strand, trans.introns))
            else:
                monoExon.setdefault(trans.chrom, []).append((trans.txStart, trans.txEnd))
    return chains, monoExon


def _is_known(read, chains, monoExon):
    """A multi-exon read with a reference intron chain, or a single-exon read covering at least half of a reference
    single-exon transcript of its chromosome (BUG-115: none there -> novel)."""
    if read.introns:
        return _intron_chain(read.chrom, read.strand, read.introns) in chains      # BUG-196: per strand
    span = (read.chromStart, read.chromEnd)
    for tx in monoExon.get(read.chrom, ()):
        # BUG-111 (fixed at CP3): the real overlap; legacy used a hard-coded 30-bp pair
        if isOverlap(tx, span) and float(getBlockLength(getOverlapOfTuple([span], [tx]))) / getBlockLength([tx]) >= 0.5:
            return True
    return False


def identifyNovelIsoformsByJunctions(gpeFile, bedFile, anno="annoIsoform.bed", novel="novelIsoform.bed"):
    chains, monoExon = _known_structures(gpeFile)
    reads = BedFile(bedFile, type="bed12+").reads
    with open(anno, "w") as annoOut, open(novel, "w") as novelOut:
        for read in reads.values():
            print(read, file=annoOut if _is_known(read, chains, monoExon) else novelOut)


def strandAdjust(genomeFasta, refGPE, bedFile, juncDiffScore, minCoverage, strandAdjust=None, strandConfirmed=None):
    tmpBed = "tmp.bed"
    GenePredObj(refGPE).toBed(outFile=tmpBed)

    strandAdjustOut = open(strandAdjust, "w")
    strandConfirmedOut = open(strandConfirmed, "w")

    bedObj = BedFile(bedFile, type="bed12+")
    gpeBedObj = pybedtools.bedtool.BedTool(tmpBed)
    singleExonReadBedList = []
    multiExonReadDinucleotideBedList = []
    for readName in bedObj.reads:
        read = bedObj.reads[readName]
        chrom, start, end, strand = read.chrom, read.chromStart, read.chromEnd, read.strand
        readExonStarts, readExonEnds = read.exonStarts, read.exonEnds
        readRegionBed = " ".join(map(str, [chrom, start, end, readName, ".", strand]))
        blockCount = read.blockCount
        if blockCount == 1:
            singleExonReadBedList.append(readRegionBed)
        else:
            readJuncStarts, readJuncEnds = getIntrons(readExonStarts, readExonEnds)
            for i in range(len(readJuncStarts)):
                leftDinucleotidePos = " ".join(map(str, [chrom, readJuncStarts[i], readJuncStarts[i] + 2, ":".join([readName, "junction"+str(i), "left"]), ".", strand]))
                rightDinucleotidePos = " ".join(map(str, [chrom, readJuncEnds[i] - 2, readJuncEnds[i], ":".join([readName, "junction"+str(i), "right"]), ".", strand]))
                multiExonReadDinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])

    # single exon process
    singleExonReadBedObj = pybedtools.bedtool.BedTool("\n".join(singleExonReadBedList), from_string=True)
    singleExonIntersectRes = gpeBedObj.intersect(singleExonReadBedObj, wa=True, wb=True)
    singleExonStrandDict = {}
    for i in singleExonIntersectRes:
        infoList = str(i).strip("\n").split("\t")
        transBed12 = Bed12("\t".join(infoList[:12]))
        transExonStarts, transExonEnds = transBed12.exonStarts, transBed12.exonEnds
        transExons = [(transExonStarts[x], transExonEnds[x]) for x in range(len(transExonStarts))]
        readExonStarts, readExonEnds = bedObj.reads[infoList[-3]].exonStarts, bedObj.reads[infoList[-3]].exonEnds
        readExons = [(readExonStarts[x], readExonEnds[x]) for x in range(len(readExonStarts))]
        overlapLen = getBlockLength(getOverlapOfTuple(readExons, transExons))
        coverageOnRead = overlapLen / float(getBlockLength((readExons)))
        coverageOnTrans = overlapLen / float(getBlockLength(transExons))
        if coverageOnRead >= minCoverage or coverageOnTrans >= minCoverage:
            if infoList[-3] not in singleExonStrandDict:
                singleExonStrandDict[infoList[-3]] = [transBed12.strand]
            else:
                singleExonStrandDict[infoList[-3]].append(transBed12.strand)
        else:
            singleExonStrandDict.setdefault(infoList[-3], [bedObj.reads[infoList[-3]].strand])

    # multi Exon process
    multiExonReadDinucleotideBedObj = pybedtools.BedTool("\n".join(multiExonReadDinucleotideBedList), from_string=True)
    multiExonReadDinucleotideBedRes = multiExonReadDinucleotideBedObj.sequence(genomeFasta, name=True, tab=True)
    multiExonReadStrandDict = {}
    for i in str(open(multiExonReadDinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
        infoList = str(i).strip("\n").split("\t")
        rName, dinucleotideType = infoList[0].split(":")[0], infoList[0].split(":")[2]
        if rName not in multiExonReadStrandDict:
            multiExonReadStrandDict[rName] = {"fwdCanonicalCounter": 0, "revCanonicalCounter": 0}
        if dinucleotideType == "left":
            flag = isCanonicalSite("+", "donor", infoList[1])
            if flag == 1: multiExonReadStrandDict[rName]["fwdCanonicalCounter"] += 1
            flag = isCanonicalSite("-", "accept", infoList[1])
            if flag == 1: multiExonReadStrandDict[rName]["revCanonicalCounter"] += 1
        if dinucleotideType == "right":
            flag = isCanonicalSite("+", "accept", infoList[1])
            if flag == 1: multiExonReadStrandDict[rName]["fwdCanonicalCounter"] += 1
            flag = isCanonicalSite("-", "donor", infoList[1])
            if flag == 1: multiExonReadStrandDict[rName]["revCanonicalCounter"] += 1

    multiExonAmbiStrandBedList = []
    for i in bedObj.reads:
        if i in singleExonStrandDict and len(set(singleExonStrandDict[i])) == 1:
            if bedObj.reads[i].strand == singleExonStrandDict[i][0]:
                print(bedObj.reads[i], file=strandConfirmedOut)
            bedObj.reads[i].strand = singleExonStrandDict[i][0]
            print(bedObj.reads[i], file=strandAdjustOut)
        elif i in multiExonReadStrandDict:
            if multiExonReadStrandDict[i]["fwdCanonicalCounter"] - multiExonReadStrandDict[i]["revCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "+":
                    print(bedObj.reads[i], file=strandConfirmedOut)
                bedObj.reads[i].strand = "+"
                print(bedObj.reads[i], file=strandAdjustOut)
            elif multiExonReadStrandDict[i]["revCanonicalCounter"] - multiExonReadStrandDict[i]["fwdCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "-":
                    print(bedObj.reads[i], file=strandConfirmedOut)
                bedObj.reads[i].strand = "-"
                print(bedObj.reads[i], file=strandAdjustOut)
            else:
                chrom, start, end, strand = bedObj.reads[i].chrom, bedObj.reads[i].chromStart, bedObj.reads[i].chromEnd, bedObj.reads[i].strand
                multiExonAmbiStrandBedList.append(" ".join(map(str, [chrom, start, end, i, ".", strand])))
        else:
            print(bedObj.reads[i], file=strandAdjustOut)
            print(bedObj.reads[i], file=strandConfirmedOut)
    multiExonAmbiStrandBedObj = pybedtools.bedtool.BedTool("\n".join(multiExonAmbiStrandBedList), from_string=True)
    multiExonAmbiStrandIntersectRes = gpeBedObj.intersect(multiExonAmbiStrandBedObj, wa=True, wb=True)
    multiExonAmbiStrandDict = {}
    for i in multiExonAmbiStrandIntersectRes:
        infoList = str(i).strip("\n").split("\t")
        if infoList[-3] not in multiExonAmbiStrandDict:
            multiExonAmbiStrandDict[infoList[-3]] = [Bed12(str(i)).strand]
        else:
            multiExonAmbiStrandDict[infoList[-3]].append(Bed12(str(i)).strand)
    for i in multiExonAmbiStrandBedList:
        infoList = i.split(" ")
        if infoList[3] in multiExonAmbiStrandDict and len(set(multiExonAmbiStrandDict[infoList[3]])) == 1:
            if bedObj.reads[infoList[3]].strand == multiExonAmbiStrandDict[infoList[3]][0]:
                print(bedObj.reads[infoList[3]], file=strandConfirmedOut)
            bedObj.reads[infoList[3]].strand = multiExonAmbiStrandDict[infoList[3]][0]
            print(bedObj.reads[infoList[3]], file=strandAdjustOut)
        else:
            print(bedObj.reads[infoList[3]], file=strandAdjustOut)
            print(bedObj.reads[infoList[3]], file=strandConfirmedOut)
    strandAdjustOut.close()
    strandConfirmedOut.close()
    removeFiles(os.getcwd(), [tmpBed])


def readsAssign(transBedFile, readsBedFile, offset=10, minConsenesusIntronN=1, minCoverageOnRead=0.9, transColNum=13,
                readsColNum=13, outPrefix="readsAssign", group=True):
    """Assign every read/isoform to the reference transcripts it overlaps on the same strand (legacy ``readsAssign``).

    ``<outPrefix>.bed12+``: the read, its type (``E`` exon-matched, ``I`` overlapping only, ``IG`` intergenic) and the
    matched transcripts / genes / consensus-intron counts / coverages. With ``group`` the ``.unambi.bed12+`` file adds
    the unique gene of ``E`` reads, or a ``chrom:strand:n`` novel cluster id for the rest."""
    if transColNum < 13:
        raise Exception("The colmun count of your transcript bed should large than 12")
    transBedObj = pybedtools.BedTool(transBedFile)
    readsBedObj = pybedtools.BedTool(readsBedFile)
    intersectRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    notintersectRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    matched = {}            # read name -> {"readsBed", "trans": [(transcript, gene, covOnRead, covOnTrans, extra)]}
    readClass = Bed12Plus if readsColNum > 12 else Bed12
    readKey, parsedTrans = None, {}     # parse each read once (its hits are consecutive rows) and each transcript once
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        if infoList[:readsColNum] != readKey:
            readKey = infoList[:readsColNum]
            readsBed = readClass("\t".join(readKey))
        transKey = "\t".join(infoList[readsColNum:readsColNum + transColNum])
        transBed = parsedTrans.get(transKey)
        if transBed is None:
            transBed = parsedTrans[transKey] = Bed12Plus(transKey)
        overlapLength = getBlockLength(getOverlapOfTuple(readsBed.exons, transBed.exons))
        coverageOnRead = overlapLength / float(getBlockLength(readsBed.exons))
        coverageOnTrans = overlapLength / float(getBlockLength(transBed.exons))
        if len(readsBed.exons) == 1:
            readStart, readEnd = readsBed.chromStart, readsBed.chromEnd
            inEndExon = int(transBed.exonStarts[0] < readEnd <= transBed.exonEnds[0]
                            or transBed.exonStarts[-1] <= readStart < transBed.exonEnds[-1])
            extra = (len(transBed.exons), inEndExon)
        else:
            consensusIntronN = getConsensusIntronN(readsBed.exonStarts, readsBed.exonEnds, transBed.exonStarts, transBed.exonEnds, offset=offset)
            readsJuncs, transJuncs = readsBed.juncChain.split(";"), transBed.juncChain.split(";")
            flank = int(transJuncs[-len(readsJuncs):] == readsJuncs or transJuncs[:len(readsJuncs)] == readsJuncs)
            extra = (consensusIntronN, flank)
        entry = matched.setdefault(readsBed.name, {"readsBed": readsBed, "trans": []})
        entry["trans"].append((transBed.name, transBed.record[12], coverageOnRead, coverageOnTrans, extra))

    outList = []
    for readName, entry in matched.items():
        readsBed, trans = entry["readsBed"], entry["trans"]
        if len(readsBed.exons) == 1:
            kept = [t for t in trans if not (t[4][0] > 1 and t[4][1] == 0 and t[2] < minCoverageOnRead)]
            consensus = "NA"
        else:
            readIntronN = len(readsBed.introns)
            kept = [t for t in trans if t[4][1] == 1
                    or (readIntronN < minConsenesusIntronN and readIntronN == t[4][0])
                    or (readIntronN >= minConsenesusIntronN and t[4][0] >= minConsenesusIntronN)
                    or t[2] >= minCoverageOnRead]
            consensus = None
        readType = "E" if kept else "I"
        rows = kept or trans                                       # "I" reads report every overlapping transcript
        outList.append(["\t".join(readsBed.record), readType, ",".join(t[0] for t in rows), ",".join(t[1] for t in rows),
                        consensus if consensus else ",".join(str(t[4][0]) for t in rows),
                        ",".join(str(t[2]) for t in rows), ",".join(str(t[3]) for t in rows)])

    for i in notintersectRes:
        outList.append([str(i).strip("\n"), "IG", "NA", "NA", "NA", "NA", "NA"])

    sortedOutList = sorted(outList, key=lambda x: (x[0].split("\t")[0], int(x[0].split("\t")[1])))
    with open(outPrefix + ".bed12+", "w") as assignOut:
        for item in sortedOutList:
            print("\t".join(map(str, item)), file=assignOut)

    if group:
        unambiOut = open(outPrefix + ".unambi.bed12+", "w")
        novelList = []
        for item in sortedOutList:
            if item[3] == "NA" or item[1] != "E":
                novelList.append(item[0].split("\t"))
            else:
                uniqGenes = list(set(item[3].split(",")))
                if len(uniqGenes) == 1:
                    print(item[0] + "\t" + uniqGenes[0], file=unambiOut)
                else:
                    novelList.append(item[0].split("\t"))
        sortedNovelList = sorted(novelList, key=lambda x: (x[0], x[5], int(x[1])))
        inc, clusterEnd, clusterKey = 0, 0, None
        for novelItem in sortedNovelList:
            myStart, myEnd = int(novelItem[1]), int(novelItem[2])
            if (novelItem[0], novelItem[5]) != clusterKey:          # new chromosome/strand: start over
                clusterKey, clusterEnd, inc = (novelItem[0], novelItem[5]), 0, 0
            if myStart < clusterEnd:
                if myEnd > clusterEnd:
                    clusterEnd = myEnd
            else:
                inc += 1
                clusterEnd = myEnd
            print("\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc])), file=unambiOut)
        unambiOut.close()


def filterJunctionFile(junctionFile, bedFile):
    import pandas as pd
    juncObj = pybedtools.BedTool(junctionFile)
    bedObj = pybedtools.BedTool(bedFile)
    intersectRes = juncObj.intersect(bedObj, wa=True, wb=True, s=True, loj=True)
    gene2junc = {}
    myList = []
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")

        geneName = infoList[24]
        if geneName == ".": continue
        if geneName not in gene2junc:
            myList.append(infoList[0:12] + [geneName])
            gene2junc[geneName] = {infoList[3]: infoList[0:12] + [geneName]}
        elif infoList[3] not in gene2junc[geneName]:
            myList.append(infoList[0:12] + [geneName])
            gene2junc[geneName].update({infoList[3]: infoList[0:12] + [geneName]})
    columns = ["chrom", "chromStart", "chromEnd", "iso", "coverage", "strand", "thickStart", "thickEnd", "rgb",
               "exonNum", "blockSizes", "blockStarts", "gene"]
    myDf = pd.DataFrame(myList, columns=columns)
    myDf.coverage = myDf.coverage.astype(int)
    # per gene, drop junctions covered by < 5 % of max(median, mean) coverage
    maxMiddle = myDf.groupby("gene").coverage.transform(lambda c: max(c.median(), c.mean()))
    myDf = myDf[~(myDf.coverage < maxMiddle * 0.05)]
    myDf.to_csv("filtered.junction.bed", index=False, header=False, sep="\t", columns=columns[:-1])
    return "filtered.junction.bed"


def _covered_by(intersectRes):
    """Junction name (column 4) -> names of the junctions it overlaps (column 10) in a ``-wa -wb`` intersection."""
    covered = {}
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        covered.setdefault(infoList[3], []).append(infoList[9])
    return covered


def refineWithJunc(bedFile=None, junctionFile=None, outFile=None, genomeFasta=None):
    """Replace each LR-only junction that overlaps no other LR junction by the NGS junctions covering it, in every read
    that has it. A read keeps its old structure when a covering junction is non-canonical or its exon length moves by
    more than 5 %, and is dropped when the new chain leaves an exon of length <= 0."""
    filteredJunctionFile = filterJunctionFile(junctionFile, bedFile)
    bedJunc = BedFile(bedFile, type="bed12+")
    junction = BedFile(filteredJunctionFile, type="bed12")
    bedJuncDict = bedJunc.getAllJuncDict(onebase=False, withStrand=True)
    bedJunc2trans = bedJunc.getJunc2trans(onebase=False, withStrand=True)
    juncDict = junction.getAllJuncDict(onebase=False, withStrand=True)
    allJunc2motif = getSpliceMotifFromJuncDict(dict(bedJuncDict, **juncDict), genomeFasta, withStrand=True)

    bedJuncObj = pybedtools.BedTool("\n".join(bedJuncDict.values()), from_string=True)
    juncObj = pybedtools.BedTool("\n".join(juncDict.values()), from_string=True)
    bedJuncIsCovered = _covered_by(bedJuncObj.intersect(juncObj, wa=True, wb=True, s=True))
    bedJuncSelfIsCovered = _covered_by(bedJuncObj.intersect(bedJuncObj, wa=True, wb=True, s=True))

    unvalidReads = set()
    for i in bedJuncDict:
        if i in juncDict or len(set(bedJuncSelfIsCovered[i])) != 1 or i not in bedJuncIsCovered:
            continue
        juncPattern = bedJuncDict[i].split("\t")[3].split(":")[1]
        covering = sorted(bedJuncIsCovered[i], key=lambda x: int(re.split(":|-", x)[1]))
        juncSub = ";".join(x.split(":")[1] for x in covering)
        canonical = all([allJunc2motif[j] in ("GT-AG", "GC-AG", "AT-AC") for j in covering])
        for read in bedJunc2trans[i]:
            readObj = bedJunc.reads[read]
            # BUG-182: replace the junction itself; a regex replace also rewrote other junctions containing its text
            # (100-200 inside 1100-2005 became 1100-1955)
            juncChain = ";".join(juncSub if junc == juncPattern else junc for junc in readObj.juncChain.split(";"))
            if not validJunction(readObj.chromStart, readObj.chromEnd, juncChain):
                unvalidReads.add(read)
                continue
            originReadObj = copy.copy(readObj)
            readObj.juncChain = juncChain
            readObj.getBedFromJuncChain(otherInfo="\t".join(readObj.otherList))
            if originReadObj.get_all_exons_len() * 0.95 > sum(readObj.blockSizes) or \
                    readObj.get_all_exons_len() * 1.05 < sum(readObj.blockSizes) or not canonical:
                bedJunc.reads[read] = originReadObj
            else:
                bedJunc.reads[read] = readObj

    with open(outFile, "w") as out:
        for read, readObj in bedJunc.reads.items():
            if read not in unvalidReads:
                print(readObj, file=out)


# ---------------------------------------------------------------------------
# Driver (LR core of legacy refineJunc, explicit paths)
# ---------------------------------------------------------------------------


def refine_isoforms(collapsed_gff, ref_genome, ref_gpe, ref_bed, out_dir, *,
                    ngs_junctions=None, adjust=True, runner=None):
    """Refine cupcake-collapsed isoforms into the grouped high-quality isoform set.

    Mirrors the LR core of legacy ``refineJunc`` (gpe→bed→strandAdjust→[refineWithJunc]→
    readsAssign→isoformGrouped→anno/novel). Runs inside ``<out_dir>/refine`` (the legacy
    stage cwd) so the ported relative temp files (tmp.bed, filtered.junction.bed) stay
    contained. ``ngs_junctions`` enables the (behavior-sensitive) NGS junction refinement;
    when None the TGS-only path is taken. Returns a dict of output paths.
    """
    collapsed_gff = os.path.abspath(collapsed_gff)
    ref_genome = os.path.abspath(ref_genome)
    ref_gpe = os.path.abspath(ref_gpe)
    ref_bed = os.path.abspath(ref_bed)
    if ngs_junctions:
        ngs_junctions = os.path.abspath(ngs_junctions)
    refine_dir = os.path.join(os.path.abspath(out_dir), "refine")
    runner = runner or ToolRunner()

    print(getCurrentTime() + " Refine the collapsed isoforms...")
    prev_dir = os.getcwd()
    resolveDir(refine_dir)   # mkdir + chdir, matching legacy stage behavior
    try:
        # cupcake GFF -> genePredExt -> bed12+ (gpe2bed.pl -g replaced by GenePredObj.toBed)
        runner.run(["gtfToGenePred", collapsed_gff, "tofu.collapsed.gpe", "-genePredExt"])
        with open("tofu.collapsed.bed12+", "w") as fh:
            for line in GenePredObj("tofu.collapsed.gpe", bincolumn=False).toBed(gene=True):
                fh.write(str(line) + "\n")

        if adjust:
            strandAdjust(ref_genome, ref_gpe, "tofu.collapsed.bed12+", 2, 0.8,
                         strandAdjust="tofu.strandAdjusted.bed12+",
                         strandConfirmed="tofu.strandConfirm.bed12+")
        else:
            makeLink("tofu.collapsed.bed12+", "tofu.strandConfirm.bed12+")

        if ngs_junctions:
            refineWithJunc(bedFile="tofu.strandConfirm.bed12+", junctionFile=ngs_junctions,
                           outFile="tofu.juncAdjusted.bed12+", genomeFasta=ref_genome)
            assignInput = "tofu.juncAdjusted.bed12+"
        elif adjust:
            assignInput = "tofu.strandConfirm.bed12+"
        else:
            assignInput = "tofu.collapsed.bed12+"

        readsAssign(ref_bed, assignInput, readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)

        # cut -f1-12,14 tofu.collapsed.assigned.unambi.bed12+ > isoformGrouped.bed12+
        with open("tofu.collapsed.assigned.unambi.bed12+") as src, open("isoformGrouped.bed12+", "w") as dst:
            for line in src:
                cols = line.rstrip("\n").split("\t")
                dst.write("\t".join(cols[:12] + [cols[13]]) + "\n")

        identifyNovelIsoformsByJunctions(ref_gpe, "isoformGrouped.bed12+",
                                         anno="isoformGrouped.anno.bed12+",
                                         novel="isoformGrouped.novel.bed12+")

        mapping_dir = os.path.join(os.path.abspath(out_dir), "mapping")
        reads_inputs = (os.path.join(mapping_dir, "flnc.processed.fa"),
                        os.path.join(mapping_dir, "flnc.addCVandID.bed12+"),
                        os.path.join(os.path.dirname(collapsed_gff), "tofu.collapsed.group.txt"))
        if all(map(os.path.isfile, reads_inputs)):
            refine_reads(ref_bed, *reads_inputs, runner=runner)
        else:
            print(getCurrentTime() + " Skip reads-level assignment (needs mapping/flnc.processed.fa, "
                  "mapping/flnc.addCVandID.bed12+ and collapse/tofu.collapsed.group.txt)")
    finally:
        os.chdir(prev_dir)

    print(getCurrentTime() + " Refine the collapsed isoforms done!")
    return {
        "isoformGrouped": os.path.join(refine_dir, "isoformGrouped.bed12+"),
        "anno": os.path.join(refine_dir, "isoformGrouped.anno.bed12+"),
        "novel": os.path.join(refine_dir, "isoformGrouped.novel.bed12+"),
        "assigned_unambi": os.path.join(refine_dir, "tofu.collapsed.assigned.unambi.bed12+"),
        "reads_assigned_unambi": os.path.join(refine_dir, "reads.assigned.unambi.bed12+"),
        "reads_gtf": os.path.join(refine_dir, "reads.unambi.gtf"),
    }


def refine_reads(ref_bed, processed_fa, processed_bed, group_file, *, runner=None):
    """Reads-level half of legacy ``refineJunc`` (runs in the refine dir): keep the FLNC reads that collapse
    used (``tofu.collapsed.group.txt``), assign them to reference genes -> ``reads.assigned.unambi.bed12+``
    (consumed by find_pa / palen_as / report) and export them as ``reads.unambi.gtf``."""
    runner = runner or ToolRunner()
    used = set()
    with open(group_file) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 1:
                used.update(cols[1].split(","))
    with open("collapsed.reads.lst", "w") as fh:
        fh.write("".join(name + "\n" for name in sorted(used)))
    runner.run(f"seqkit grep {q(processed_fa)} -f collapsed.reads.lst -w 0 > flnc.processed.ignore_id_removed.fa")
    with open(processed_bed) as src, open("flnc.processed.ignore_id_removed.bed12+", "w") as dst:   # filter.pl -2 4 -m i
        for line in src:
            if line.split("\t", 4)[3] in used:
                dst.write(line)
    readsAssign(ref_bed, "flnc.processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned", group=True)
    runner.run("cut -f 1-12,15 reads.assigned.unambi.bed12+ | {}/bed2gpe.pl -b 12 -g 13 - | "
               "genePredToGtf file stdin reads.unambi.gtf -source=iFLAS".format(q(_find_utils_dir())))
