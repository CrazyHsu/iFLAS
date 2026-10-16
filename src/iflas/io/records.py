"""Read/record structures — faithfully ported from legacy commonObjs.py.

Classes: BedFile (load BED into a name->record dict + junction/exon views),
Gene2Reads (per-gene read aggregation), ReadLineStruc (a BED12+ read line).
pybedtools is imported lazily (only BedFile.getSpliceMotif needs it).
"""

from __future__ import annotations

import sys

from ..constants import AStypes
from .intervals import Bed6, Bed6Plus, Bed12, Bed12Plus

MAXINT = sys.maxsize


def junc_name(chrom, start, end, onebase=True, strand=None):
    """``chrom:start-end[:strand]`` key of an intron; ``onebase`` reports the 1-based start used by the NGS junction
    files, otherwise the 0-based BED coordinate."""
    name = "{}:{}-{}".format(chrom, start + 1 if onebase else start, end)
    return "{}:{}".format(name, strand) if strand else name


class BedFile(object):
    def __init__(self, bedFile, type=None):
        self.bedFile = bedFile
        self.reads = self.getReadsInfo(type)

    def getReadsInfo(self, bedType):
        readsDict = {}
        with open(self.bedFile) as f:
            for line in f:
                if bedType == "bed12":
                    b = Bed12(line)
                    readsDict.__setitem__(b.name, b)
                elif bedType == "bed12+":
                    b = Bed12Plus(line)
                    readsDict.__setitem__(b.name, b)
                elif bedType == "bed6":
                    b = Bed6(line)
                    readsDict.__setitem__(b.name, b)
                elif bedType == "bed6+":
                    b = Bed6Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict

    def getGenePos(self, bedType=None, geneCol=13):
        genePos = {}
        if bedType == "bed12+":
            for r in self.reads:
                gene = self.reads[r].otherList[geneCol - 13]
                if gene not in genePos:
                    genePos[gene] = [self.reads[r].chrom, self.reads[r].thickStart, self.reads[r].thickEnd]
                else:
                    if self.reads[r].thickStart < genePos[gene][1] or self.reads[r].thickEnd > genePos[gene][2]:
                        genePos[gene] = [self.reads[r].chrom, self.reads[r].thickStart, self.reads[r].thickEnd]

            return genePos

    def getJuncChainDict(self, withStrand=False):
        """``chrom:chain`` (``chrom:chain:strand`` with ``withStrand``) -> first record with that intron chain."""
        juncChainDict = {}
        for i in self.reads:
            if len(self.reads[i].exons) > 1:
                juncChainInfo = "{}:{}".format(self.reads[i].chrom, self.reads[i].juncChain)
                if withStrand:
                    juncChainInfo += ":" + self.reads[i].strand
                if juncChainInfo not in juncChainDict:
                    juncChainDict[juncChainInfo] = self.reads[i]
        return juncChainDict

    def getAllJuncDict(self, onebase=True, withStrand=False):
        allJuncDict = {}
        for i in self.reads:
            if len(self.reads[i].exons) > 1:
                for junc in self.reads[i].introns:
                    juncName = junc_name(self.reads[i].chrom, junc[0], junc[1], onebase, self.reads[i].strand if withStrand else None)
                    if juncName not in allJuncDict:
                        allJuncDict[juncName] = "\t".join(map(str, [self.reads[i].chrom, junc[0], junc[1], juncName, self.reads[i].score, self.reads[i].strand]))
        return allJuncDict

    def getAllExonDict(self):
        allExonDict = {}
        for i in self.reads:
            for exon in self.reads[i].exons:
                exonName = "{}:{}-{}".format(self.reads[i].chrom, exon[0], exon[1])
                if exonName not in allExonDict:
                    allExonDict[exonName] = "\t".join(map(str, [self.reads[i].chrom, exon[0], exon[1], exonName, self.reads[i].score, self.reads[i].strand]))
        return allExonDict

    def getJunc2trans(self, onebase=True, withStrand=False):
        junc2trans = {}
        for i in self.reads:
            if len(self.reads[i].exons) > 1:
                for junc in self.reads[i].introns:
                    juncName = junc_name(self.reads[i].chrom, junc[0], junc[1], onebase, self.reads[i].strand if withStrand else None)
                    if juncName not in junc2trans:
                        junc2trans[juncName] = [self.reads[i].name]
                    else:
                        junc2trans[juncName].append(self.reads[i].name)
        return junc2trans

    def getSpliceMotif(self, genomeFasta=None, withStrand=False):
        import pybedtools  # lazy: only this method needs pybedtools
        dinucleotideBedList = []
        for i in self.reads:
            bedObj = self.reads[i]
            for j in bedObj.juncChain.split(":")[1].split(";"):
                chrom = bedObj.chrom
                strand = bedObj.strand
                juncStart, juncEnd = list(map(int, j.split("-")))
                juncName = junc_name(chrom, juncStart, juncEnd, True, strand if withStrand else None)
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


class Gene2Reads(object):
    def __init__(self, geneName):
        self.minpos = MAXINT
        self.maxpos = 0
        self.geneName = geneName
        self.readNames = []
        self.trans = {}
        self.reads = {}
        self.asDict = {as_type: {} for as_type in AStypes}

    def update(self, readStruc):
        self.chrom = readStruc.chrom
        self.strand = readStruc.strand
        self.minpos = min(self.minpos, readStruc.chromStart)
        self.maxpos = max(self.maxpos, readStruc.chromEnd)
        self.readNames.append(readStruc.name)
        self.reads.update({readStruc.name: readStruc})

    def __str__(self):
        return "%s:%s-%s (%s) (%s) (%d reads)" % \
               (self.chrom, self.minpos, self.maxpos, self.geneName, self.strand, len(self.reads))

    __repr__ = __str__


class ReadLineStruc(Bed12):
    """
    get Bed12+ read structure
    """
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.record = line.strip().split("\t")
        self.exonLen = sum(self.blockSizes)
        self.readLen = abs(self.chromStart - self.chromEnd)
        self.geneName = self.record[12]



    def __str__(self):
        return "%s:%s:%d-%d (%s)" % (self.name, self.chrom, self.chromStart, self.chromEnd, self.geneName)

    __repr__ = __str__
