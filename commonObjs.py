#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: commonObjs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from sys import maxint as MAXINT
from multiprocessing import Pool
import multiprocessing, pybedtools
import multiprocessing.pool

AStypes = ["IR", "SE", "A3SS", "A5SS"]

class Number:
    def __init__(self, count=0):
        self.count = count


class Bed6(object):
    "BED6 format gene structure"

    class BedError(Exception):
        "Error in manipulating bed12 structures"

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    def __init__(self, line=None):
        if line:
            self.record = line.strip().split("\t")
            (self.chrom, self.chromStart, self.chromEnd, self.name,
             self.score, self.strand) = self.record[:6]
            self.chromStart = int(self.chromStart)
            self.chromEnd = int(self.chromEnd)
            try:
                self.score = int(float(self.score))
            except ValueError:
                pass
            self.pos = (str(self.chrom) + ":" + str(self.chromStart) +
                        "-" + str(self.chromEnd))
        else:
            self.empty()

    def empty(self):
        (self.chrom, self.chromStart, self.chromEnd, self.name,
         self.score, self.strand) = ("", 0, 0, "", 0, "")

    def toBed12(self):
        line = "\t".join([repr(self), repr(self.chromStart),
                          repr(self.chromEnd), "255,0,0", "1",
                          repr(self.chromEnd - self.chromStart), "0"])
        return Bed12(line)

    def __repr__(self):
        "return a line of bed6 format, without newline ending"
        fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
                  self.name, str(self.score), self.strand]
        return "\t".join(fields)


class Bed12(object):
    "BED12 format gene structure."

    class BedError(Exception):
        "Error in manipulating Bed12 structures"

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    def __init__(self, line=None):
        if line:
            self.record = line.strip().split("\t")
            (self.chrom, self.chromStart, self.chromEnd, self.name,
             self.score, self.strand, self.thickStart, self.thickEnd,
             self.itemRgb, self.blockCount, self.blockSizes,
             self.blockStarts) = self.record[:12]

            self.chromStart = int(self.chromStart)
            self.chromEnd = int(self.chromEnd)
            self.score = int(float(self.score))
            self.thickStart = int(self.thickStart)
            self.thickEnd = int(self.thickEnd)
            self.blockCount = int(self.blockCount)

            self.blockSizes = self.blockSizes.strip(",").split(",")
            self.blockStarts = self.blockStarts.strip(",").split(",")

            assert len(self.blockStarts) == len(self.blockSizes)
            assert len(self.blockStarts) == self.blockCount

            for i in range(self.blockCount):
                self.blockSizes[i] = int(self.blockSizes[i])
                self.blockStarts[i] = int(self.blockStarts[i])

            self.exonStarts = []
            self.exonEnds = []
            for i in range(self.blockCount):
                self.exonStarts.append(self.chromStart + self.blockStarts[i])
                self.exonEnds.append(self.exonStarts[i] + self.blockSizes[i])

            self.exons = self.parse_exon()
            self.introns = self.parse_intron()
            self.exonChain = ';'.join(map(lambda x: str(self.exons[x][0])+"-"+str(self.exons[x][1]), range(len(self.exons))))
            self.juncChain = ";".join(map(lambda x: str(self.introns[x][0])+"-"+str(self.introns[x][1]), range(len(self.introns))))
        else:
            self.empty()

    def empty(self):
        "return an empty class with all values None, '', or []."
        self.chrom = ""
        self.chromStart = self.chromEnd = 0
        self.name = ""
        self.score = 0
        self.strand = ""
        self.thickStart = self.thickEnd = 0
        self.itemRgb = "0,0,0"
        self.blockCount = 0
        self.blockSizes = []
        self.blockStarts = []

    def __len__(self):
        "the length of transcript"
        return sum(self.blockSizes)

    def cds_len(self):
        "return cds/thick length"
        return sum([min(ed, self.thickEnd) - max(st, self.thickStart)
                    for (st, ed) in zip(self.exonStarts, self.exonEnds)
                    if (ed - self.thickStart) * (st - self.thickEnd) < 0])

    def parse_exon(self):
        "return a list of exon pos [(st, ed), (st, ed) , ... ]"
        exons = []
        for i in range(self.blockCount):
            st = self.exonStarts[i]
            ed = self.exonEnds[i]
            exons.append((st, ed))
        return exons

    def parse_intron(self):
        "return a list of intron pos [(st, ed], (st, ed], ... ]"
        introns = []
        for i in range(self.blockCount - 1):
            st = self.exonEnds[i]
            ed = self.exonStarts[i + 1]
            introns.append((st, ed))
        return introns

    def get_all_exons_len(self):
        return sum(map(lambda x: int(x[1]) - int(x[0]), self.exons))

    def get_all_introns_len(self):
        return sum(map(lambda x: int(x[1]) - int(x[0]), self.introns))

    def utr_5_len(self):
        "return the length of 5'UTR"
        if self.strand == "+":
            utr_5_len = sum([min(ed, self.thickStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.thickStart])
        else:
            utr_5_len = sum([ed - max(st, self.thickEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.thickEnd])
        return utr_5_len

    def utr_3_len(self):
        "return the length of 3'UTR"
        if self.strand == "-":
            utr_3_len = sum([min(ed, self.thickStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.thickStart])
        else:
            utr_3_len = sum([ed - max(st, self.thickEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.thickEnd])
        return utr_3_len

    def has_intron(self, intron):
        if intron.chrom == self.chrom and intron.strand == self.strand:
            for i in range(self.blockCount - 1):
                if (self.exonEnds[i] == intron.start and
                        self.exonStarts[i + 1] == intron.end):
                    return True
        return False

    def getBedFromJuncChain(self, otherInfo=""):
        juncStarts = map(int, [i.split("-")[0] for i in self.juncChain.split(";")])
        juncEnds = map(int, [i.split("-")[1] for i in self.juncChain.split(";")])
        self.blockSizes = [juncStarts[0] - self.chromStart]
        self.blockStarts = [0]
        for i in range(len(juncStarts)-1):
            self.blockSizes.append(juncStarts[i+1] - juncEnds[i])
        self.blockSizes.append(self.chromEnd - juncEnds[-1])

        for i in range(len(juncStarts)):
            self.blockStarts.append(juncEnds[i] - self.chromStart)

        self.blockCount = len(self.blockStarts)
        fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
                  self.name, str(self.score), self.strand, str(self.thickStart),
                  str(self.thickEnd), self.itemRgb, str(len(self.blockSizes)), ",".join(map(str, self.blockSizes)),
                  ",".join(map(str, self.blockStarts)), otherInfo]
        return "\t".join(fields)

    def __repr__(self):
        "return a line of bed12 format, without newline ending"
        fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
                  self.name, str(self.score), self.strand, str(self.thickStart),
                  str(self.thickEnd), self.itemRgb, str(self.blockCount)]

        blockSizesline = ','.join(repr(i) for i in self.blockSizes)
        blockStartsline = ','.join(repr(i) for i in self.blockStarts)
        fields.extend((blockSizesline, blockStartsline))
        return "\t".join(fields)

    def go_to_gpe(self):
        "return a gpe with same structure"
        gpe = GenePredExtLine()
        gpe.bin = None
        gpe.transName = self.name
        gpe.chrom = self.chrom
        gpe.strand = self.strand
        gpe.txStart = self.chromStart
        gpe.txEnd = self.chromEnd
        gpe.cdsStart = self.thickStart
        gpe.cdsEnd = self.thickEnd
        gpe.exonCount = self.blockCount
        gpe.exonStarts = self.exonStarts
        gpe.exonEnds = self.exonEnds
        gpe.score = self.score
        gpe.geneName = "."
        gpe.cdsStartStat = "."
        gpe.cdsEndStat = "."
        gpe.exonFrames = ["."]
        return gpe

    def cds_bed(self):
        "return a new Bed12 object of cds region of self"
        assert self.thickStart < self.thickEnd
        newbed = Bed12(self.__repr__())
        newBlockSizes = []
        newBlockStarts = []
        newBlockCount = newbed.blockCount
        for i in range(newbed.blockCount):
            if newbed.thickEnd < newbed.exonStarts[i]:
                newBlockCount -= 1
            elif (newbed.thickStart <= newbed.exonStarts[i] and
                  newbed.thickEnd <= newbed.exonEnds[i]):
                newBlockSizes.append(newbed.thickEnd - newbed.exonStarts[i])
                newBlockStarts.append(newbed.exonStarts[i] - newbed.thickStart)
            elif (newbed.thickStart >= newbed.exonStarts[i] and
                  newbed.thickEnd <= newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.exonStarts[i])
                newBlockStarts.append(0)
            elif (newbed.thickStart >= newbed.exonStarts[i] and
                  newbed.thickEnd > newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.thickStart)
                newBlockStarts.append(0)
            elif (newbed.thickStart < newbed.exonStarts[i] and
                  newbed.thickEnd > newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.exonStarts[i])
                newBlockStarts.append(newbed.exonStarts[i] - newbed.thickStart)
            elif newbed.thickStart > newbed.exonEnds[i]:
                newBlockCount -= 1
            else:
                raise self.BedError("Un-expected transcript structure.")
        assert len(newBlockSizes) == len(newBlockStarts) == newBlockCount
        newbed.blockSizes = newBlockSizes
        newbed.blockStarts = newBlockStarts
        newbed.blockCount = newBlockCount
        newbed.chromStart = newbed.thickStart
        newbed.chromEnd = newbed.thickEnd
        # renew exonStarts and exonEnds, in case further use
        return Bed12(newbed.__repr__())


class Bed6Plus(Bed6):
    def __init__(self, line=None):
        if line:
            Bed6.__init__(self, line)
            self.otherList = self.record[6:]

    def __str__(self):
        return Bed6.__str__(self) + "\t" + "\t".join(self.otherList)


class Bed12Plus(Bed12):
    def __init__(self, line=None):
        if line:
            Bed12.__init__(self, line)
            self.otherList = self.record[12:]

    def __str__(self):
        return Bed12.__str__(self) + "\t" + "\t".join(self.otherList)


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

    def getJuncChainDict(self):
        juncChainDict = {}
        annoSingleExonList = []
        for i in self.reads:
            if len(self.reads[i].exons) > 1:
                juncChainInfo = "{}:{}".format(self.reads[i].chrom, self.reads[i].juncChain)
                if juncChainInfo not in juncChainDict:
                    juncChainDict[juncChainInfo] = self.reads[i]
            # else:
            #     singleExon = "\t".join([self.reads[i].chrom, self.reads[i].chromStart, self.reads[i].chromEnd, self.reads[i].name, ".", self.reads[i].strand])
            #     annoSingleExonList.append(singleExon)
        # return juncDict, annoSingleExonList
        return juncChainDict

    def getAllJuncDict(self, onebase=True, withStrand=False):
        allJuncDict = {}
        for i in self.reads:
            if len(self.reads[i].exons) > 1:
                for junc in self.reads[i].introns:
                    if onebase:
                        juncName = "{}:{}-{}".format(self.reads[i].chrom, junc[0] + 1, junc[1])
                    else:
                        juncName = "{}:{}-{}".format(self.reads[i].chrom, junc[0], junc[1])
                    if withStrand:
                        juncName = "{}:{}".format(juncName, self.reads[i].strand)
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
                    if onebase:
                        juncName = "{}:{}-{}".format(self.reads[i].chrom, junc[0] + 1, junc[1])
                    else:
                        juncName = "{}:{}-{}".format(self.reads[i].chrom, junc[0], junc[1])
                    if withStrand:
                        juncName = "{}:{}".format(juncName, self.reads[i].strand)
                    if juncName not in junc2trans:
                        junc2trans[juncName] = [self.reads[i].name]
                    else:
                        junc2trans[juncName].append(self.reads[i].name)
        return junc2trans


    def getSpliceMotif(self, genomeFasta=None, withStrand=False):
        dinucleotideBedList = []
        for i in self.reads:
            bedObj = self.reads[i]
            for j in bedObj.juncChain.split(":")[1].split(";"):
                chrom = bedObj.chrom
                strand = bedObj.strand
                juncStart, juncEnd = map(int, j.split("-"))
                if not withStrand:
                    juncName = "{}:{}-{}".format(chrom, juncStart + 1, juncEnd)
                else:
                    juncName = "{}:{}-{}:{}".format(chrom, juncStart + 1, juncEnd, strand)
                leftDinucleotidePos = "\t".join(map(str, [chrom, juncStart, juncStart + 2, ":".join([juncName, "left"]), ".", strand]))
                rightDinucleotidePos = "\t".join(map(str, [chrom, juncEnd - 2, juncEnd, ":".join([juncName, "right"]), ".", strand]))
                dinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])
        dinucleotideBedObj = pybedtools.BedTool("\n".join(dinucleotideBedList), from_string=True)
        dinucleotideBedRes = dinucleotideBedObj.sequence(genomeFasta, name=True, tab=True, s=True)

        juncMotifDict = {}
        for i in str(open(dinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            # s[s.find("(") + 1:s.find(")")]
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
        self.asDict = dict.fromkeys(AStypes, {})

    def update(self, readStruc):
        self.chrom = readStruc.chrom
        self.strand = readStruc.strand
        self.minpos = min(self.minpos, readStruc.chromStart)
        self.maxpos = max(self.maxpos, readStruc.chromEnd)
        self.readNames.append(readStruc.name)
        self.reads.update({readStruc.name: readStruc})

    # def updateFromGene2Reads(self, gene2ReadsObj):
    #     self.chrom = gene2ReadsObj.chrom
    #     self.strand = gene2ReadsObj.strand
    #     self.minpos = min(self.minpos, gene2ReadsObj.minpos)
    #     self.maxpos = max(self.maxpos, gene2ReadsObj.maxpos)
    #     self.trans.update(gene2ReadsObj.trans)
    #     self.reads.update(gene2ReadsObj.reads)
    #     self.readNames.append(gene2ReadsObj.readNames)

    def __str__(self):
        return "%s:%s-%s (%s) (%s) (%d reads)" % \
               (self.chrom, self.minpos, self.maxpos, self.geneName, self.strand, len(self.reads))

    __repr__ = __str__


class GenePredExtLine(object):
    ' a line of gpe file'

    def __init__(self, line="", bincolumn=True):
        'initialize each field; attribute blockSizes and blockStarts as BED12.'
        if line:
            self.record = line.strip().split("\t")
            self.bin = None
            if bincolumn == True:
                self.bin = self.record.pop(0)
            self.transName = self.record[0]
            self.chrom = self.record[1]
            self.strand = self.record[2]
            self.txStart = int(self.record[3])
            self.txEnd = int(self.record[4])
            self.cdsStart = int(self.record[5])
            self.cdsEnd = int(self.record[6])
            self.exonCount = int(self.record[7])
            self.exonStarts = [int(i) for i in self.record[8].strip(',').split(',')]
            self.exonEnds = [int(i) for i in self.record[9].strip(',').split(',')]
            self.score = int(float(self.record[10]))
            self.geneName = self.record[11]
            self.cdsStartStat = self.record[12]
            self.cdsEndStat = self.record[13]
            self.exonFrames = [int(c) for c in
                               self.record[14].strip(",").split(",")]
            self.blockSizes = []
            self.blockStarts = []
            for i in range(self.exonCount):
                self.blockSizes.append(self.exonEnds[i] - self.exonStarts[i])
                self.blockStarts.append(self.exonStarts[i] - self.txStart)

            self.exons = self.parse_exon()
            self.introns = self.parse_intron()
        else:
            self.empty()

    def empty(self):
        "construct an empty gpe instance with all fields None, [], or 0"
        self.bin = None
        self.transName = ""
        self.chrom = ""
        self.strand = ""
        self.txStart = 0
        self.txEnd = 0
        self.cdsStart = 0
        self.cdsEnd = 0
        self.exonCount = 0
        self.exonStarts = []
        self.exonEnds = []
        self.score = 0
        self.geneName = ""
        self.cdsStartStat = ""
        self.cdsEndStat = ""
        self.exonFrames = []

    def __len__(self):
        "return total length of transcript"
        return sum([ed - st for st, ed in zip(self.exonStarts, self.exonEnds)])

    def copy(self):
        "return a new object of self"
        if self.bin:
            has_bin = True
        else:
            has_bin = False
        return GenePredExtLine(repr(self), bincolumn=has_bin)

    def cds_len(self):
        "return cds length"
        return sum([min(ed, self.cdsEnd) - max(st, self.cdsStart)
                    for (st, ed) in zip(self.exonStarts, self.exonEnds)
                    if (ed - self.cdsStart) * (st - self.cdsEnd) < 0])

    def utr_5_len(self):
        "return the length of 5'UTR"
        if self.strand == "+":
            utr_5_len = sum([min(ed, self.cdsStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.cdsStart])
        else:
            utr_5_len = sum([ed - max(st, self.cdsEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.cdsEnd])
        return utr_5_len

    def utr_3_len(self):
        "return the length of 3'UTR"
        assert self.is_standard()
        if self.strand == "-":
            utr_3_len = sum([min(ed, self.cdsStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.cdsStart])
        else:
            utr_3_len = sum([ed - max(st, self.cdsEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.cdsEnd])
        return utr_3_len

    def is_standard(self):
        "check if all fields in gpe are standard (i.e. no abnormal positions)"
        "the standards might be modified to accommodate specific filters"
        if (self.txStart < self.txEnd and
                self.cdsStart < self.cdsEnd and
                self.exonCount > 0):
            return True
        else:
            return False

    def is_complete(self):
        "return true if cdsStartStat and cdsEndStat are cmpl, else False"
        if self.cdsStartStat == self.cdsEndStat == "cmpl":
            return True
        else:
            return False

    def parse_exon(self):
        "return a list of exon pos [(st, ed), (st, ed) , ... ]"
        exons = []
        for i in range(self.exonCount):
            st = self.exonStarts[i]
            ed = self.exonEnds[i]
            exons.append((st, ed))
        return exons

    def parse_intron(self):
        "return a list of intron pos [(st, ed], (st, ed], ... ]"
        introns = []
        for i in range(self.exonCount - 1):
            st = self.exonEnds[i]
            ed = self.exonStarts[i + 1]
            introns.append((st, ed))
        return introns

    def has_intron(self, intron):
        "determine if the argument is one of self's introns"
        if self.chrom == intron.chrom and self.strand == intron.strand:
            for pos in self.introns:
                if pos[0] == intron.st and pos[1] == intron.ed:
                    return True
        return False

    def go_to_bed(self, plus=False, gene=False):
        "return a Bed12 object of same gene structure."
        if plus == True or gene == True:
            bed = Bed12Plus()
        else:
            bed = Bed12()
        bed.chrom = self.chrom
        bed.chromStart = self.txStart
        bed.chromEnd = self.txEnd
        bed.name = self.transName
        bed.score = self.score
        bed.itemRgb = "0,0,0"
        bed.strand = self.strand
        bed.thickStart = self.cdsStart
        bed.thickEnd = self.cdsEnd
        bed.blockCount = self.exonCount
        bed.blockSizes = self.blockSizes
        bed.blockStarts = self.blockStarts
        bed.exonStarts = self.exonStarts
        bed.exonEnds = self.exonEnds
        if gene == True:
            bed.otherList = [self.geneName]
        return bed

    def __repr__(self):
        "return the line generating this gpe object without the last newline."
        outputlist = [self.transName, self.chrom, self.strand, repr(self.txStart),
                      repr(self.txEnd), repr(self.cdsStart), repr(self.cdsEnd),
                      repr(self.exonCount)]
        exonStarts_seq = ",".join([repr(i) for i in self.exonStarts])
        exonEnds_seq = ",".join([repr(i) for i in self.exonEnds])
        outputlist.extend((exonStarts_seq, exonEnds_seq))
        exonFrames_seq = ",".join([repr(i) for i in self.exonFrames])
        outputlist.extend((repr(self.score), self.geneName, self.cdsStartStat,
                           self.cdsEndStat, exonFrames_seq))
        if self.bin:
            return self.bin + "\t" + "\t".join(outputlist)
        else:
            return "\t".join(outputlist)


class GenePredObj(object):
    """
    build gene prediction object from gpe file
    """
    def __init__(self, gpeFile, bincolumn=False):
        self.gpeFile = gpeFile
        self.bincolumn = bincolumn
        self.geneName2gpeObj = {}
        self.genePredDict = self.buildObj()

    def buildObj(self, geneId=None):
        tmpDict = {}
        with open(self.gpeFile) as f:
            for i in f:
                gpeObj = GenePredExtLine(i, bincolumn=self.bincolumn)
                chrom, strand, start, end = gpeObj.chrom, gpeObj.strand, gpeObj.txStart, gpeObj.txEnd
                if geneId != None:
                    gpeObj.transName = gpeObj.transName.replace(gpeObj.geneName, geneId)
                    gpeObj.geneName = geneId
                geneName = gpeObj.geneName

                if geneName not in self.geneName2gpeObj:
                    self.geneName2gpeObj[geneName] = [gpeObj]
                else:
                    self.geneName2gpeObj[geneName].append(gpeObj)

                if chrom not in tmpDict:
                    tmpDict[chrom] = {strand: [[start, end, gpeObj]]}
                elif strand not in tmpDict[chrom]:
                    tmpDict[chrom][strand] = [[start, end, gpeObj]]
                else:
                    tmpDict[chrom][strand].append([start, end, gpeObj])
            for k, v in tmpDict.iteritems():
                for s in v:
                    v[s] = sorted(v[s], key=lambda x: (x[0], x[1]))
        return tmpDict

    def changeGeneId(self, geneId):
        self.geneName2gpeObj = {}
        self.genePredDict = self.buildObj(geneId=geneId)

    def gpeObj2file(self, fout, geneName):
        f = open(fout, "w")
        gps = self.geneName2gpeObj[geneName]
        for i in gps:
            print >>f, i
        f.close()

    def getBlockLength(self, blockList):
        return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))

    def getGeneExonLength(self):
        exonLength = 0
        for g in self.geneName2gpeObj:
            exonList = []
            for l in self.geneName2gpeObj[g]:
                exonList.extend(l.exons)
            exonList = set(exonList)
            exonLength += self.getBlockLength(exonList)
        return exonLength

    def getGeneIntronLength(self):
        intronLength = 0
        for g in self.geneName2gpeObj:
            intronList = []
            for l in self.geneName2gpeObj[g]:
                intronList.extend(l.introns)
            intronList = set(intronList)
            intronLength += self.getBlockLength(intronList)
        return intronLength

    def toBed(self, plus=False, gene=False, outFile=None):
        if outFile:
            out = open(outFile, "w")
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bed = gp[2].go_to_bed()
                        print >>out, bed
            out.close()
        else:
            bedList = []
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bedList.append(str(gp[2].go_to_bed(gene=gene)))
            return bedList


class Junction(Bed12):
    """
    junction generated by tophat, with 5' donor and 3' acceptor attribute,
    (donor, acceptor) resembles the intron inside junction.
    """

    def __init__(self, line):
        Bed12.__init__(self, line)
        self.donor = self.chromStart + self.blockSizes[0]
        self.acceptor = self.chromStart + self.blockStarts[1]
        self.jPos = "%s:%d-%d" % (self.chrom, self.donor, self.acceptor)
        self.dn = self.donor if self.strand == "+" else self.acceptor
        self.ac = self.acceptor if self.strand == "+" else self.donor
        self.dnPos = "%s:%d:%s" % (self.chrom, self.dn, self.strand)


class NoDaemonProcess(multiprocessing.Process):
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


class MergedTgsSample(object):
    def __init__(self):
        self.project_name = None
        self.sample_name = None
        self.ref_strain = None
        self.strain = None
        # self.condition = None

        self.tgs_plat = None
        self.strategy = None
        self.data_location = None
        self.data_processed_location = None
        # self.primer = None
        self.polya_location = None

        self.ngs_left_reads = None
        self.ngs_right_reads = None
        self.ngs_reads_paired = "paired"
        self.ngs_reads_length = None
        self.ngs_junctions = None

        self.use_fmlrc2 = True
        self.single_run_threads = 1

    def __str__(self):
        return "%s:%s:%s, %s:%s:%s, %s:%s:%s:%s:%s, %d" % \
               (self.project_name, self.sample_name, self.ref_strain, self.tgs_plat, self.strategy,
                self.data_processed_location, self.ngs_left_reads, self.ngs_right_reads, self.ngs_reads_paired,
                self.ngs_reads_length, self.ngs_junctions, self.single_run_threads)

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

    def getExonFramesForGPE(self, gpeObj):
        exonFrame, exonFrames = 0, []
        if gpeObj.strand == "+":
            for i in range(gpeObj.exonCount):
                if gpeObj.exonStarts[i] < gpeObj.cdsStart:
                    exonFrame = -1 if gpeObj.exonEnds[i] < gpeObj.cdsStart else 0
                elif gpeObj.exonStarts[i] == gpeObj.cdsStart:
                    exonFrame = 0
                else:
                    if gpeObj.exonStarts[i] < gpeObj.cdsEnd:
                        if gpeObj.exonStarts[i - 1] < gpeObj.cdsStart:
                            exonFrame = (gpeObj.exonEnds[i - 1] - gpeObj.cdsStart) % 3
                        else:
                            exonFrame = (gpeObj.exonEnds[i - 1] - gpeObj.exonStarts[i - 1] + exonFrames[i - 1] - 3) % 3
                    else:
                        exonFrame = -1
                exonFrames.append(exonFrame)
        else:
            for i in range(gpeObj.exonCount - 1, -1, -1):
                if gpeObj.exonEnds[i] > gpeObj.cdsEnd:
                    exonFrame = -1 if gpeObj.exonStarts[i] > gpeObj.cdsEnd else 0
                elif gpeObj.exonEnds[i] == gpeObj.cdsEnd:
                    exonFrame = -1 if gpeObj.cdsStart == gpeObj.cdsEnd else 0
                else:
                    if gpeObj.exonEnds[i] > gpeObj.cdsStart:
                        if gpeObj.exonEnds[i + 1] > gpeObj.cdsEnd:
                            exonFrame = (gpeObj.cdsEnd - gpeObj.exonStarts[i + 1]) % 3
                        else:
                            exonFrame = (gpeObj.exonEnds[i + 1] - gpeObj.exonStarts[i + 1] + exonFrames[
                                gpeObj.exonCount - i - 2] - 3) % 3
                    else:
                        exonFrame = -1
                exonFrames.append(exonFrame)
        return exonFrames

    def go_to_gpe(self):
        gpe = GenePredExtLine()
        gpe.bin = None
        gpe.transName = self.name
        gpe.chrom = self.chrom
        gpe.strand = self.strand
        gpe.txStart = self.chromStart
        gpe.txEnd = self.chromEnd
        gpe.cdsStart = self.thickStart
        gpe.cdsEnd = self.thickEnd
        gpe.exonCount = self.blockCount
        gpe.exonStarts = self.exonStarts
        gpe.exonEnds = self.exonEnds
        gpe.score = 0
        gpe.geneName = self.geneName
        gpe.cdsStartStat = "unk"
        gpe.cdsEndStat = "unk"
        gpe.exonFrames = self.getExonFramesForGPE(gpe)
        return gpe

    def __str__(self):
        return "%s:%s:%d-%d (%s)" % (self.name, self.chrom, self.chromStart, self.chromEnd, self.geneName)

    __repr__ = __str__


