"""GenePred (GPE) parsers — faithfully ported from legacy commonObjs.py (behavior preserved).

Classes: GenePredExtLine (one gpe line), GenePredObj (build from a gpe file).
Depends only on the BED parsers in .intervals; no heavy imports.
"""

from __future__ import annotations

from .intervals import Bed12, Bed12Plus


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
            for k, v in tmpDict.items():
                for s in v:
                    v[s] = sorted(v[s], key=lambda x: (x[0], x[1]))
        return tmpDict



    def getBlockLength(self, blockList):
        return sum([int(x[1]) - int(x[0]) for x in blockList])



    def toBed(self, plus=False, gene=False, outFile=None):
        if outFile:
            out = open(outFile, "w")
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bed = gp[2].go_to_bed()
                        print(bed, file=out)
            out.close()
        else:
            bedList = []
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bedList.append(str(gp[2].go_to_bed(gene=gene)))
            return bedList
