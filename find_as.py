#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: find_as.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import pybedtools, itertools
from commonFuncs import *
from commonObjs import *

def getCnsSite(intronStarts, transIntronStarts, offset):
    '''
    get consensus intron boundary site between reference and reads info
    :return: count of consensus site num
    '''
    s2, cnsN = 0, 0
    for intronStart1 in intronStarts:
        for i in range(s2, len(transIntronStarts)):
            intronStart2 = transIntronStarts[i]
            if intronStart2 > intronStart1 + offset: break
            if intronStart1 - offset <= intronStart2 and intronStart2 <= intronStart1 + offset:
                cnsN += 1
                s2 += 1
    return cnsN


def getConsensusIntronNfromAS(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    totalCnsN = 0
    for asType in AStypes:
        if asType in ["IR", "SE"]:
            j, consensusN = 0, 0
            for i in range(len(intronStarts1)):
                for k in range(j, len(intronStarts2)):
                    if intronStarts2[k] > intronEnds1[i]: break
                    if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronStarts1[i] + offset and \
                            intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                        consensusN += 1
                        j += 1
            totalCnsN += consensusN
        elif asType in ["A3SS", "A5SS"]:
            cnsIntronStartN = getCnsSite(intronStarts1, intronStarts2, offset)
            cnsIntronEndN = getCnsSite(intronEnds1, intronEnds2, offset)
            totalCnsN += cnsIntronStartN + cnsIntronEndN
    return totalCnsN


def filterIrByJunc(irBed, junctionBed, confidentIrBed):
    juncDict = {}
    with open(junctionBed) as f:
        for line in f.readlines():
            juncDict.update({Junction(line.strip("\n")).jPos: ""})

    out = open(confidentIrBed, "w")
    with open(irBed) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            juncPos = "{}:{}-{}".format(infoList[0], infoList[1], infoList[2])
            if juncPos in juncDict:
                print >> out, line.strip("\n")
    out.close()


def getAnnotation(annoBedRes, novelBedRes, offset=0, reference=0):
    annoDict, novelDict = {}, {}
    isNovelDict = {}
    for line in annoBedRes:
        record = str(line).strip("\n").split("\t")
        read = Bed12("\t".join(record[:12]))
        trans = ReadLineStruc("\t".join(record[12:]))
        if reference and record[12] != record[-1]: continue
        if read.name not in isNovelDict:
            isNovelDict[read.name] = {"read": read, "isNovel": 1}
        transStarts, transEnds = trans.exonStarts, trans.exonEnds
        readStarts, readEnds, blockCount = read.exonStarts, read.exonEnds, read.blockCount
        consensusIntronN = getConsensusIntronNfromAS(readStarts, readEnds, transStarts, transEnds, offset)
        if consensusIntronN >= 1:
            transName, geneName = trans.name, trans.geneName
            if geneName not in annoDict:
                gene2reads = Gene2Reads(geneName)
                gene2reads.update(read)
                gene2reads.trans = {transName: trans}
                annoDict[geneName] = gene2reads
            else:
                gene2reads = annoDict[geneName]
                gene2reads.update(read)
                gene2reads.trans.update({transName: trans})
                annoDict[geneName] = gene2reads
            isNovelDict[read.name]["isNovel"] = 0
        else:
            if blockCount == 1:
                transName, geneName = trans.name, trans.geneName
                if isOverlap((read.chromStart, read.chromEnd), (trans.chromStart, trans.chromEnd)):
                    overlap = getOverlapOfTuple([(read.chromStart, read.chromEnd)], [(trans.chromStart, trans.chromEnd)])
                    sumOverlap = sum([x[1] - x[0]for x in overlap])
                    if sumOverlap/float(read.chromEnd-read.chromStart) >= 0.5:
                        if geneName not in annoDict:
                            gene2reads = Gene2Reads(geneName)
                            gene2reads.update(read)
                            gene2reads.trans = {transName: trans}
                            annoDict[geneName] = gene2reads
                        else:
                            gene2reads = annoDict[geneName]
                            gene2reads.update(read)
                            gene2reads.trans.update({transName: trans})
                            annoDict[geneName] = gene2reads
                        isNovelDict[read.name]["isNovel"] = 0

    for line in novelBedRes:
        read = Bed12(str(line))
        if read.chrom not in novelDict:
            novelDict[read.chrom] = {read.strand: [read]}
        elif read.strand not in novelDict[read.chrom]:
            novelDict[read.chrom][read.strand] = [read]
        else:
            novelDict[read.chrom][read.strand].append(read)

    for i in isNovelDict:
        if isNovelDict[i]["isNovel"]:
            read = isNovelDict[i]["read"]
            if read.chrom not in novelDict:
                novelDict[read.chrom] = {read.strand: [read]}
            elif read.strand not in novelDict[read.chrom]:
                novelDict[read.chrom][read.strand] = [read]
            else:
                novelDict[read.chrom][read.strand].append(read)
    return annoDict, novelDict


def findIR(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        irDict = {}
        for rName1 in rNames:
            exonStarts, exonEnds = readsDict[rName1].exonStarts, readsDict[rName1].exonEnds
            for exon1 in range(len(exonStarts) - 1):
                IR = "%d-%d" % (exonEnds[exon1], exonStarts[exon1 + 1])
                if IR in irDict: continue
                irDict[IR] = {}
                if "spliced" not in irDict[IR]:
                    irDict[IR].__setitem__("spliced", [readsDict[rName1]])
                else:
                    irDict[IR]["spliced"].append(readsDict[rName1])
                for rName2 in rNames:
                    if rName1 == rName2: continue
                    exons2 = readsDict[rName2].exonChain.split(";")
                    for exon2 in range(len(exons2)):
                        exon2Start, exon2End = [int(x) for x in exons2[exon2].split("-")]
                        if exon2Start < exonEnds[exon1] and exonStarts[exon1 + 1] < exon2End \
                                and abs(exon2Start - exonEnds[exon1]) > offset \
                                and abs(exonStarts[exon1 + 1] - exon2End) > offset:
                            if "retention" not in irDict[IR]:
                                irDict[IR].__setitem__("retention", [readsDict[rName2]])
                                break
                            else:
                                irDict[IR]["retention"].append(readsDict[rName2])
                        if exon2End == exonEnds[exon1]:
                            if exon2 + 1 >= len(exons2): break
                            exon2NStart, exon2NEnd = [int(x) for x in exons2[exon2 + 1].split("-")]
                            if exon2NStart == exonStarts[exon1 + 1]:
                                irDict[IR]["spliced"].append(readsDict[rName2])
                        if exon2Start > exonEnds[exon1 + 1]: break
        newIrDict = {}
        for ir in irDict:
            if "retention" in irDict[ir]:
                newIrDict[ir] = irDict[ir]
                if outAS:
                    irStart, irEnd = [int(x) for x in ir.split("-")]
                    if anno:
                        transDict = gene2ReadsDict[gName].trans
                        overlapWithGene = 0
                        for tran in transDict:
                            if irStart < transDict[tran].chromEnd and irEnd > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            retentionReads = [x.name for x in irDict[ir]["retention"]]
                            splicedReads = [x.name for x in irDict[ir]["spliced"]]
                            psi = float(len(retentionReads))/(len(retentionReads)+len(splicedReads))*1000
                            if isoform2reads:
                                rawRetentionReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in retentionReads]))
                                rawSplicedReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in splicedReads]))
                                psi = float(len(rawRetentionReads))/(len(rawRetentionReads)+len(rawSplicedReads))*1000
                            print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName+":"+ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
                    else:
                        retentionReads = [x.name for x in irDict[ir]["retention"]]
                        splicedReads = [x.name for x in irDict[ir]["spliced"]]
                        psi = float(len(retentionReads)) / (len(retentionReads) + len(splicedReads)) * 1000
                        if isoform2reads:
                            rawRetentionReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in retentionReads]))
                            rawSplicedReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in splicedReads]))
                            psi = float(len(rawRetentionReads))/(len(rawRetentionReads)+len(rawSplicedReads))*1000
                        print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName + ":" + ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
        gene2ReadsDict[gName].asDict.update({"IR": newIrDict})


def findSE(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None, restrict=False):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        exonDict = {}
        for rName1 in rNames:
            exons = readsDict[rName1].exonChain.split(";")
            for l in range(len(exons) - 2):
                for r in range(l + 2, len(exons)):
                    lExonEnd = int(exons[l].split("-")[1])
                    rExonStart = int(exons[r].split("-")[0])
                    SEs = ";".join(exons[l + 1:r])
                    seChain = "%d@%s@%d" % (lExonEnd, SEs, rExonStart)
                    if seChain in exonDict: continue
                    exonDict[seChain] = {}
                    if "keep" not in exonDict[seChain]:
                        exonDict[seChain].__setitem__("keep", [readsDict[rName1]])
                    else:
                        exonDict[seChain]["keep"].append(readsDict[rName1])
                    for rName2 in rNames:
                        if rName2 == rName1: continue
                        skipPat = "-%d;%d-" % (lExonEnd, rExonStart)
                        keepPat = "-%d;%s;%d-" % (lExonEnd, SEs, rExonStart)
                        if re.search(skipPat, readsDict[rName2].exonChain):
                            if "skip" not in exonDict[seChain]:
                                exonDict[seChain].__setitem__("skip", [readsDict[rName2]])
                            else:
                                exonDict[seChain]["skip"].append(readsDict[rName2])
                        elif re.search(keepPat, readsDict[rName2].exonChain):
                            exonDict[seChain]["keep"].append(readsDict[rName2])
        newExonDict = {}
        for seChain1 in exonDict.keys():
            if "skip" in exonDict[seChain1]:
                newExonDict[seChain1] = exonDict[seChain1]
                if outAS:
                    keepReads = [x.name for x in exonDict[seChain1]["keep"]]
                    skipReads = [x.name for x in exonDict[seChain1]["skip"]]
                    keepReadsCount = len(keepReads)
                    skipReadsCount = len(skipReads)
                    lSpliceSite, SEs1, rSpliceSite = seChain1.split("@")
                    SElist = SEs1.split(";")
                    starts, ends = [], []
                    for exon in SElist:
                        exonStart, exonEnd = listStr2Int(exon.split("-"))
                        starts.append(exonStart)
                        ends.append(exonEnd)
                    if anno:
                        overlapWithGene = 0
                        transDict = gene2ReadsDict[gName].trans
                        for tran in transDict:
                            if starts[0] < transDict[tran].chromEnd and ends[-1] > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                            blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                            psi = float(keepReadsCount)/(keepReadsCount + skipReadsCount) * 1000
                            if isoform2reads:
                                rawKeepReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in keepReads]))
                                rawSkipReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in skipReads]))
                                psi = float(len(rawKeepReads)) / (len(rawKeepReads) + len(rawSkipReads)) * 1000
                            print >> out, "\t".join(map(str, [chrom, str(starts[0]), str(ends[-1]),
                                                              "%s:%s" % (gName, seChain1), psi, strand, starts[0],
                                                              ends[-1], "0,0,0", len(starts), blockSizes,
                                                              blockRelStarts, gName, lSpliceSite, rSpliceSite,
                                                              ",".join(keepReads), len(keepReads), ",".join(skipReads),
                                                              len(skipReads)]))
                    else:
                        blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                        blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                        psi = float(keepReadsCount) / (keepReadsCount + skipReadsCount) * 1000
                        if isoform2reads:
                            rawKeepReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in keepReads]))
                            rawSkipReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in skipReads]))
                            psi = float(len(rawKeepReads)) / (len(rawKeepReads) + len(rawSkipReads)) * 1000
                        print >> out, "\t".join(map(str,
                                                    [chrom, str(starts[0]), str(ends[-1]), "%s:%s" % (gName, seChain1),
                                                     psi, strand, starts[0], ends[-1], "0,0,0", len(starts), blockSizes,
                                                     blockRelStarts, gName, lSpliceSite, rSpliceSite,
                                                     ",".join(keepReads), len(keepReads), ",".join(skipReads),
                                                     len(skipReads)]))
        gene2ReadsDict[gName].asDict.update({"SE": newExonDict})


def findA3SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None, restrict=False):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 1 if strand == "+" else 0
        for r1 in range(len(rNames) - 1):
            exonsChain1 = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain1.split(";")
            indexEnd1 = len(exons1) if strand == "+" else len(exons1) - 1
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].exonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) if strand == "+" else len(exons2) - 1
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "-":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon2End - int(exons1[i + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon1End - int(exons2[j + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon1Start - int(exons2[j - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon2Start - int(exons1[i - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            if outAS:
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
                incJunc = altDict[alt]["incJunc"]
                excJunc = altDict[alt]["excJunc"]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                if isoform2reads:
                    rawIncReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in incReads]))
                    rawExcReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in excReads]))
                    psi = float(len(rawIncReads)) / (len(rawIncReads) + len(rawExcReads)) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount), \
                                          "inc:{}-{}".format(incJunc[0], incJunc[1]), \
                                          "exc:{}-{}".format(excJunc[0], excJunc[1])]))

        gene2ReadsDict[gName].asDict.update({"A3SS": altDict})


def findA5SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None, restrict=False):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 0 if strand == "+" else 1
        for r1 in range(len(rNames) - 1):
            exonsChain1 = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain1.split(";")
            indexEnd1 = len(exons1) - 1 if strand == "+" else len(exons1)
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].exonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) - 1 if strand == "+" else len(exons2)
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "+":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon2End - int(exons1[i + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon1End - int(exons2[j + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon1Start - int(exons2[j - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon2Start - int(exons1[i - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            incJunc = altDict[alt]["incJunc"]
            excJunc = altDict[alt]["excJunc"]
            if outAS:
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                if isoform2reads:
                    rawIncReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in incReads]))
                    rawExcReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in excReads]))
                    psi = float(len(rawIncReads)) / (len(rawIncReads) + len(rawExcReads)) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount), \
                                          "inc:{}-{}".format(incJunc[0], incJunc[1]), \
                                          "exc:{}-{}".format(excJunc[0], excJunc[1])]))
        gene2ReadsDict[gName].asDict.update({"A5SS": altDict})


def find_single_as(asType="IR", gene2ReadsDict=None, anno=False, outFile=None, outAppend=False, isoform2reads=None):
    if outAppend:
        out = open(outFile, "a")
    else:
        out = open(outFile, "w")

    findASType = {"IR": findIR, "SE": findSE, "A3SS": findA3SS, "A5SS": findA5SS}
    if anno:
        if asType:
            if asType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[asType](gene2ReadsDict, asType, out, offset=5, isoform2reads=isoform2reads)
        else:
            for asType in findASType:
                findASType[asType](gene2ReadsDict, offset=5, isoform2reads=isoform2reads)
    else:
        newGene2ReadsDict = {}
        for myChr, myValue in gene2ReadsDict.iteritems():
            for myStrand in myValue:
                sortedStrandV = sorted(myValue[myStrand], key=lambda x: x.chromStart)
                count = 1
                fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                gene2reads = Gene2Reads(fakeGeneName)
                gene2reads.update(sortedStrandV[0])
                clusterEnd = sortedStrandV[0].chromEnd
                for i in range(1, len(sortedStrandV)):
                    myStart, myEnd = sortedStrandV[i].chromStart, sortedStrandV[i].chromEnd
                    if myStart < clusterEnd:
                        gene2reads.update(sortedStrandV[i])
                        if myEnd > clusterEnd:
                            clusterEnd = myEnd
                    else:
                        newGene2ReadsDict[fakeGeneName] = gene2reads
                        count += 1
                        fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                        gene2reads = Gene2Reads(fakeGeneName)
                        gene2reads.update(sortedStrandV[i])
                        clusterEnd = sortedStrandV[i].chromEnd
                newGene2ReadsDict[fakeGeneName] = gene2reads
        if asType:
            if asType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[asType](newGene2ReadsDict, asType, out, offset=5, isoform2reads=isoform2reads, anno=False)
        else:
            for asType in findASType:
                findASType[asType](newGene2ReadsDict, offset=5, isoform2reads=isoform2reads)

    out.close()

    # out = open(outFile, "w")
    # findAS(annoDict, outASType=asType, anno=True, out=out, isoform2reads=isoform2reads)
    # findAS(novelDict, outASType=asType, anno=False, out=out, isoform2reads=isoform2reads)
    # out.close()


def find_all_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)

    workDir = os.path.join(baseDir, "as_events", "ordinary_as")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(workDir)
    # isoformBed = os.path.join(baseDir, "refine", "isoformGrouped.bed12+")
    isoformBed = os.path.join(baseDir, "hqIsoforms", "hq.collapsed.bed12+")
    tofuGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    if dataObj.ngs_junctions == None and (dataObj.ngs_left_reads or dataObj.ngs_right_reads):
        dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")

    # find_python(isoformBed, tofuGroupFile, dataObj=dataObj, refParams=refParams)
    if not os.path.exists("LR"):
        os.makedirs("LR")
    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    tmpReads = {}
    isoform2reads = {}
    with open(tofuGroupFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            tmpReads[isoform] = reads.split(",")
    with open(isoformBed) as f:
        for line in f.readlines():
            readStruc = Bed12Plus(line)
            isoform2reads[readStruc.name] = list(itertools.chain.from_iterable([tmpReads[i] for i in readStruc.otherList[2].split(",")]))
            if len(isoform2reads[readStruc.name]) < 2: continue
            readsBedList.append("\t".join(readStruc.record[:12]))

    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAnnotation(annoBedRes, novelBedRes, offset=0)
    find_single_as(asType="IR", gene2ReadsDict=annoDict, anno=True, outFile="LR/IR.bed6+", outAppend=False,
                   isoform2reads=isoform2reads)
    find_single_as(asType="IR", gene2ReadsDict=novelDict, anno=False, outFile="LR/IR.bed6+", outAppend=True,
                   isoform2reads=isoform2reads)
    find_single_as(asType="SE", gene2ReadsDict=annoDict, anno=True, outFile="LR/SE.bed12+", outAppend=False,
                   isoform2reads=isoform2reads)
    find_single_as(asType="SE", gene2ReadsDict=novelDict, anno=False, outFile="LR/SE.bed12+", outAppend=True,
                   isoform2reads=isoform2reads)
    find_single_as(asType="A5SS", gene2ReadsDict=annoDict, anno=True, outFile="LR/A5SS.bed6+", outAppend=False,
                   isoform2reads=isoform2reads)
    find_single_as(asType="A5SS", gene2ReadsDict=novelDict, anno=False, outFile="LR/A5SS.bed6+", outAppend=True,
                   isoform2reads=isoform2reads)
    find_single_as(asType="A3SS", gene2ReadsDict=annoDict, anno=True, outFile="LR/A3SS.bed6+", outAppend=False,
                   isoform2reads=isoform2reads)
    find_single_as(asType="A3SS", gene2ReadsDict=novelDict, anno=False, outFile="LR/A3SS.bed6+", outAppend=True,
                   isoform2reads=isoform2reads)

    if os.path.islink("LR/IR.confident.bed6+") or os.path.exists("LR/IR.confident.bed6+"):
        os.remove("LR/IR.confident.bed6+")
    if os.path.islink("LR/SE.confident.bed12+") or os.path.exists("LR/SE.confident.bed12+"):
        os.remove("LR/SE.confident.bed12+")
    if os.path.islink("LR/A5SS.confident.bed12+") or os.path.exists("LR/A5SS.confident.bed6+"):
        os.remove("LR/A5SS.confident.bed6+")
    if os.path.islink("LR/A3SS.confident.bed6+") or os.path.exists("LR/A3SS.confident.bed6+"):
        os.remove("LR/A3SS.confident.bed6+")

    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        if not os.path.exists("NGS"):
            os.makedirs("NGS")
        try:
            # from scanESbyNGS import scanEsByNGS
            from scanSeByNGS import scanSeByNGS, assignGeneNameToSE
            scanSeByNGS(fref=refParams.ref_gpe, has_bin=False, fjunc=dataObj.ngs_junctions,
                        flr=isoformBed, outFile="SEsupportByNGS.bed12+")
            assignGeneNameToSE(fref=refParams.ref_gpe, has_bin=False, fes="SEsupportByNGS.bed12+",
                               outFile="NGS/SE.bed12+", errorOutFile="NGS/novelSE2gene.log")
        except Exception as e:
            print e
        filterFile(originFile="NGS/SE.bed12+", targetFile="LR/SE.bed12+", originField=4, targetField=4,
                   outFile="LR/SE.NGS.bed12+")
        filterFile(originFile="LR/SE.bed12+", targetFile="NGS/SE.bed12+", originField=4, targetField=4,
                   outFile="NGS/SE.LR.bed12+")
        cmd = "{}/juncAssign.pl -g {} {} >junction.assigned.bed12+".format(utilDir, refParams.ref_gpe, dataObj.ngs_junctions)
        subprocess.call(cmd, shell=True)
        cmd = '''
                    {}/AnSSconfirmByJunc.pl -t 5 -j junction.assigned.bed12+ LR/A5SS.bed6+ >NGS/A5SS.LR.bed6+ 2>LR/A5SS.NGS.bed6+
                    {}/AnSSconfirmByJunc.pl -t 3 -j junction.assigned.bed12+ LR/A3SS.bed6+ >NGS/A3SS.LR.bed6+ 2>LR/A3SS.NGS.bed6+
            '''.format(utilDir, utilDir)
        subprocess.call(cmd, shell=True)

        # cmd = '''
        #             awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.LR.bed6+ | hist.R -p=NGS/A5SS.LR.size.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.LR.bed6+ | hist.R -p=NGS/A5SS.LR.InclusionRatio.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $3-$2}' NGS/A3SS.LR.bed6+ | hist.R -p=NGS/A3SS.LR.size.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $5/1000}' NGS/A3SS.LR.bed6+ | hist.R -p=NGS/A3SS.LR.InclusionRatio.pdf 2>/dev/null
        #     '''
        # subprocess.call(cmd, shell=True)
        cmd = '''
                    awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | %s/getNgsA5SS.pl -e 0 >NGS/A5SS.known.bed6+ 2>NGS/A5SS.novel.bed6+
                    awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | %s/getNgsA3SS.pl -e 0 >NGS/A3SS.known.bed6+ 2>NGS/A3SS.novel.bed6+
            ''' % (utilDir, utilDir)

        subprocess.call(cmd, shell=True)
        # cmd = '''
        #             awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.size.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.InclusionRatio.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $3-$2}' NGS/A3SS.{known,novel}.bed6+ | hist.R -p=NGS/A3SS.size.pdf 2>/dev/null
        #             awk '$9>1 && $11>1{print $5/1000}' NGS/A3SS.{known,novel}.bed6+ | hist.R -p=NGS/A3SS.InclusionRatio.pdf 2>/dev/null
        #     '''
        # subprocess.call(cmd, shell=True, executable="/bin/bash")
        filterIrByJunc("LR/IR.bed6+", dataObj.ngs_junctions, "LR/IR.NGS.bed6+")
        makeLink("IR.NGS.bed6+", "LR/IR.confident.bed6+")
        makeLink("SE.NGS.bed12+", "LR/SE.confident.bed12+")
        makeLink("A5SS.NGS.bed6+", "LR/A5SS.confident.bed6+")
        makeLink("A3SS.NGS.bed6+", "LR/A3SS.confident.bed6+")
    else:
        makeLink("IR.bed+", "LR/IR.confident.bed12+")
        makeLink("SE.bed12+", "LR/SE.confident.bed12+")
        makeLink("A5SS.bed6+", "LR/A5SS.confident.bed6+")
        makeLink("A3SS.bed6+", "LR/A3SS.confident.bed6+")

    cmd = '''(cut -f 8,10 --output-delimiter=',' LR/A3SS.confident.bed6+ LR/A5SS.confident.bed6+ LR/IR.confident.bed6+;
                  cut -f 16,18 --output-delimiter=',' LR/SE.confident.bed12+) | tr ',' '\n' | sort -u |
                  {}/filter.pl -o - {} -2 4 -m i > isoformGrouped.AS.confident.bed12+'''.format(utilDir, isoformBed)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    # cmd = '''
    #             awk '{print $3-$2}' LR/A5SS.confident.bed6+ | hist.R -p=LR/A5SS.size.pdf 2>/dev/null
    #             awk '$9>1 && $11>1{print $5/1000}' LR/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=LR/A5SS.InclusionRatio.pdf 2>/dev/null
    #             awk '{print $3-$2}' LR/A3SS.confident.bed6+ | hist.R -p=LR/A3SS.size.pdf 2>/dev/null
    #             awk '$9>1 && $11>1{print $5/1000}' LR/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=LR/A3SS.InclusionRatio.pdf 2>/dev/null
    #         '''
    # subprocess.call(cmd, shell=True, executable="/bin/bash")

    os.chdir(prevDir)
    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {} done!".format(projectName, sampleName)

