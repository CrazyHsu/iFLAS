#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: refine.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from commonFuncs import *
from commonObjs import *
import pybedtools, copy

def isCanonicalSite(strand, dinucleotideType, dinucleotide):
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

def identifyNovelIsoformsByJunctions(gpeFile, bedFile, anno="annoIsoform.bed", novel="novelIsoform.bed"):
    annoOut = open(anno, "w")
    novelOut = open(novel, "w")
    gpeObj = GenePredObj(gpeFile, bincolumn=False)
    bedObj = BedFile(bedFile, type="bed12+")
    knownJuncDict = {}
    for gene in gpeObj.geneName2gpeObj:
        for trans in gpeObj.geneName2gpeObj[gene]:
            chrom = trans.chrom
            if len(trans.introns) == 0:
                if "monoExon" not in knownJuncDict:
                    knownJuncDict["monoExon"] = {chrom: [(trans.txStart, trans.txEnd)]}
                elif chrom not in knownJuncDict["monoExon"]:
                    knownJuncDict["monoExon"][chrom] = [(trans.txStart, trans.txEnd)]
                else:
                    knownJuncDict["monoExon"][chrom].append((trans.txStart, trans.txEnd))
            else:
                transIntronsChain = chrom + ":" + ";".join(map(lambda x: str(trans.introns[x][0])+"-"+str(trans.introns[x][1]), range(len(trans.introns))))
                if transIntronsChain not in knownJuncDict:
                    knownJuncDict[transIntronsChain] = ''
    for i in bedObj.reads:
        readsIntrons = bedObj.reads[i].introns
        chrom = bedObj.reads[i].chrom
        if len(readsIntrons) == 0:
            readStart = bedObj.reads[i].chromStart
            readEnd = bedObj.reads[i].chromEnd
            if chrom not in knownJuncDict["monoExon"]:
                print >> novelOut, bedObj.reads[i]
            else:
                for txStart, txEnd in knownJuncDict["monoExon"][chrom]:
                    if isOverlap((txStart, txEnd), (readStart, readEnd)):
                        overlapLen = getBlockLength(getOverlapOfTuple([(10, 50)], [(20, 60)]))
                        if float(overlapLen)/getBlockLength([(txStart, txEnd)]) >= 0.5:
                            print >> annoOut, bedObj.reads[i]
                            break
                else:
                    print >> novelOut, bedObj.reads[i]
        else:
            readsIntronsChain = chrom + ":" + ";".join(map(lambda x: str(readsIntrons[x][0])+"-"+str(readsIntrons[x][1]), range(len(readsIntrons))))
            if readsIntronsChain in knownJuncDict:
                print >> annoOut, bedObj.reads[i]
            else:
                print >> novelOut, bedObj.reads[i]
    annoOut.close()
    novelOut.close()


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
            singleExonStrandDict[infoList[-3]] = [bedObj.reads[infoList[-3]].strand]

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

    # print
    multiExonAmbiStrandBedList = []
    for i in bedObj.reads:
        if i in singleExonStrandDict and len(set(singleExonStrandDict[i])) == 1:
            if bedObj.reads[i].strand == singleExonStrandDict[i][0]:
                print >>strandConfirmedOut, bedObj.reads[i]
            bedObj.reads[i].strand = singleExonStrandDict[i][0]
            print >>strandAdjustOut, bedObj.reads[i]
        elif i in multiExonReadStrandDict:
            if multiExonReadStrandDict[i]["fwdCanonicalCounter"] - multiExonReadStrandDict[i]["revCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "+":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "+"
                print >> strandAdjustOut, bedObj.reads[i]
            elif multiExonReadStrandDict[i]["revCanonicalCounter"] - multiExonReadStrandDict[i]["fwdCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "-":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "-"
                print >> strandAdjustOut, bedObj.reads[i]
            else:
                chrom, start, end, strand = bedObj.reads[i].chrom, bedObj.reads[i].chromStart, bedObj.reads[i].chromEnd, bedObj.reads[i].strand
                multiExonAmbiStrandBedList.append(" ".join(map(str, [chrom, start, end, i, ".", strand])))
        else:
            print >>strandAdjustOut, bedObj.reads[i]
            print >>strandConfirmedOut, bedObj.reads[i]
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
        if infoList[3] in multiExonAmbiStrandDict and set(multiExonAmbiStrandDict[infoList[3]]) == 1:
            if bedObj.reads[infoList[3]].strand == multiExonAmbiStrandDict[infoList[3]][0]:
                print >>strandConfirmedOut, bedObj.reads[infoList[3]]
            bedObj.reads[infoList[3]].strand = multiExonAmbiStrandDict[infoList[3]][0]
            print >>strandAdjustOut, bedObj.reads[infoList[3]]
        else:
            print >>strandAdjustOut, bedObj.reads[infoList[3]]
            print >>strandConfirmedOut, bedObj.reads[infoList[3]]
    strandAdjustOut.close()
    strandConfirmedOut.close()
    removeFiles(os.getcwd(), [tmpBed])


def readsAssign(transBedFile, readsBedFile, offset=10, minConsenesusIntronN=1, minCoverageOnRead=0.9, singleLine=True, transColNum=13, readsColNum=13, outPrefix="readsAssign", group=True):
    transBedObj = pybedtools.BedTool(transBedFile)
    readsBedObj = pybedtools.BedTool(readsBedFile)
    intersectRes = readsBedObj.intersect(transBedObj, wa=True, wb=True)
    notintersectRes = readsBedObj.intersect(transBedObj, v=True)
    matchedDict = {}
    outList = []
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        if readsColNum > 12:
            readsBed = Bed12Plus("\t".join(infoList[:readsColNum]))
        else:
            readsBed = Bed12("\t".join(infoList[:readsColNum]))
        if transColNum < 13:
            raise Exception("The colmun count of your transcript bed should large than 12")
        transBed = Bed12Plus("\t".join(infoList[readsColNum:readsColNum+transColNum]))
        readStart, readEnd, readsExons, transExons = readsBed.chromStart, readsBed.chromEnd, readsBed.exons, transBed.exons
        overlapLength = getBlockLength(getOverlapOfTuple(readsExons, transExons))
        coverageOnRead = overlapLength/float(getBlockLength(readsExons))
        coverageOnTrans = overlapLength/float(getBlockLength(transExons))
        transGene = transBed.record[12]
        if len(readsBed.exons) == 1:
            if transBed.exonStarts[0] < readEnd and readEnd <= transBed.exonEnds[0] or \
                    transBed.exonStarts[-1] <= readStart and readStart < transBed.exonEnds[-1]:
                inEndExonAndExtension = 1
            else:
                inEndExonAndExtension = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "matchedTransExonCount": [len(transBed.exons)],
                                              "inEndExonAndExtension": [inEndExonAndExtension], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedTransExonCount"].append(len(transBed.exons))
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["inEndExonAndExtension"].append(inEndExonAndExtension)
        else:
            consensusIntronN = getConsensusIntronN(readsBed.exonStarts, readsBed.exonEnds, transBed.exonStarts, transBed.exonEnds, offset=offset)
            readsJuncChain, transJuncChain = readsBed.juncChain, transBed.juncChain
            if re.search(readsJuncChain + "$", transJuncChain) or re.search("^" + readsJuncChain, transJuncChain):
                juncChainFlank = 1
            else:
                juncChainFlank = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "consensusIntronN": [consensusIntronN],
                                              "juncChainFlank": [juncChainFlank], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["consensusIntronN"].append(consensusIntronN)
                matchedDict[readsBed.name]["juncChainFlank"].append(juncChainFlank)

    for readName in matchedDict:
        if singleLine:
            if "juncChainFlank" in matchedDict[readName]:
                readType = "I"
                readIntronN = len(matchedDict[readName]["readsBed"].introns)
                realMatchedTrans, realConsensusIntronN, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                        (readIntronN < minConsenesusIntronN and readIntronN == matchedDict[readName]["consensusIntronN"][j]) or \
                        (readIntronN >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN) or \
                        (matchedDict[readName]["matchedCovOnReads"][j] >= minCoverageOnRead):
                        readType = "E"
                        realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                        realConsensusIntronN.append(matchedDict[readName]["consensusIntronN"][j])
                        realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                        realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                        realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(realMatchedTrans), ",".join(realMatchedGene),
                                    ",".join(map(str, realConsensusIntronN)), ",".join(map(str, realMatchedCovOnReads)),
                                    ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]),
                                    ",".join(map(str, matchedDict[readName]["consensusIntronN"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
            else:
                readType = "I"
                realMatchedTrans, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        continue
                    readType = "E"
                    realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                    realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                    realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                    realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(map(str, realMatchedTrans)), ",".join(map(str, realMatchedGene)), "NA",
                                    ",".join(map(str, realMatchedCovOnReads)), ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]), "NA",
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
        else:
            if "juncChainFlank" in matchedDict[readName]:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    readType = "I"
                    exonNum = len(matchedDict[readName]["readsBed"].exons)
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                            (exonNum < minConsenesusIntronN and exonNum - 1 == matchedDict[readName]["consensusIntronN"][j]) or \
                            (exonNum >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN):
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j], matchedDict[readName]["matchedGenes"][j],
                                    matchedDict[readName]["consensusIntronN"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])
            else:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        readType = "I"
                    else:
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j], "NA",
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])

    for i in notintersectRes:
        outList.append([str(i).strip("\n"), "IG", "NA", "NA", "NA", "NA", "NA"])

    assignOut = open(outPrefix + ".bed12+", "w")
    sortedOutList = sorted(outList, key=lambda x: (x[0].split("\t")[0], int(x[0].split("\t")[1])))
    for item in sortedOutList:
        print >> assignOut, "\t".join(map(str, item))
    assignOut.close()

    if group:
        unambiOut = open(outPrefix + ".unambi.bed12+", "w")
        novelList = []
        for item in sortedOutList:
            if item[3] == "NA" or item[1] != "E":
                novelList.append(item[0].split("\t"))
            else:
                uniqGenes = list(set(item[3].split(",")))
                if len(uniqGenes) == 1:
                    print >> unambiOut, item[0] + "\t" + uniqGenes[0]
                else:
                    novelList.append(item[0].split("\t"))
        sortedNovelList = sorted(novelList, key=lambda x: int(x[1]))
        inc, clusterEnd = 0, 0
        for novelItem in sortedNovelList:
            myStart, myEnd = int(novelItem[1]), int(novelItem[2])
            if myStart < clusterEnd:
                if myEnd > clusterEnd:
                    clusterEnd = myEnd
            else:
                inc += 1
                clusterEnd = myEnd
            print >> unambiOut, "\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc]))
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

    myDfGp = myDf.groupby("gene")

    juncCheck = pd.Series()
    for group in myDfGp.groups:
        tmpDf = myDfGp.get_group(group)
        maxMiddle = max([tmpDf.coverage.median(), tmpDf.coverage.mean()])
        q_abnormal_L = tmpDf.coverage < maxMiddle * 0.05
        juncCheck = pd.concat([juncCheck, q_abnormal_L])

    myDf = myDf[~juncCheck]
    myDf.to_csv("filtered.junction.bed", index=False, header=False, sep="\t", columns=columns[:-1])
    return "filtered.junction.bed"


def refineWithJunc(bedFile=None, junctionFile=None, outFile=None, genomeFasta=None, refGpe=None):
    filteredJunctionFile = filterJunctionFile(junctionFile, bedFile)
    bedJunc = BedFile(bedFile, type="bed12+")
    junction = BedFile(filteredJunctionFile, type="bed12")

    bedJuncDict = bedJunc.getAllJuncDict(onebase=False, withStrand=True)
    bedJunc2trans = bedJunc.getJunc2trans(onebase=False, withStrand=True)
    juncDict = junction.getAllJuncDict(onebase=False, withStrand=True)

    allJuncDict = dict(bedJuncDict, **juncDict)
    allJunc2motif = getSpliceMotifFromJuncDict(allJuncDict, genomeFasta, withStrand=True)

    bedJuncObj = pybedtools.BedTool("\n".join(bedJuncDict.values()), from_string=True)
    juncObj = pybedtools.BedTool("\n".join(juncDict.values()), from_string=True)

    intersectRes = bedJuncObj.intersect(juncObj, wa=True, wb=True, s=True)
    selfIntersectRes = bedJuncObj.intersect(bedJuncObj, wa=True, wb=True, s=True)

    bedJuncIsCovered = {}
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        if infoList[3] not in bedJuncIsCovered:
            bedJuncIsCovered[infoList[3]] = [infoList[9]]
        else:
            bedJuncIsCovered[infoList[3]].append(infoList[9])

    bedJuncSelfIsCovered = {}
    for i in selfIntersectRes:
        infoList = str(i).strip("\n").split("\t")
        if infoList[3] not in bedJuncSelfIsCovered:
            bedJuncSelfIsCovered[infoList[3]] = [infoList[9]]
        else:
            bedJuncSelfIsCovered[infoList[3]].append(infoList[9])

    gpeObj = GenePredObj(refGpe, bincolumn=False)
    knownJuncDict = {}
    for gene in gpeObj.geneName2gpeObj:
        for trans in gpeObj.geneName2gpeObj[gene]:
            chrom = trans.chrom
            if len(trans.introns) == 0:
                if "monoExon" not in knownJuncDict:
                    knownJuncDict["monoExon"] = {chrom: [(trans.txStart, trans.txEnd)]}
                elif chrom not in knownJuncDict["monoExon"]:
                    knownJuncDict["monoExon"][chrom] = [(trans.txStart, trans.txEnd)]
                else:
                    knownJuncDict["monoExon"][chrom].append((trans.txStart, trans.txEnd))
            else:
                transIntronsChain = chrom + ":" + ";".join(
                    map(lambda x: str(trans.introns[x][0]) + "-" + str(trans.introns[x][1]), range(len(trans.introns))))
                if transIntronsChain not in knownJuncDict:
                    knownJuncDict[transIntronsChain] = ''

    unvalidReads = []
    for i in bedJuncDict:
        if i not in juncDict and len(set(bedJuncSelfIsCovered[i])) == 1 and i in bedJuncIsCovered:
            for read in bedJunc2trans[i]:
                readObj = bedJunc.reads[read]
                originReadObj = copy.copy(readObj)
                juncPattern = bedJuncDict[i].split("\t")[3].split(":")[1]

                tmpbedJuncIsCovered = sorted(bedJuncIsCovered[i], key=lambda x: re.split(":|-", x)[1])
                juncSub = ";".join([x.split(":")[1] for x in tmpbedJuncIsCovered])

                canonicalFlag = True
                for j in tmpbedJuncIsCovered:
                    if allJunc2motif[j] not in ["GT-AG", "GC-AG", "AT-AC"]:
                        canonicalFlag = False

                if not validJunction(readObj.chromStart, readObj.chromEnd, re.sub(juncPattern, juncSub, readObj.juncChain)):
                    unvalidReads.append(read)
                    continue
                readObj.juncChain = re.sub(juncPattern, juncSub, readObj.juncChain)
                readObj.getBedFromJuncChain(otherInfo="\t".join(readObj.otherList))
                if originReadObj.get_all_exons_len() * 0.95 > sum(readObj.blockSizes) or \
                        readObj.get_all_exons_len() * 1.05 < sum(readObj.blockSizes) or not canonicalFlag:
                    bedJunc.reads[read] = originReadObj
                else:
                    bedJunc.reads[read] = readObj

    out = open(outFile, "w")
    for read in bedJunc.reads:
        if read not in unvalidReads:
            print >>out, bedJunc.reads[read]
    out.close()


def refineJunc(dataObj=None, refParams=None, dirSpec=None, refine=True, adjust=True):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Refine the collapsed isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    refineDir = os.path.join(baseDir, "refine")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(refineDir)
    processedFa = os.path.join(baseDir, "mapping", "flnc.processed.fa")
    processedBed = os.path.join(baseDir, "mapping", "flnc.addCVandID.bed12+")

    collapsedGff = os.path.join(baseDir, "collapse", "tofu.collapsed.good.gff")
    collapsedGroup = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    cmd = "gtfToGenePred {} tofu.collapsed.gpe -genePredExt".format(collapsedGff)
    subprocess.call(cmd, shell=True)
    cmd = "{}/gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+".format(utilDir)
    subprocess.call(cmd, shell=True)

    if refine:
        if adjust:
            strandAdjust(refParams.ref_genome, refParams.ref_gpe, "tofu.collapsed.bed12+", 0.8, 2,
                         strandAdjust="tofu.strandAdjusted.bed12+", strandConfirmed="tofu.strandConfirm.bed12+")
        else:
            makeLink("tofu.collapsed.bed12+", "tofu.strandConfirm.bed12+")

        if dataObj.ngs_junctions == None:
            dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")
        refineWithJunc(bedFile="tofu.strandConfirm.bed12+", junctionFile=dataObj.ngs_junctions,
                       outFile="tofu.juncAdjusted.bed12+", genomeFasta=refParams.ref_genome, refGpe=refParams.ref_gpe)
        # if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        #     if dataObj.ngs_junctions == None:
        #         dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")
        #     juncScoringParams = "-r {} tofu.strandConfirm.bed12+ -j {}".format(refParams.ref_gpe, dataObj.ngs_junctions)
        # else:
        #     juncScoringParams = "-r {} tofu.strandConfirm.bed12+".format(refParams.ref_gpe)
        # cmd = "{}/juncConsensus.pl -s <({}/juncScoring.pl {}) -l 5 tofu.strandConfirm.bed12+ > tofu.juncAdjusted.bed12+".format(utilDir, utilDir, juncScoringParams)
        # subprocess.call(cmd, shell=True, executable="/bin/bash")

        readsAssign(refParams.ref_bed, "tofu.juncAdjusted.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
    else:
        readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
    cmd = "cut -f1-12,14 tofu.collapsed.assigned.unambi.bed12+ > isoformGrouped.bed12+"
    subprocess.call(cmd, shell=True)
    identifyNovelIsoformsByJunctions(refParams.ref_gpe, "isoformGrouped.bed12+", anno="isoformGrouped.anno.bed12+",
                                     novel="isoformGrouped.novel.bed12+")

    ######## for reads
    cmd = '''seqkit grep {} -f <(cut -f 2 {} | tr ',' '\n') -w 0 > flnc.processed.ignore_id_removed.fa'''.format(processedFa, collapsedGroup)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''{}/filter.pl -o <(cut -f 2 {} | tr ',' '\n') {} -2 4 -m i > flnc.processed.ignore_id_removed.bed12+'''.format(utilDir, collapsedGroup, processedBed)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    readsAssign(refParams.ref_bed, "flnc.processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned", group=True)

    cmd = "cut -f 1-12,15 reads.assigned.unambi.bed12+ | {}/bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS".format(utilDir)
    subprocess.call(cmd, shell=True)
    print getCurrentTime() + " Refine the collapsed isoforms for project {} sample {} done!".format(projectName, sampleName)
