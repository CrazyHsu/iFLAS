#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: find_pa.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import copy, itertools, pybedtools
import numpy as np
import pandas as pd
from commonFuncs import *
from commonObjs import *
from collections import Counter


def getPAC(depth, distance=20, windowSize=2, paSup=5, rpkmPAC=0, allReadCount=1000000, effLen=1):
    N = len(depth)
    # peaks = []
    # counts = []
    currPeaks = np.zeros(len(depth))
    for i in xrange(N):
        if np.sum(depth[max(0, i - windowSize):min(N, i + windowSize + 1)]) >= paSup:
            currPeaks[i] = depth[i] * 2 + np.mean(depth[max(0, i - windowSize):min(N, i + windowSize + 1)])
            # currPeaks[i] = sum(depth[max(0, i - windowSize - 1):min(N, i + windowSize)])
    currPeaksBak = copy.copy(currPeaks)
    pacDict = {}
    count = 0
    while np.max(currPeaks):
        cp = np.argmax(currPeaks)
        overlapFlag = 0
        maxDiff = 1000000
        closestPAC = ""
        for pac in pacDict:
            if isOverlap((cp-distance, cp+distance+1), (min(pacDict[pac]), max(pacDict[pac]))):
                overlapFlag = 1
                diff2pac = min(abs(cp - min(pacDict[pac])), abs(cp - max(pacDict[pac])))
                if diff2pac < maxDiff:
                    maxDiff = diff2pac
                    closestPAC = pac
        if overlapFlag:
            pacDict[closestPAC].append(cp)
        else:
            count += 1
            pacName = "PAC_{}".format(count)
            pacDict[pacName] = [cp, max([0, cp-windowSize]), min([cp+windowSize+1, N])]
        currPeaks[cp] = 0

    if len(pacDict) == 0: return []
    newPacDict = {pacDict.keys()[0]: copy.copy(pacDict[pacDict.keys()[0]])}
    for pac1 in pacDict.keys()[1:]:
        pac1Start = min(pacDict[pac1])
        pac1End = max(pacDict[pac1])
        for pac2 in newPacDict:
            pac2Start = min(newPacDict[pac2])
            pac2End = max(newPacDict[pac2])
            if isOverlap((pac1Start-distance, pac1End+distance+1), (pac2Start, pac2End)):
                newPacDict[pac2].extend(pacDict[pac1])
                break
        else:
            newPacDict[pac1] = copy.copy(pacDict[pac1])

    finalPeak2Count = []
    for pac in newPacDict:
        pacStart = min(newPacDict[pac])
        pacEnd = max(newPacDict[pac])
        peakPosInPac = int(np.argmax(currPeaksBak[pacStart: pacEnd]))
        peak = range(pacStart, pacEnd)[peakPosInPac]
        count = sum(depth[pacStart:pacEnd])
        if rpkmPAC:
            rpkm = (count * 1000000) / float(allReadCount*effLen)
            if rpkm >= rpkmPAC:
                finalPeak2Count.append((peak, count, [x for x in range(pacStart, pacEnd) if depth[x] > 0]))
        else:
            finalPeak2Count.append((peak, count, [x for x in range(pacStart, pacEnd) if depth[x] > 0]))
    if finalPeak2Count:
        finalPeak2Count = sorted(finalPeak2Count, key=lambda pair: pair[1], reverse=True)
    return finalPeak2Count


def paCluster(readsBedList, distance=20, windowSize=3, outHandle=None, allReadCount=1000000, paSup=5, rpkmPAC=0, confidentPaDict=None):
    chroms = [i.chrom for i in readsBedList]
    strands = [i.strand for i in readsBedList]
    refGene = Counter([i.otherList[2] for i in readsBedList]).most_common()[0][0]
    if len(set(strands)) != 1: return
    chrom = list(set(chroms))[0]
    strand = list(set(strands))[0]
    if strand == "+":
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromEnd)
        paRange = [sortedPAs[0].chromEnd, sortedPAs[-1].chromEnd]
        effectiveLength = abs(paRange[0] - paRange[1]) + 1
        depth = np.zeros(paRange[1] - paRange[0] + 1, int)
        depth2reads = {}
        offest = paRange[0]
        for read in sortedPAs:
            depth[read.chromEnd-offest] += 1
            if read.chromEnd-offest not in depth2reads:
                depth2reads[read.chromEnd-offest] = []
            depth2reads[read.chromEnd-offest].append(read.name)
        finalPeak2Count = getPAC(depth, distance=distance, windowSize=windowSize, paSup=paSup, rpkmPAC=rpkmPAC, allReadCount=allReadCount, effLen=effectiveLength)
    else:
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromStart)
        paRange = [sortedPAs[0].chromStart, sortedPAs[-1].chromStart]
        effectiveLength = abs(paRange[0] - paRange[1]) + 1
        depth = np.zeros(paRange[1] - paRange[0] + 1, int)
        depth2reads = {}
        offest = paRange[0]
        for read in sortedPAs:
            depth[read.chromStart-offest] += 1
            if read.chromStart-offest not in depth2reads:
                depth2reads[read.chromStart-offest] = []
            depth2reads[read.chromStart-offest].append(read.name)
        finalPeak2Count = getPAC(depth, distance=distance, windowSize=windowSize, paSup=paSup, rpkmPAC=rpkmPAC, allReadCount=allReadCount, effLen=effectiveLength)
    for peak, count, pac in finalPeak2Count:
        currPaSite = "{}_{}".format(chrom, peak+offest)
        if confidentPaDict:
            if currPaSite in confidentPaDict:
                annotation = "Known"
            else:
                annotation = "Novel"
        else:
            annotation = "Unknown"
        readNames = list(itertools.chain.from_iterable([depth2reads[x] for x in range(min(pac), max(pac)+1) if x in depth2reads]))
        print >> outHandle, "\t".join(map(str, [chrom, offest + min(pac), offest + max(pac), ",".join(readNames),
                                                len(readNames), strand, offest + peak-1, offest + peak, refGene,
                                                annotation]))


def getPaCluster(readsBed=None, tofuGroup=None, filterByCount=0, threads=None, paClusterOut=None, paDist=20, windowSize=5, rpkmPAC=0, confidentPa=None):
    gene2reads = {}
    readsDict = BedFile(readsBed, type="bed12+").reads
    allReadCount = len(readsDict)
    with open(tofuGroup) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            gene = ".".join(infoList[0].split(".")[:2])
            if gene not in gene2reads:
                gene2reads[gene] = {"allReads": []}
            gene2reads[gene]["allReads"].extend(infoList[1].split(","))
    confidentPaDict = {}
    if confidentPa and validateFile(confidentPa):
        with open(confidentPa) as f:
            for i in f.readlines():
                infoList = i.strip("\n").split("\t")
                for x in xrange(max(0, int(infoList[2])-paDist-1), int(infoList[2])+paDist):
                    paSite = "{}_{}".format(infoList[0], x)
                    confidentPaDict.update({paSite: ""})
    outHandle = open(paClusterOut, "w")
    for z in gene2reads:
        if filterByCount:
            if len(gene2reads[z]["allReads"]) < filterByCount:
                continue
        readsList = gene2reads[z]["allReads"]
        readsBedList = [readsDict[r] for r in readsList if r in readsDict]
        if len(readsBedList) == 0: continue
        chroms = [r.chrom for r in readsBedList]
        mostChrom = Counter(chroms).most_common(1)[0][0]
        readsBedList = [r for r in readsBedList if r.chrom == mostChrom]
        paCluster(readsBedList, distance=paDist, windowSize=windowSize, outHandle=outHandle,
                       allReadCount=allReadCount, paSup=filterByCount, rpkmPAC=rpkmPAC, confidentPaDict=confidentPaDict)
        # paCluster(readsBedList, outHandle=outHandle, distance=paDist)
    outHandle.close()


def motifAroundPA(bed6plus=None, up1=100, down1=100, up2=100, down2=100, refFasta=None, chrLenFile=None):
    singleNucleotideMotif = ["A", "T", "C", "G"]
    sixNucleotideMotif = ["AATAAA", "AAATAA", "ATAAAA", "ATTAAA", "ATAAAT", "TAATAA",
                          "ATAAAG", "AAAATA", "CAATAA", "ATAAAC", "AAAAAA", "AAAAAG"]
    chrLenDict = {}
    with open(chrLenFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            chrLenDict[infoList[0]] = int(infoList[1])
    with open(bed6plus) as f:
        singleNucleotideUpBedList = []
        singleNucleotideDownBedList = []
        sixNucleotideUpBedList = []
        sixNucleotideDownBedList = []
        for i in f.readlines():
            bedObj = Bed6Plus(i)
            chrom, start, end, strand = bedObj.chrom, bedObj.chromStart, bedObj.chromEnd, bedObj.strand
            name = bedObj.name
            if strand == "+":
                if start > up1 and end + down1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up1, start, name, ".", strand])))
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down1, name, ".", strand])))
                if start > up2 and end + down2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up2, start + 5, name, ".", strand])))
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down2 + 5, name, ".", strand])))
            else:
                if start > down1 and end + up1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up1, name, ".", strand])))
                    downLeft = start - down1 - 1 if start - down1 - 1 > 0 else 0
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, downLeft, start-1, name, ".", strand])))
                if start > down2 and end + up2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up2 + 5, name, ".", strand])))
                    downLeft = start - down2 - 5 if start - down2 - 5 > 0 else 0
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, downLeft, start, name, ".", strand])))
        singleNucleotideUpBedObj = pybedtools.BedTool("\n".join(singleNucleotideUpBedList), from_string=True)
        singleNucleotideDownBedObj = pybedtools.BedTool("\n".join(singleNucleotideDownBedList), from_string=True)
        sixNucleotideUpBedObj = pybedtools.BedTool("\n".join(sixNucleotideUpBedList), from_string=True)
        sixNucleotideDownBedObj = pybedtools.BedTool("\n".join(sixNucleotideDownBedList), from_string=True)

        singleNucleotideUpBedFastaRes = singleNucleotideUpBedObj.sequence(refFasta, name=True, tab=True, s=True)
        singleNucleotideDownBedFastaRes = singleNucleotideDownBedObj.sequence(refFasta, name=True, tab=True, s=True)
        sixNucleotideUpBedFastaRes = sixNucleotideUpBedObj.sequence(refFasta, name=True, tab=True, s=True)
        sixNucleotideDownBedFastaRes = sixNucleotideDownBedObj.sequence(refFasta, name=True, tab=True, s=True)

        singleNucleotideUpBedFastaDict = {}
        singleNucleotideDownBedFastaDict = {}
        sixNucleotideUpBedFastaDict = {}
        sixNucleotideDownBedFastaDict = {}
        for res in str(open(singleNucleotideUpBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            singleNucleotideUpBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i+1] for i in range(0, len(infoList[1]), 1)]
        for res in str(open(singleNucleotideDownBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            singleNucleotideDownBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 1] for i in range(0, len(infoList[1]), 1)]
        for res in str(open(sixNucleotideUpBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            sixNucleotideUpBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 6] for i in range(0, len(infoList[1]) - 5, 1)]
        for res in str(open(sixNucleotideDownBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            sixNucleotideDownBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 6] for i in range(0, len(infoList[1]) - 5, 1)]
        singleNucleotideUpDf = pd.DataFrame.from_dict(singleNucleotideUpBedFastaDict, orient="index")
        singleNucleotideDownDf = pd.DataFrame.from_dict(singleNucleotideDownBedFastaDict, orient="index")
        sixNucleotideUpDf = pd.DataFrame.from_dict(sixNucleotideUpBedFastaDict, orient="index")
        sixNucleotideDownDf = pd.DataFrame.from_dict(sixNucleotideDownBedFastaDict, orient="index")

        for i in singleNucleotideMotif:
            singleNucleotideOut = open(i+".nucleotide", "w")
            singleNucleotideUpDf1 = singleNucleotideUpDf == i
            upPercent = singleNucleotideUpDf1.sum() / float(len(singleNucleotideUpDf1))
            for j in range(len(upPercent)):
                print >>singleNucleotideOut, "\t".join(map(str, [j - up1, upPercent[j]]))

            singleNucleotideDownDf1 = singleNucleotideDownDf == i
            downPercent = singleNucleotideDownDf1.sum() / float(len(singleNucleotideDownDf1))
            for j in range(len(downPercent)):
                print >>singleNucleotideOut, "\t".join(map(str, [j + 1, downPercent[j]]))
            singleNucleotideOut.close()

        for i in sixNucleotideMotif:
            sixNucleotideOut = open(i + ".PAS", "w")
            sixNucleotideUpDf1 = sixNucleotideUpDf == i
            upPercent = sixNucleotideUpDf1.sum() / float(len(sixNucleotideUpDf1))
            for j in range(len(upPercent)):
                print >>sixNucleotideOut, "\t".join(map(str, [j - up2, upPercent[j]]))

            sixNucleotideDownDf1 = sixNucleotideDownDf == i
            downPercent = sixNucleotideDownDf1.sum() / float(len(sixNucleotideDownDf1))
            for j in range(len(downPercent)):
                print >>sixNucleotideOut, "\t".join(map(str, [j + 1, downPercent[j]]))
            sixNucleotideOut.close()


def find_pa(dataObj=None, refParams=None, dirSpec=None, filterPaByRPKM=0, filterPaByCount=5, confidentPa=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Polyadenylation analysis for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    paDir = os.path.join(baseDir, "as_events", "pa")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(paDir)
    makeLink(os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt"), "tofu.collapsed.group.txt")
    makeLink(os.path.join(baseDir, "refine", "reads.assigned.unambi.bed12+"), "reads.assigned.unambi.bed12+")
    makeLink(os.path.join(baseDir, "refine", "isoformGrouped.bed12+"), "isoformGrouped.bed12+")

    with open("reads.assigned.unambi.bed12+") as f:
        cleavageList = {}
        count = 0
        cleavageOut = open("cleavage.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            chrom, start, end, strand = lineInfo[0], int(lineInfo[1]), int(lineInfo[2]), lineInfo[5]
            if strand == "+":
                keyStr = "{}:{}-{}:{}".format(chrom, end-1, end, strand)
                if keyStr not in cleavageList:
                    count += 1
                    cleavageList.update({keyStr: 1})
                    print >>cleavageOut, "\t".join(map(str, [chrom, end - 1, end, "ClevageSite{}".format(count), 0, strand]))
            elif strand == "-":
                keyStr = "{}:{}-{}:{}".format(chrom, start, start + 1, strand)
                if keyStr not in cleavageList:
                    count += 1
                    cleavageList.update({keyStr: 1})
                    print >>cleavageOut, "\t".join(map(str, [chrom, start, start + 1, "ClevageSite{}".format(count), 0, strand]))
        cleavageOut.close()

    cmd = '''{}/bedClosest.pl cleavage.bed6 | tee adjacentCleavage.tsv | cut -f13 | {}/distrCurve.R -w=10 -xl=2 -d -x='Binned Distance between Adjacent Cleavage Site' -p=adjacentCleavageDistr.pdf 2>/dev/null'''.format(utilDir, utilDir)
    subprocess.call(cmd, shell=True)

    # pa(polish_flnc_cluster="clustered_unclustered.merged_report.csv", bed12="target.transcript.correlated.flnc.sorted.bed12+",
    #    tofu_group="tofu.collapsed.group.txt", threads=refParams.threads, out="paCluster.bed8+")
    getPaCluster(readsBed="reads.assigned.unambi.bed12+", tofuGroup="tofu.collapsed.group.txt",
                 threads=dataObj.single_run_threads, paClusterOut="paCluster.bed8+", paDist=24, windowSize=3,
                 rpkmPAC=filterPaByRPKM, filterByCount=filterPaByCount, confidentPa=confidentPa)
    # cmd = '''awk '$5>1{print $3-$2}' paCluster.bed8+ |
    #          distrCurve.R -d -x='Cluster Size (limited in 1-100)' -y='Density' -m='Distribution of Cluster Size' -x1=0 -x2=100 -b=1 -p=paClusterSize.pdf 2>/dev/null
    # '''
    # subprocess.call(cmd, shell=True)
    paBed6 = open("PA.bed6+", "w")
    with open("paCluster.bed8+") as f:
        for line in f.readlines():
            lineInfo = line.strip("\n").split("\t")
            if int(lineInfo[4]) >= filterPaByCount:
                print >> paBed6, "\t".join(map(str, [lineInfo[0], lineInfo[6], lineInfo[7]] + lineInfo[3:6] + lineInfo[9:]))
    paBed6.close()
    cmd = "{}/paGroup.pl isoformGrouped.bed12+ >isoform.paGrouped.tsv 2>isoform.paGrouped.bed6".format(utilDir)
    subprocess.call(cmd, shell=True)
    cmd = '''{}/3endRevise.pl -p PA.bed6+ <(cut -f 1-12,15 reads.assigned.unambi.bed12+) > reads.3endRevised.bed12+'''.format(utilDir)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''cut -f 4 PA.bed6+ | tr ',' '\n' | {}/filter.pl -o - reads.3endRevised.bed12+ -2 4 -m i | {}/paGroup.pl >reads.paGrouped.tsv 2>reads.paGrouped.bed6'''.format(utilDir, utilDir)
    subprocess.call(cmd, shell=True)

    # cmd = "PAClassbyRead.pl -a reads.assigned.unambi.bed12+ <(cut -f1-8 paCluster.bed8+) >paCluster.type.bed8+ 2>singleExonReadWithExonInMEread.bed12+"
    # subprocess.call(cmd, shell=True, executable="/bin/bash")
    #
    # deSingleExon = open("paCluster.deSingleExonRead.bed8+", "w")
    # singleExon = open("paCluster.singleExonRead.bed8+", "w")
    # with open("paCluster.type.bed8+") as f:
    #     for line in f:
    #         lineInfo = line.strip("\n").split("\t")
    #         if lineInfo[8] != "SE":
    #             print >> deSingleExon, line.strip("\n")
    #         elif lineInfo[8] == "SE":
    #             print >> singleExon, line.strip("\n")
    # deSingleExon.close()
    # singleExon.close()

    resolveDir("motif")
    motifAroundPA("../PA.bed6+", up1=100, down1=100, up2=50, down2=0, refFasta=refParams.ref_genome, chrLenFile=refParams.ref_size)

    cmd = '''%s/lines.R -w=10 -y=Frequency -x='Distance Relative to PA' -m='Distribution of Nucleotide Frequency around PA' -p=nucleotide.pdf *.nucleotide 2>/dev/null''' %(utilDir)
    subprocess.call(cmd, shell=True)
    cmd = '''%s/lines.R -p=PAS.1.pdf {AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=10 2>/dev/null &&
             %s/lines.R -p=PAS.2.pdf {ATAAAG,AAAATA,CAATAA,ATAAAC,AAAAAA,AAAAAG}.PAS -x1=-50 -x2=0 -w=10 2>/dev/null
    ''' % (utilDir, utilDir)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    otherMotifFileList = ["ATAAAC.PAS", "ATAAAG.PAS", "CAATAA.PAS", "AAAAAA.PAS", "AAAAAG.PAS", "AAAATA.PAS"]
    motifPercentDF = [pd.read_csv(i, sep="\t", index_col=0, header=None) for i in otherMotifFileList]
    motifPercentSummedDF = pd.concat(motifPercentDF, axis=1).sum(axis=1)
    motifPercentSummedDF.to_frame().to_csv("Other.PAS", sep="\t", header=None)

    cmd = '''%s/lines.R -p=PAS.pdf {Other,AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=10 2>/dev/null''' % (utilDir)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print getCurrentTime() + " Polyadenylation analysis for project {} entry {} done!".format(projectName, sampleName)
