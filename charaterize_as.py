#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: charaterize_as.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import pybedtools
from collections import Counter
from commonFuncs import *
from commonObjs import *

def drawSSmotif(asMotif=None, outPrefix=None, asType="IR"):
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    with open(asMotif) as f:
        lineList = f.readlines()
        mySum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList])
        tmp = open("{}.tmp.txt".format(outPrefix), "w")
        if len(lineList) < 4:
            for i in lineList:
                lineInfo = i.strip("\n").split("\t")
                if asType == "A3SS":
                    spliceSite = "-{}".format(lineInfo[0])
                elif asType == "A5SS":
                    spliceSite = "{}-".format(lineInfo[0])
                else:
                    spliceSite = lineInfo[0]
                print >> tmp, "\t".join([spliceSite, str(float(lineInfo[1])/mySum)])
        else:
            for i in lineList[0:3]:
                lineInfo = i.strip("\n").split("\t")
                if asType == "A3SS":
                    spliceSite = "-{}".format(lineInfo[0])
                elif asType == "A5SS":
                    spliceSite = "{}-".format(lineInfo[0])
                else:
                    spliceSite = lineInfo[0]
                print >> tmp, "\t".join([spliceSite, str(float(lineInfo[1]) / mySum)])
            otherSum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList[3:]])
            print >> tmp, "\t".join(["Other", str(float(otherSum)/mySum)])
        tmp.close()
        cmd = "cat {}.tmp.txt | {}/bar.R -fillV=V1 -fp -lgPos=top -w=10 -p={}.ssMotif.pdf 2>/dev/null".format(outPrefix, utilDir, outPrefix)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        # os.remove("tmp.txt")
        # os.remove(asMotif)


def getAnnoASList(inFile, outFile, PA=False, append=False, uniq=False):
    with open(inFile) as f:
        if not append:
            out = open(outFile, "w")
        else:
            out = open(outFile, "w+")
        asLst = []
        if not PA:
            for line in f:
                if "Novel" in line: continue
                lineInfo = line.strip("\n").split("\t")
                asEvent = "{}:{}:{}".format(lineInfo[0], lineInfo[3], lineInfo[5])
                if uniq:
                    if asEvent not in asLst:
                        asLst.append(asEvent)
                    else:
                        continue
                print >>out, asEvent
        else:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                paPos = lineInfo[3].split(",")
                if len(paPos) > 1:
                    asEvent = ":".join(lineInfo[0:4])
                    if uniq:
                        if asEvent not in asLst:
                            asLst.append(asEvent)
                        else:
                            continue
                    print >>out, asEvent
        out.close()


def getASstatistics(asType="IR", asFile=None, annoFile=None, novelFile=None, outFile=None):
    out = open(outFile, "w")
    if asType == "PA":
        print >> out, "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
        with open(annoFile) as f1:
            for line in f1:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Known", lineInfo[2], lineInfo[3]])
        with open(novelFile) as f2:
            for line in f2:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Novel", lineInfo[2], lineInfo[3]])
    else:
        asDict = {}
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                key = ":".join([lineInfo[0], lineInfo[3], lineInfo[5]])
                asDict[key] = "\t".join([key, line.strip("\n")])
        annoList = filterFile(originFile=annoFile, targetFile=asDict, returnFlag=True)
        novelList = filterFile(originFile=novelFile, targetFile=asDict, returnFlag=True)
        os.remove("filtered.txt")
        geneCol = 7
        if asType == "IR":
            print >> out, "##AS ID is composed of Gene:Retained intron start-Retained intron end"
        elif asType == "SE":
            geneCol = 13
            print >> out, "##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start"
            print >> out, "##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]"
        elif asType == "A3SS":
            print >> out, "##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end"
        elif asType == "A5SS":
            print >> out, "##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end"
        print >> out, "#Chr\tStrand\tKnown or Novel\tAS ID\tGene"

        for i in annoList:
            tmpList = asDict[i].strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Known", tmpList[4], tmpList[geneCol]])
        for i in novelList:
            tmpList = asDict[i].strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Novel", tmpList[4], tmpList[geneCol]])

    out.close()


def getDist2TTS(refParams=None, paGroup=None):
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    with open(paGroup) as f:
        out = open("lrPA.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if lineInfo[5] == "+":
                print >>out, "\t".join(map(str, [lineInfo[0], int(lineInfo[2])-1, lineInfo[2]] + lineInfo[3:]))
            else:
                print >>out, "\t".join(map(str, [lineInfo[0], lineInfo[1], int(lineInfo[1])+1] + lineInfo[3:]))
        out.close()
    cmd = '''
        bedtools closest -a <(sort -k1,1 -k2,2n lrPA.bed6) -b <({}/gpeFeature.pl --tts {}|
        sort -k1,1 -k2,2n) -s -D a | {}/select.pl -i 13,4 | sort -u | tee lrPA2TTS.tsv | 
        cut -f1 | {}/box.R -w=10 -height=6 -ng -nJ -no -y='Distance to TTS' -p=lrPA2TTS.pdf
    '''.format(utilDir, refParams.ref_gpe, utilDir, utilDir)
    subprocess.call(cmd, shell=True, executable="/bin/bash")


def getSpliceSite(asType=None, asFile=None, outFile=None):
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    if asType == "IR":
        out = open(outFile, "w")
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split(":")
                posList = lineInfo[2].split("-")
                print >>out, "\t".join([lineInfo[0], posList[0], posList[1], lineInfo[3]])
        out.close()
    elif asType == "SE":
        cmd = "{}/seDecompose.pl confident.SE.lst >SE.inc.splicesite 2>SE.exc.splicesite".format(utilDir)
        subprocess.call(cmd, shell=True)
    elif asType == "A3SS":
        cmd = "{}/anssDecompose.pl -n 5 confident.A5SS.lst >A5SS.inc.splicesite 2>A5SS.exc.splicesite".format(utilDir)
        subprocess.call(cmd, shell=True)
    elif asType == "A5SS":
        cmd = "{}/anssDecompose.pl -n 3 confident.A3SS.lst >A3SS.inc.splicesite 2>A3SS.exc.splicesite".format(utilDir)
        subprocess.call(cmd, shell=True)


def splicesite2seq(refFasta, spliceSite, noFiveEnd=False, noThreeEnd=False, outFile=None):
    donorList, acceptorList = [], []
    exonDist, intronDist = 0, 2
    with open(spliceSite) as f:
        count = 0
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            chrom, leftSite, rightSite, strand = infoList[0], infoList[1], infoList[2], infoList[3]
            count += 1
            juncName = "junction_" + str(count)
            if infoList[3] == "+":
                if not noFiveEnd:
                    leftSite = int(leftSite)
                    donorBed = "\t".join(map(str, [chrom, leftSite - exonDist, leftSite + intronDist, juncName, ".", strand]))
                    donorList.append(donorBed)
                if not noThreeEnd:
                    rightSite = int(rightSite)
                    acceptorBed = "\t".join(map(str, [chrom, rightSite - intronDist, rightSite + exonDist, juncName, ".", strand]))
                    acceptorList.append(acceptorBed)
            else:
                if not noFiveEnd:
                    rightSite = int(rightSite)
                    donorBed = "\t".join(map(str, [chrom, rightSite - intronDist, rightSite + exonDist, juncName, ".", strand]))
                    donorList.append(donorBed)
                if not noThreeEnd:
                    leftSite = int(leftSite)
                    acceptorBed = "\t".join(map(str, [chrom, leftSite - exonDist, leftSite + intronDist, juncName, ".", strand]))
                    acceptorList.append(acceptorBed)
    donorBedObj = pybedtools.BedTool("\n".join(donorList), from_string=True)
    acceptorBedObj = pybedtools.BedTool("\n".join(acceptorList), from_string=True)
    junc2seq = {}
    if not noFiveEnd and not noThreeEnd:
        donorBedGetfastaRes = donorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        acceptorBedGetfastaRes = acceptorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(donorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
        for i in str(open(acceptorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = junc2seq[junc] + "-" + baseSeq
    elif not noFiveEnd and noThreeEnd:
        donorBedGetfastaRes = donorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(donorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
    elif noFiveEnd and not noThreeEnd:
        acceptorBedGetfastaRes = acceptorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(acceptorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
    d_tmp = [(v, k) for k, v in Counter(junc2seq.values()).iteritems()]
    d_tmp.sort(reverse=True)
    out = open(outFile, "w")
    for v, k in d_tmp:
        print >> out, "\t".join(map(str, [k, v]))
    out.close()


def splicesite2figure(refParams=None):
    refGenome = refParams.ref_genome
    splicesite2seq(refGenome, "IR.splicesite", outFile="IR.tmp")
    drawSSmotif(asMotif="IR.tmp", outPrefix="IR", asType="IR")
    splicesite2seq(refGenome, "IR.anno.splicesite", outFile="IR.anno.tmp")
    drawSSmotif(asMotif="IR.anno.tmp", outPrefix="IR.anno", asType="IR")
    splicesite2seq(refGenome, "IR.novel.splicesite", outFile="IR.novel.tmp")
    drawSSmotif(asMotif="IR.novel.tmp", outPrefix="IR.novel", asType="IR")
    splicesite2seq(refGenome, "SE.inc.splicesite", outFile="SE.inc.tmp")
    drawSSmotif(asMotif="SE.inc.tmp", outPrefix="SE.inc", asType="SE")
    splicesite2seq(refGenome, "SE.exc.splicesite", outFile="SE.exc.tmp")
    drawSSmotif(asMotif="SE.exc.tmp", outPrefix="SE.exc", asType="SE")
    splicesite2seq(refGenome, "A5SS.inc.splicesite", noThreeEnd=True, outFile="A5SS.inc.tmp")
    drawSSmotif(asMotif="A5SS.inc.tmp", outPrefix="A5SS.inc", asType="A5SS")
    splicesite2seq(refGenome, "A5SS.exc.splicesite", noThreeEnd=True, outFile="A5SS.exc.tmp")
    drawSSmotif(asMotif="A5SS.exc.tmp", outPrefix="A5SS.exc", asType="A5SS")
    splicesite2seq(refGenome, "A3SS.inc.splicesite", noFiveEnd=True, outFile="A3SS.inc.tmp")
    drawSSmotif(asMotif="A3SS.inc.tmp", outPrefix="A3SS.inc", asType="A3SS")
    splicesite2seq(refGenome, "A3SS.exc.splicesite", noFiveEnd=True, outFile="A3SS.exc.tmp")
    drawSSmotif(asMotif="A3SS.exc.tmp", outPrefix="A3SS.exc", asType="A3SS")


def charaterize_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " ASE Characterization for project {} entry {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    characAsDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "as_events", "characterization")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(characAsDir)
    getAnnoASList("../ordinary_as/LR/SE.bed12+", "LR.SE.lst")
    getAnnoASList("../ordinary_as/LR/A5SS.bed6+", "LR.A5SS.lst")
    getAnnoASList("../ordinary_as/LR/A3SS.bed6+", "LR.A3SS.lst")
    getAnnoASList("../ordinary_as/LR/IR.bed6+", "LR.IR.lst")
    getAnnoASList("../pa/reads.paGrouped.tsv", "LR.APA.lst", PA=True)

    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        if dataObj.ngs_junctions == None:
            dataObj.ngs_junctions = os.path.join(dirSpec.out_dir, projectName, sampleName, "mapping", "rna-seq", "reassembly", "junctions.bed")
        getAnnoASList("../ordinary_as/NGS/SE.bed12+", "NGS.SE.lst")
        getAnnoASList("../ordinary_as/NGS/A5SS.known.bed6+", "NGS.A5SS.lst")
        getAnnoASList("../ordinary_as/NGS/A5SS.LR.bed6+", "NGS.A5SS.lst", append=True, uniq=True)
        getAnnoASList("../ordinary_as/NGS/A3SS.known.bed6+", "NGS.A3SS.lst")
        getAnnoASList("../ordinary_as/NGS/A3SS.LR.bed6+", "NGS.A3SS.lst", append=True, uniq=True)

        filterFile(originFile="NGS.SE.lst", targetFile="LR.SE.lst", outFile="Common.SE.lst")
        filterFile(originFile="NGS.A5SS.lst", targetFile="LR.A5SS.lst", outFile="Common.A5SS.lst")
        filterFile(originFile="NGS.A3SS.lst", targetFile="LR.A3SS.lst", outFile="Common.A3SS.lst")

        makeLink("Common.SE.lst", "confident.SE.lst")
        makeLink("Common.A5SS.lst", "confident.A5SS.lst")
        makeLink("Common.A3SS.lst", "confident.A3SS.lst")

        cmd = '''
                %s/venn.R {NGS,LR}.SE.lst   -rD=90 -cP=180,0 -p=SE.pdf 1>/dev/null 2>&1
                %s/venn.R {NGS,LR}.A5SS.lst -rD=90 -cP=180,0 -p=A5SS.pdf 1>/dev/null 2>&1
                %s/venn.R {NGS,LR}.A3SS.lst -rD=90 -cP=180,0 -p=A3SS.pdf 1>/dev/null 2>&1
        ''' % (utilDir, utilDir, utilDir)
        subprocess.call(cmd, shell=True)

        filterFile(originFile="NGS.SE.lst", targetFile="LR.SE.lst", outFile="LR.specific.SE.lst", mode="e")
        filterFile(originFile="NGS.A5SS.lst", targetFile="LR.A5SS.lst", outFile="LR.specific.A5SS.lst", mode="e")
        filterFile(originFile="NGS.A3SS.lst", targetFile="LR.A3SS.lst", outFile="LR.specific.A3SS.lst", mode="e")

    else:
        makeLink("LR.SE.lst", "confident.SE.lst")
        makeLink("LR.A5SS.lst", "confident.A5SS.lst")
        makeLink("LR.A3SS.lst", "confident.A3SS.lst")

    # Classify APE into annotated or novel according to reference gene model
    cmd = '''awk '{{print $0"\t"$12}}' {} | {}/gpe2bed.pl -p >reference.bed12+'''.format(refParams.ref_gpe, utilDir)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    with open("reference.bed12+") as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            readsBedList.append("\t".join(readStruc.record[:12]))
    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    from find_as import getAnnotation, find_single_as
    annoDict, novelDict = getAnnotation(annoBedRes, novelBedRes, offset=0)
    find_single_as(asType="IR", gene2ReadsDict=annoDict, anno=True, outFile="IR.reference.bed6+", outAppend=False)
    find_single_as(asType="IR", gene2ReadsDict=novelDict, anno=False, outFile="IR.reference.bed6+", outAppend=True)
    find_single_as(asType="SE", gene2ReadsDict=annoDict, anno=True, outFile="SE.reference.bed6+", outAppend=False)
    find_single_as(asType="SE", gene2ReadsDict=novelDict, anno=False, outFile="SE.reference.bed6+", outAppend=True)
    find_single_as(asType="A5SS", gene2ReadsDict=annoDict, anno=True, outFile="A5SS.reference.bed6+", outAppend=False)
    find_single_as(asType="A5SS", gene2ReadsDict=novelDict, anno=False, outFile="A5SS.reference.bed6+", outAppend=True)
    find_single_as(asType="A3SS", gene2ReadsDict=annoDict, anno=True, outFile="A3SS.reference.bed6+", outAppend=False)
    find_single_as(asType="A3SS", gene2ReadsDict=novelDict, anno=False, outFile="A3SS.reference.bed6+", outAppend=True)


    getAnnoASList("IR.reference.bed6+", "IR.reference.lst")
    getAnnoASList("SE.reference.bed6+", "SE.reference.lst")
    getAnnoASList("A3SS.reference.bed6+", "A3SS.reference.lst")
    getAnnoASList("A5SS.reference.bed6+", "A5SS.reference.lst")
    cmd = "{}/paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6".format(utilDir)
    subprocess.call(cmd, shell=True)
    cmd = "{}/paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../pa/reads.paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv".format(utilDir)
    subprocess.call(cmd, shell=True)

    filterFile(originFile="IR.reference.lst", targetFile="LR.IR.lst", outFile="IR.anno.lst", mode="i")
    filterFile(originFile="IR.reference.lst", targetFile="LR.IR.lst", outFile="IR.novel.lst", mode="e")
    filterFile(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.anno.lst", mode="i")
    filterFile(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.novel.lst", mode="e")
    filterFile(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.anno.lst", mode="i")
    filterFile(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.novel.lst", mode="e")
    filterFile(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.anno.lst", mode="i")
    filterFile(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.novel.lst", mode="e")

    # Get Statistic of APE
    getASstatistics(asType="IR", asFile="../ordinary_as/LR/IR.bed6+", annoFile="IR.anno.lst", novelFile="IR.novel.lst",
                    outFile="statistic.IR.tsv")
    getASstatistics(asType="PA", annoFile="paGroup.anno.tsv", novelFile="paGroup.novel.tsv",
                    outFile="statistic.APA.tsv")
    getASstatistics(asType="SE", asFile="../ordinary_as/LR/SE.confident.bed12+", annoFile="SE.anno.lst",
                    novelFile="SE.novel.lst", outFile="statistic.SE.tsv")
    getASstatistics(asType="A3SS", asFile="../ordinary_as/LR/A3SS.confident.bed6+", annoFile="A3SS.anno.lst",
                    novelFile="A3SS.novel.lst", outFile="statistic.A3SS.tsv")
    getASstatistics(asType="A5SS", asFile="../ordinary_as/LR/A5SS.confident.bed6+", annoFile="A5SS.anno.lst",
                    novelFile="A5SS.novel.lst", outFile="statistic.A5SS.tsv")

    # Event Decompose and Motif Evaluation
    getSpliceSite(asType="IR", asFile="LR.IR.lst", outFile="IR.splicesite")
    getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
    getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
    getSpliceSite(asType="SE", asFile="confident.SE.lst")
    getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst")
    getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst")
    splicesite2figure(refParams=refParams)

    getDist2TTS(refParams=refParams, paGroup="../pa/reads.paGrouped.bed6")

    os.chdir(prevDir)
    print getCurrentTime() + " ASE Characterization for project {} entry {} done!".format(projectName, sampleName)
