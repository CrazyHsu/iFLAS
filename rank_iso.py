#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: rank_iso.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-20
Last modified: 2022-01-20
'''

# 2. 样本内高质量异构体
# 3. CPC验证编码能力
# 4. salmon查看isoform表达量
# 1. 样本特异异构体
import pandas as pd

from commonFuncs import *
from commonObjs import *


def procAsFile(asFile, asType="IR"):
    isoDict = {}
    if asType == "SE":
        incIndex = 15
        excIndex = 17
    else:
        incIndex = 7
        excIndex = 9

    with open(asFile) as f:
        for i in f.readlines():
            infoList = i.strip().split("\t")
            incIsos = infoList[incIndex]
            excIsos = infoList[excIndex]
            for j in incIsos.split(","):
                if j not in isoDict:
                    isoDict[j] = {}
            for j in excIsos.split(","):
                if j not in isoDict:
                    isoDict[j] = {}
    return isoDict


def getIsoTPM(quant_sf):
    iso2tpmDict = {}
    with open(quant_sf) as f:
        for line in f.readlines()[1:]:
            iso = line.strip("\n").split("\t")[0]
            tpm = line.strip("\n").split("\t")[3]
            iso2tpmDict.update({iso: float(tpm)})
    return iso2tpmDict


def cpc2eval(longestIso):
    from CPC2 import calculate_potential
    calculate_potential(longestIso, "+", 0, "cpc2out")
    codingPotential = {}
    with open("cpc2out.txt") as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip().split("\t")
            isoId, coding_potential = infoList[0], infoList[-1]
            if isoId not in codingPotential:
                codingPotential[isoId] = coding_potential
    return codingPotential


def salmonQuant(dataObj, dirSpec):
    cmd = "salmon index -t longestIso.fa -i longestIso_index --keepDuplicates -p {} 1>/dev/null 2>&1".format(dataObj.single_run_threads)
    subprocess.call(cmd, shell=True)

    from preprocess import renameNGSdata2fastp, processRnaseq
    processRnaseq(dataObj=dataObj, threads=dataObj.single_run_threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
    renameNGSdata2fastp(dataObj=dataObj)
    quant_sf_list = []
    iso2tpm = {}
    if dataObj.ngs_reads_paired == "paired":
        leftReadsRepeats = dataObj.ngs_left_reads.split(";")
        rightReadsRepeats = dataObj.ngs_right_reads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = " ".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = " ".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            salmonOut = "repeat{}.salmon_quant".format(i)
            cmd = "salmon quant -l A -i longestIso_index -1 {} -2 {} -p {} -o {} --consensusSlack 0.5 --preMergeChainSubThresh 0.9 1>/dev/null 2>&1"
            cmd = cmd.format(leftReads, rightReads, dataObj.single_run_threads, salmonOut)
            subprocess.call(cmd, shell=True)
            quant_sf_list.append(os.path.abspath("{}/quant.sf".format(salmonOut)))
            iso2tpm.update({salmonOut: getIsoTPM("{}/quant.sf".format(salmonOut))})
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            singleReadsRepeats = dataObj.ngs_left_reads.split(";")
        elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
            singleReadsRepeats = dataObj.ngs_right_reads.split(";")
        else:
            raise Exception(
                "The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            salmonOut = "repeat{}.salmon_quant".format(i)
            cmd = "salmon quant -l A -i longestIso_index -r {} -p {} -o {} --consensusSlack 0.5 --preMergeChainSubThresh 0.9 1>/dev/null 2>&1"
            cmd = cmd.format(singleReads, dataObj.single_run_threads, salmonOut)
            subprocess.call(cmd, shell=True)
            quant_sf_list.append(os.path.abspath("{}/quant.sf".format(salmonOut)))
            iso2tpm.update({salmonOut: getIsoTPM("{}/quant.sf".format(salmonOut))})
    return iso2tpm

def enumAsIsos(isoDict, isoformBed, collapsedTrans2reads, annoIsoformFile, dataObj, refParams, dirSpec):
    isoform2reads = getDictFromFile(collapsedTrans2reads, sep="\t", inlineSep=",", valueCol=2)

    annoJuncDict = {}
    with open(annoIsoformFile) as f:
        for line in f.readlines():
            juncChain = Bed12Plus(line.strip()).juncChain
            if juncChain not in annoJuncDict:
                annoJuncDict[juncChain] = ""

    juncDict = {}
    gene2isoDict = {}
    for iso in isoDict:
        if isoformBed[iso].juncChain not in juncDict:
            isoLength = getBlockLength(isoformBed[iso].exons)
            juncDict[isoformBed[iso].juncChain] = {"iso": [iso], "gene": [isoformBed[iso].otherList[0]], "longest": [iso, isoLength]}
        else:
            juncDict[isoformBed[iso].juncChain]["iso"].append(iso)
            juncDict[isoformBed[iso].juncChain]["gene"].append(isoformBed[iso].otherList[0])
            if getBlockLength(isoformBed[iso].exons) > juncDict[isoformBed[iso].juncChain]["longest"][1]:
                juncDict[isoformBed[iso].juncChain]["longest"][0] = iso
                juncDict[isoformBed[iso].juncChain]["longest"][1] = getBlockLength(isoformBed[iso].exons)

        gene = isoformBed[iso].otherList[0]
        if gene not in gene2isoDict:
            gene2isoDict[gene] = {"isos": [iso], "count": len(isoform2reads[iso])}
        else:
            gene2isoDict[gene]["isos"].append(iso)
            gene2isoDict[gene]["count"] += len(isoform2reads[iso])

    longestIsoOut = open("longestIso.bed", "w")
    for junc in juncDict:
        longestIso = juncDict[junc]["longest"][0]
        print >> longestIsoOut, str(isoformBed[longestIso])
    longestIsoOut.close()

    cmd = '''cut -f 1-12 longestIso.bed | bedtools getfasta -fi {} -bed - -name -split -s | seqkit replace -w 0 -p "(.*?):(.*)" -r '$1' > longestIso.fa'''.format(refParams.ref_genome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    salmonIsoTMP = salmonQuant(dataObj, dirSpec)
    cpc2IsoCoding = cpc2eval("longestIso.fa")

    isoEnumOut = open("isoEnumerate.txt", "w")
    print >>isoEnumOut, "\t".join(["gene", "longestIso", "similarIsos", "chrom", "junction", "annotation", "readSupport",
                                   "allReadsCount", "readsFreq", "minTPM", "maxTPM", "meanTPM", "codingPotential"])
    for junc in juncDict:
        if len(set(juncDict[junc]["gene"])) != 1: continue
        longestIso = juncDict[junc]["longest"][0]
        similarIsos = ",".join(juncDict[junc]["iso"])
        gene = juncDict[junc]["gene"][0]
        juncSupport = sum([len(isoform2reads[x]) for x in juncDict[junc]["iso"]])
        geneSupport = gene2isoDict[gene]["count"]
        freq = float(juncSupport) / geneSupport
        annotation = "annotated" if junc in annoJuncDict else "novel"
        tpmList = [salmonIsoTMP[x][longestIso] for x in salmonIsoTMP]
        minTPM = min(tpmList)
        maxTPM = max(tpmList)
        meanTPM = sum(tpmList) / float(len(tpmList))
        print >> isoEnumOut, "\t".join(map(str, [gene, longestIso, similarIsos, isoformBed[longestIso].chrom, junc,
                                                 annotation, juncSupport, geneSupport, freq, minTPM, maxTPM,
                                                 meanTPM, cpc2IsoCoding[longestIso]]))
    isoEnumOut.close()
    return "isoEnumerate.txt"


def filterIsos(isoEnumerate, isoformBed, args):
    isoEnumerateDF = pd.read_csv(isoEnumerate, sep="\t")
    # filter by coding potential
    if args.coding:
        isoEnumerateDF = isoEnumerateDF.loc[isoEnumerateDF["codingPotential"] == "coding"]
    # filter by salmon tpm
    isoEnumerateDF = isoEnumerateDF.loc[isoEnumerateDF["minTPM"] >= float(args.min_tpm)]
    # filter by reads frequency
    isoEnumerateDF = isoEnumerateDF.loc[isoEnumerateDF["readsFreq"] >= float(args.reads_freq)]
    # filter by reads count
    isoEnumerateDF = isoEnumerateDF.loc[isoEnumerateDF["readSupport"] >= int(args.read_support)]

    # group filter
    isoEnumerateGroup = isoEnumerateDF.groupby("gene")
    # isoSupportDict = {k: f.groupby('annotation')['readSupport'].apply(list).to_dict() for k, f in isoEnumerateGroup}

    hqGeneDF = isoEnumerateGroup.agg({"annotation": "nunique"}).transform(lambda x: x == 2)
    hqGeneDF = hqGeneDF[hqGeneDF["annotation"] == True]
    # isoEnumerateNew = isoEnumerateGroup.apply(lambda x: x[x["gene"].isin(hqGeneDF.index)])
    # isoEnumerateNew.index = isoEnumerateNew.index.droplevel()
    isoInfoDict = {}
    with open(isoEnumerate) as f:
        for line in f.readlines()[1:]:
            infoList = line.strip().split("\t")
            gene, longestIso, annotation, readSupport = infoList[0], infoList[1], infoList[5], int(infoList[6])
            if gene not in hqGeneDF.index: continue
            if gene not in isoInfoDict:
                isoInfoDict[gene] = {annotation: {longestIso: readSupport}}
            elif annotation not in isoInfoDict[gene]:
                isoInfoDict[gene].update({annotation: {longestIso: readSupport}})
            else:
                isoInfoDict[gene][annotation].update({longestIso: readSupport})

    hqNovelIsos = []
    hqAnnoIsos = []
    for gene in isoInfoDict:
        annoIsos = isoInfoDict[gene]["annotated"]
        novelIsos = isoInfoDict[gene]["novel"]
        annoMaxReadCount = max(annoIsos.values())
        annoMaxReadCountIso = max(annoIsos, key=annoIsos.get)
        for iso in novelIsos:
            if novelIsos[iso] >= annoMaxReadCount * 1.5:
                hqNovelIsos.append(iso)
                hqAnnoIsos.append(annoMaxReadCountIso)
            elif annoMaxReadCount < novelIsos[iso] and novelIsos[iso] < 1.5 * annoMaxReadCount:
                if iso not in isoEnumerateDF.longestIso.to_list() or annoMaxReadCountIso not in isoEnumerateDF.longestIso: continue
                if float(isoEnumerateDF[isoEnumerateDF.longestIso==iso].minTPM) > float(isoEnumerateDF[isoEnumerateDF.longestIso==annoMaxReadCountIso].minTPM):
                    hqNovelIsos.append(iso)
                    hqAnnoIsos.append(annoMaxReadCountIso)

    isoEnumerateDF[isoEnumerateDF.longestIso.isin(hqNovelIsos+hqAnnoIsos)].to_csv("hqIso.txt", sep="\t", index=False, header=True)


def rank_iso(dataObj=None, dirSpec=None, refParams=None, args=None, rawDataObjs=None, optionTools=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Start finding high quality isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    prevDir = os.getcwd()
    isoIdentityDir = os.path.join(baseDir, "isoIdentity")
    resolveDir(isoIdentityDir)
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "IR.confident.bed6+")
    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "SE.confident.bed12+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A3SS.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A5SS.confident.bed6+")

    isoDict = {}
    isoDict.update(procAsFile(irFile, asType="IR"))
    isoDict.update(procAsFile(seFile, asType="SE"))
    isoDict.update(procAsFile(a5ssFile, asType="A5SS"))
    isoDict.update(procAsFile(a3ssFile, asType="A3SS"))

    isoformFile = os.path.join(baseDir, "refine", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    annoIsoformFile = os.path.join(baseDir, "refine", "isoformGrouped.anno.bed12+")

    isoformBed = BedFile(isoformFile, type="bed12+").reads
    isoEnumerate = enumAsIsos(isoDict, isoformBed, collapsedTrans2reads, annoIsoformFile, dataObj, refParams, dirSpec)
    # isoEnumerate = "isoEnumerate.txt"
    filterIsos(isoEnumerate, isoformBed, args)

    if optionTools and (optionTools.merge_data_from_same_strain or args.merge):
        from tissue_spec_iso import tissue_spec_iso
        tissue_spec_iso(dataObj, rawDataObjs, dirSpec, isoformBed, isoDict)
    os.chdir(prevDir)
    print getCurrentTime() + " End finding high quality isoforms for project {} sample {}!".format(projectName, sampleName)
