#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: diff_as.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-19
Last modified: 2022-01-19
'''

import copy
import pandas as pd
from itertools import combinations
from commonFuncs import *
from commonObjs import *

def getSample2Bam(samples, compCondFile, dirSpec=None):
    sample2bamDict = {"cond2bam": {}}
    compCondHeader = ["project", "sampleName", "condition", "repeat", "bamFile", "paired", "readsLength"]
    if os.path.isfile(compCondFile):
        compCond = pd.read_csv(compCondFile, sep="\t")
        if len(compCond.columns) == 7 and len(set(compCond.columns) & set(compCondHeader)) == 7:
            for i, row in compCond.iterrows():
                bamFile = row.bamFile
                if row.condition not in sample2bamDict["cond2bam"]:
                    sample2bamDict["cond2bam"][row.condition] = [(row.project + "_" + row.condition, bamFile, row.paired, row.readsLength)]
                else:
                    sample2bamDict["cond2bam"][row.condition].append((row.project + "_" + row.condition, bamFile, row.paired, row.readsLength))
        else:
            for dataObj in samples:
                alignDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "mapping", "rna-seq", "alignment")
                if dataObj.ngs_reads_paired == "paired":
                    leftReadsRepeats = dataObj.ngs_left_reads.split(";")
                    for i in range(len(leftReadsRepeats)):
                        repeatName = "repeat" + str(i)
                        bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                        if dataObj.condition not in sample2bamDict["cond2bam"]:
                            sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                        else:
                            sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
                else:
                    if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
                        singleReadsRepeats = dataObj.ngs_left_reads.split(";")
                    elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
                        singleReadsRepeats = dataObj.ngs_right_reads.split(";")
                    else:
                        raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

                    for i in range(len(singleReadsRepeats)):
                        repeatName = "repeat" + str(i)
                        bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                        if dataObj.condition not in sample2bamDict["cond2bam"]:
                            sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                        else:
                            sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
    else:
        for dataObj in samples:
            alignDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "mapping", "rna-seq", "alignment")
            if dataObj.ngs_reads_paired == "paired":
                leftReadsRepeats = dataObj.ngs_left_reads.split(";")
                for i in range(len(leftReadsRepeats)):
                    repeatName = "repeat" + str(i)
                    bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                    if dataObj.condition not in sample2bamDict["cond2bam"]:
                        sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                    else:
                        sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
            else:
                if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
                    singleReadsRepeats = dataObj.ngs_left_reads.split(";")
                elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
                    singleReadsRepeats = dataObj.ngs_right_reads.split(";")
                else:
                    raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

                for i in range(len(singleReadsRepeats)):
                    repeatName = "repeat" + str(i)
                    bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                    if dataObj.condition not in sample2bamDict["cond2bam"]:
                        sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                    else:
                        sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
    return sample2bamDict


def validateSamples(compCondFile, dataToProcess):
    samples = []
    conditions = set()
    if compCondFile:
        if validateFile(compCondFile):
            compCond = pd.read_csv(compCondFile, sep="\t")
            if "condition" in compCond.columns:
                for i, row in compCond.iterrows():
                    conditions.add(row.condition)
                    samples.extend([x for x in dataToProcess if x.project_name == row.project and x.condition == row.condition])
            else:
                raise Exception("Please set the condition in {}".format(compCondFile))
            samples = list(set(samples))
        else:
            compCondList = compCondFile.split(",")
            samples = [x for x in dataToProcess if x.condition in set(compCondList)]
            conditions = set(compCondList)
    else:
        samples = dataToProcess
    return samples, conditions


def mergeIsoforms(samples=None, dirSpec=None, pu_filter=False, refInfoParams=None):
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    tmpDict = {}
    mergedIso2ReadsBed = open("all_sample_merged_iso.bed", "w")
    for i in samples:
        isoformGroupedBed12 = os.path.join(dirSpec.out_dir, i.project_name, i.sample_name, "refine", "isoformGrouped.bed12+")
        if pu_filter:
            puFilterIso = os.path.join(dirSpec.out_dir, i.project_name, i.sample_name, "hqIsoforms", "hq.collapsed.bed12+")
            if validateFile(puFilterIso):
                cmd = '''{}/filter.pl -o {} {} -1 4 -2 4 -m i > as_isoform.bed12+'''.format(utilDir, puFilterIso, isoformGroupedBed12)
                subprocess.call(cmd, shell=True)
            else:
                from iso_pu import iso_pu1
                iso_pu1(dataObj=i, dirSpec=dirSpec, refParams=refInfoParams[i.ref_strain])
        else:
            aseDir = os.path.join(dirSpec.out_dir, i.project_name, i.sample_name, "as_events", "ordinary_as")
            cmd = '''(cut -f 8,10 --output-delimiter=',' {}/LR/A3SS.confident.bed6+ {}/LR/A5SS.confident.bed6+ {}/LR/IR.confident.bed6+;
                  cut -f 16,18 --output-delimiter=',' {}/LR/SE.confident.bed12+) | tr ',' '\n' | sort -u |
                  {}/filter.pl -o - {} -2 4 -m i > as_isoform.bed12+'''.format(aseDir, aseDir, aseDir, aseDir, utilDir, isoformGroupedBed12)
            subprocess.call(cmd, shell=True)
        isoBedObj = BedFile("as_isoform.bed12+", type="bed12+")
        gene2iso = {}

        for iso in isoBedObj.reads:
            isoName = "{}_{}".format(i.sample_name, iso)
            isoBedObj.reads[iso].name = isoName
            if isoBedObj.reads[iso].otherList[0] not in gene2iso:
                gene2iso[isoBedObj.reads[iso].otherList[0]] = []
            gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])

        for gene in gene2iso:
            for iso in gene2iso[gene]:
                if iso.chrom + "_" + iso.juncChain not in tmpDict:
                    tmpDict[iso.chrom + "_" + iso.juncChain] = [iso]
                else:
                    tmpDict[iso.chrom + "_" + iso.juncChain].append(iso)

    for tmp in tmpDict:
        isos = tmpDict[tmp]
        mergedIsoName = "+".join([x.name for x in isos])
        sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
        repIso = copy.copy(sortedIsos[0])
        repIso.name = mergedIsoName
        print >> mergedIso2ReadsBed, str(repIso)
    mergedIso2ReadsBed.close()
    return os.path.join(os.getcwd(), "all_sample_merged_iso.bed")


def diff_as(dataToProcess, compCondFile=None, dirSpec=None, sampleMerged=False, args=None, optionTools=None, refInfoParams=None):
    print getCurrentTime() + " Identify differential alternative spliced genes..."
    if validateFile(compCondFile):
        compCondFile = os.path.abspath(compCondFile)
    else:
        compCondFile = compCondFile.split(",")
    prevDir = os.getcwd()
    dasDir = os.path.join(dirSpec.out_dir, "das")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(dasDir)
    samples, conditions = validateSamples(compCondFile, dataToProcess)
    if len(samples) == 0:
        raise Exception("May be you provide wrong sample comparison condition, please check it!")
    mergedIsoBed = mergeIsoforms(samples=samples, dirSpec=dirSpec, pu_filter=args.pu_filter, refInfoParams=refInfoParams)
    cmd = "{}/bed2gpe.pl -g 13 {} > all_sample_merged_iso.gpe".format(utilDir, mergedIsoBed)
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file all_sample_merged_iso.gpe all_sample_merged_iso.gtf"
    subprocess.call(cmd, shell=True)
    gtfFile = os.path.join(os.getcwd(), "all_sample_merged_iso.gtf")

    sample2bamDict = getSample2Bam(samples, compCondFile, dirSpec=dirSpec)


    currentThreads = 1
    for cond1, cond2 in combinations(conditions, 2):
        comp1Name = list(set([x[0] for x in sample2bamDict["cond2bam"][cond1]]))[0]
        comp2Name = list(set([x[0] for x in sample2bamDict["cond2bam"][cond2]]))[0]
        b1File = "{}.txt".format(comp1Name)
        b2File = "{}.txt".format(comp2Name)
        b1List = list(set([x[1] for x in sample2bamDict["cond2bam"][cond1]]))
        b2List = list(set([x[1] for x in sample2bamDict["cond2bam"][cond2]]))
        out1 = open(b1File, "w")
        print >> out1, ",".join(b1List)
        out1.close()
        out2 = open(b2File, "w")
        print >> out2, ",".join(b2List)
        out2.close()

        compOutDir = comp1Name + "_vs_" + comp2Name
        cond1Paired = list(set([x[2] for x in sample2bamDict["cond2bam"][cond1]]))
        cond2Paired = list(set([x[2] for x in sample2bamDict["cond2bam"][cond2]]))
        cond1Length = list(set([x[3] for x in sample2bamDict["cond2bam"][cond1]]))
        cond2Length = list(set([x[3] for x in sample2bamDict["cond2bam"][cond2]]))
        if cond1Paired != cond2Paired or len(cond1Paired) > 1 or len(cond2Paired) > 1 or cond1Length != cond2Length or len(cond1Length) > 1 or len(cond2Length) > 1:
            print "rMATS can't resolve the situation where condition {} and {} with different ngs sequencing strategy or read length!".format(comp1Name, comp2Name)
            continue
        cmd = "rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --tstat {} --nthread {} 1>{}.rmats.log 2>&1"
        cmd = cmd.format(b1File, b2File, gtfFile, compOutDir, cond1Paired[0], cond1Length[0], currentThreads, currentThreads, compOutDir)
        subprocess.call(cmd, shell=True)

        resolveDir("{}.sigDiffAS".format(compOutDir), chdir=False)
        irDiff = pd.read_csv("{}/RI.MATS.JC.txt".format(compOutDir), sep="\t")
        seDiff = pd.read_csv("{}/SE.MATS.JC.txt".format(compOutDir), sep="\t")
        a3ssDiff = pd.read_csv("{}/A3SS.MATS.JC.txt".format(compOutDir), sep="\t")
        a5ssDiff = pd.read_csv("{}/A5SS.MATS.JC.txt".format(compOutDir), sep="\t")
        irDiff = irDiff.loc[irDiff.FDR <= 0.05]
        seDiff = seDiff.loc[seDiff.FDR <= 0.05]
        a3ssDiff = a3ssDiff.loc[a3ssDiff.FDR <= 0.05]
        a5ssDiff = a5ssDiff.loc[a5ssDiff.FDR <= 0.05]

        sigIRfile = "{}.sigDiffAS/IR.sig.txt".format(compOutDir)
        sigSEfile = "{}.sigDiffAS/SE.sig.txt".format(compOutDir)
        sigA3SSfile = "{}.sigDiffAS/A3SS.sig.txt".format(compOutDir)
        sigA5SSfile = "{}.sigDiffAS/A5SS.sig.txt".format(compOutDir)

        irDiff.to_csv(sigIRfile, sep="\t", header=True, index=False)
        seDiff.to_csv(sigSEfile, sep="\t", header=True, index=False)
        a3ssDiff.to_csv(sigA3SSfile, sep="\t", header=True, index=False)
        a5ssDiff.to_csv(sigA5SSfile, sep="\t", header=True, index=False)

        if args.go:
            tmpSigAsFiles = [sigIRfile, sigSEfile, sigA3SSfile, sigA5SSfile]
            cmd = "cat {} | grep -v 'ID' | cut -f 2 | sort -u > {}.sigDiffAS/dasg.lst".format(" ".join(tmpSigAsFiles), compOutDir)
            subprocess.call(cmd, shell=True)
            from plotRscriptStrs import plotTargetGenesGoEnrichmentStr
            # outName = compPair
            gene2goFile = args.gene2goFile if args.gene2goFile else optionTools.gene2go
            if not gene2goFile:
                print "You don't provide gene2go file, the GO enrichment will not be carried out!"
                continue

            from rpy2 import robjects
            from rpy2.rinterface import RRuntimeWarning
            import warnings
            warnings.filterwarnings("ignore", category=RRuntimeWarning)
            robjects.r(plotTargetGenesGoEnrichmentStr)
            robjects.r.plotTargetGenesGoEnrichment("{}.sigDiffAS/dasg.lst".format(compOutDir), compOutDir, gene2goFile,
                                                   "{}.sigDiffAS/sigDiff".format(compOutDir), float(args.cutoff),
                                                   args.filterBy, int(args.showCategory))
        # enrichResult = os.path.abspath("sigDiff.goEnrichResults.txt")
        # enrichPlot = convertPdf2png(inPdf=os.path.abspath(outName + ".pdf"))
        # resultDict["das"][compPair].update({"goEnrichResults": enrichResult})
        # resultDict["das"][compPair].update({"goEnrichPlot": enrichPlot})

    os.chdir(prevDir)
    print getCurrentTime() + " Identify differential alternative spliced genes done!"
