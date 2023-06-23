#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: asas.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-19
Last modified: 2022-01-19
'''

import itertools, glob, StringIO, copy, vcf, pysam
import pandas as pd
from scipy.stats import chi2_contingency
from ConfigParser import RawConfigParser
from commonFuncs import *
from commonObjs import *

min_bq = 13
ploidy = 2


def assignHaplo(haplotypes, vcfHaploDict):
    from difflib import SequenceMatcher
    haploAssigned = {}
    for haplo in haplotypes:
        tmpHaplo = ""
        tmpRatio = 0
        for vcfHaplo in vcfHaploDict:
            if SequenceMatcher(None, haplo, vcfHaplo).ratio() > tmpRatio:
                tmpHaplo = vcfHaplo
                tmpRatio = SequenceMatcher(None, haplo, vcfHaplo).ratio()
        haploAssigned[haplo] = vcfHaploDict[tmpHaplo]
    return haploAssigned


def getASpairedIsoforms(asFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=0, mergeByJunc=True):
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    isoBedObj = BedFile(isoformFile, type="bed12+")
    asPairs = {}
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            asEvent = records[3]
            if mergeByJunc:
                incJuncCombDict = {}
                excJuncCombDict = {}
                for incIso in inclusionIsos:
                    incIsoObj = isoBedObj.reads[incIso]
                    if len(incIsoObj.introns) > 1:
                        if incIsoObj.juncChain not in incJuncCombDict:
                            incJuncCombDict[incIsoObj.juncChain] = [incIso]
                        else:
                            incJuncCombDict[incIsoObj.juncChain].append(incIso)
                    else:
                        if "monoExon" not in incJuncCombDict:
                            incJuncCombDict["monoExon"] = [incIso]
                        else:
                            incJuncCombDict["monoExon"].append(incIso)
                for excIso in exclusionIsos:
                    excIsoObj = isoBedObj.reads[excIso]
                    if len(excIsoObj.introns) > 1:
                        if excIsoObj.juncChain not in excJuncCombDict:
                            excJuncCombDict[excIsoObj.juncChain] = [excIso]
                        else:
                            excJuncCombDict[excIsoObj.juncChain].append(excIso)
                    else:
                        if "monoExon" not in excJuncCombDict:
                            excJuncCombDict["monoExon"] = [excIso]
                        else:
                            excJuncCombDict["monoExon"].append(excIso)

                for item in itertools.product(incJuncCombDict.values(), excJuncCombDict.values()):
                    newInclusionIsos = []
                    newExclusionIsos = []
                    for x in item[0]:
                        countInc = sum([len(collapsedTrans2reads[y]) for y in isoBedObj.reads[x].otherList[2].split(",")])
                        if countInc >= filterByCount:
                            newInclusionIsos.append(x)
                    for x in item[1]:
                        countExc = sum([len(collapsedTrans2reads[y]) for y in isoBedObj.reads[x].otherList[2].split(",")])
                        if countExc >= filterByCount:
                            newExclusionIsos.append(x)
                    # newInclusionIsos = [x for x in item[0] if len(collapsedTrans2reads[x]) >= filterByCount]
                    # newExclusionIsos = [x for x in item[1] if len(collapsedTrans2reads[x]) >= filterByCount]
                    if len(newInclusionIsos) == 0 or len(newExclusionIsos) == 0:
                        continue
                    '''the isoform PB.x.x split by "." to determine gene PB.x'''
                    inclusionGene = [".".join(x.split(".")[:2]) for x in newInclusionIsos]
                    exclusionGene = [".".join(x.split(".")[:2]) for x in newExclusionIsos]
                    uniqGene = list(set(inclusionGene+exclusionGene))
                    if len(uniqGene) == 1:
                        if uniqGene[0] not in asPairs:
                            asPairs[uniqGene[0]] = {asEvent: [newInclusionIsos, newExclusionIsos]}
                        else:
                            asPairs[uniqGene[0]].update({asEvent: [newInclusionIsos, newExclusionIsos]})
            else:
                newInclusionIsos = []
                newExclusionIsos = []
                for x in inclusionIsos:
                    countInc = sum([len(collapsedTrans2reads[y]) for y in isoBedObj.reads[x].otherList[2].split(",")])
                    if countInc >= filterByCount:
                        newInclusionIsos.append(x)
                for x in exclusionIsos:
                    countExc = sum([len(collapsedTrans2reads[y]) for y in isoBedObj.reads[x].otherList[2].split(",")])
                    if countExc >= filterByCount:
                        newExclusionIsos.append(x)
                # newInclusionIsos = [x for x in inclusionIsos if len(collapsedTrans2reads[x]) >= filterByCount]
                # newExclusionIsos = [x for x in exclusionIsos if len(collapsedTrans2reads[x]) >= filterByCount]
                if len(newInclusionIsos) == 0 or len(newExclusionIsos) == 0:
                    continue
                inclusionGene = [".".join(x.split(".")[:2]) for x in newInclusionIsos]
                exclusionGene = [".".join(x.split(".")[:2]) for x in newExclusionIsos]
                uniqGene = list(set(inclusionGene + exclusionGene))
                if len(uniqGene) == 1:
                    if uniqGene[0] not in asPairs:
                        asPairs[uniqGene[0]] = {asEvent: [newInclusionIsos, newExclusionIsos]}
                    else:
                        asPairs[uniqGene[0]].update({asEvent: [newInclusionIsos, newExclusionIsos]})

        return asPairs


def generate_regions(fasta_index_file, size, chunks=False, chromosomes=None, outFile=None):
    from math import ceil
    if not fasta_index_file.endswith(".fai"):
        fasta_index_file = fasta_index_file + ".fai"

    out = open(outFile, "w")
    with open(fasta_index_file, "r") as fasta_index:
        for line in fasta_index:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])
            if chromosomes is not None and chrom_name not in chromosomes:
                continue
            region_start = 0
            if chunks is True:
                region_size = ceil(chrom_length / size)  # have to make sure this works
            else:
                region_size = size
            while region_start < chrom_length:
                region_end = region_start + region_size
                if region_end > chrom_length:
                    region_end = chrom_length
                start = str(region_start)
                end = str(region_end)
                print >> out, "{}:{}-{}".format(chrom_name, start, end)
                region_start = region_end


def getHaploIsoformSeq(lociDir, isoformFile, refGenome, refFa=None, altFa=None):
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    representIsoOut = open("representIso.bed", "w")
    haplo2flncCount = {}
    for i in lociDir:
        lociName = os.path.basename(i)
        partialHaplotype = os.path.join(i, "phased.partial.cleaned.human_readable.txt")
        cleanVCF = os.path.join(i, "phased.partial.cleaned.vcf")
        if os.path.exists(partialHaplotype):
            partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
            vcfHaploDict = vcfHaplo(cleanVCF)
            haplotypes = list(partialDF.index)
            # isoforms = list(partialDF.columns)
            haploAssigned = assignHaplo(haplotypes, vcfHaploDict)
            for haplo in haplotypes:
                junc2iso = {}
                isoformAssigned = partialDF.iloc[:, partialDF.loc[haplo, ].to_numpy().nonzero()[0]].columns
                for j in isoformAssigned:
                    if j not in isoformBed: continue
                    if len(isoformBed[j].introns) > 1:
                        if isoformBed[j].juncChain not in junc2iso:
                            junc2iso[isoformBed[j].juncChain] = {"iso": [j]}
                        else:
                            junc2iso[isoformBed[j].juncChain]["iso"].append(j)
                    else:
                        if "monoExon" not in junc2iso:
                            junc2iso["monoExon"] = {"iso": [j]}
                        else:
                            junc2iso["monoExon"]["iso"].append(j)
                tmpJunc = ""
                tmpCount = 0
                for y in junc2iso:
                    if partialDF.loc[haplo, junc2iso[y]["iso"]].sum() > tmpCount:
                        tmpJunc = y
                        tmpCount = partialDF.loc[haplo, junc2iso[y]["iso"]].sum()
                repIsos = [isoformBed[x] for x in junc2iso[tmpJunc]["iso"]]
                sortedRepIsos = sorted(repIsos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                repIso = copy.copy(sortedRepIsos[0])
                repIso.name = "{}_{}".format(lociName, haploAssigned[haplo])
                haplo2flncCount[repIso.name] = [tmpCount, repIso.otherList[0]]
                print >> representIsoOut, str(repIso)
    representIsoOut.close()
    cmd = '''cut -f 1-12 representIso.bed | bedtools getfasta -fi {} -bed - -name -split | seqkit replace -w 0 -p "(.*?):(.*)" -r '$1' | seqkit rmdup -n > representIso.fa'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    if not refFa or not validateFile(refFa):
        cmd = '''seqkit grep -p "REF$" -r representIso.fa > representIso.REF.fa'''
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        refFa = os.path.join(os.getcwd(), "representIso.REF.fa")
    if not altFa or not validateFile(altFa):
        cmd = '''seqkit grep -p "ALT$" -r representIso.fa > representIso.ALT.fa'''
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        altFa = os.path.join(os.getcwd(), "representIso.ALT.fa")
    return refFa, altFa, haplo2flncCount


def tagBamAndFilter(snpPosFile, inBam, taggedBam, utilDir):
    assignReadsToRef = os.path.join(utilDir, "biostar214299.jar")
    cmd = "java -jar {} -p {} {} -o {}".format(assignReadsToRef, snpPosFile, inBam, taggedBam)
    subprocess.call(cmd, shell=True)
    refReads = []
    altReads = []
    samfile = pysam.AlignmentFile(taggedBam, "rb")
    for record in samfile:
        if record.get_tag("RG") == "REF" and record.get_tag("NH") == 1:
            refReads.append(record.query_name)
        if record.get_tag("RG") == "ALT" and record.get_tag("NH") == 1:
            altReads.append(record.query_name)
    samfile.close()
    return refReads, altReads


def vcfHaplo(cleanVCF):
    # import vcf
    vcf_reader = vcf.Reader(open(cleanVCF, "r"))
    refStr, altStr = "", ""
    for record in vcf_reader:
        refStr += record.REF
        altStr += str(record.ALT[0])
    return {refStr: "REF", altStr: "ALT"}


def makeAbundanceFile(groupFile, outFile=None):
    with open(groupFile) as f:
        cidInfo = {}
        for line in f.readlines():
            pbid, members = line.strip().split('\t')
            for cid in members.split(','):
                cidInfo[cid] = pbid
        if outFile:
            out = open(outFile, "w")
            print >> out, "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print >> out, "\t".join([i, "NA", "Y", "unique", cidInfo[i]])
            out.close()
        else:
            print "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print "\t".join([i, "NA", "Y", "unique", cidInfo[i]])


def mergeASEfromSalmonOut(ref_quant_sf, alt_quant_sf, haplo2count, mergedASEfile, calASE=False):
    ref_quant = pd.read_csv(ref_quant_sf, sep="\t")
    alt_quant = pd.read_csv(alt_quant_sf, sep="\t")
    aseDict = {}
    for row in ref_quant.iteritems():
        gene = row.Name.split("_")[0]
        if gene not in aseDict:
            aseDict[gene] = {"tgs": {"REF": 0, "ALT": 0}, "ngs": {"REF": [0, 0], "ALT": [0, 0]}, "refGene": ""}
        aseDict[gene]["tgs"]["REF"] = float(haplo2count[row.Name][0])
        aseDict[gene]["ngs"]["REF"] = [float(row.TPM), float(row.NumReads)]
        aseDict[gene]["refGene"] = haplo2count[row.Name][1]
    for row in alt_quant.iteritems():
        gene = row.Name.split("_")[0]
        if gene not in aseDict:
            aseDict[gene] = {"tgs": {"REF": 0, "ALT": 0}, "ngs": {"REF": [0, 0], "ALT": [0, 0]}, "refGene": ""}
        aseDict[gene]["tgs"]["ALT"] = float(haplo2count[row.Name][0])
        aseDict[gene]["ngs"]["ALT"] = [float(row.TPM), float(row.NumReads)]
        aseDict[gene]["refGene"] = haplo2count[row.Name][1]

    out = open(mergedASEfile, "w")
    for gene in aseDict:
        refGene = aseDict[gene]["refGene"]
        tgsRefCount = aseDict[gene]["tgs"]["REF"]
        tgsAltCount = aseDict[gene]["tgs"]["ALT"]
        ngsRefTPM, ngsRefCount = aseDict[gene]["ngs"]["REF"]
        ngsAltTPM, ngsAltCount = aseDict[gene]["ngs"]["ALT"]
        if calASE:
            tgsPvalue = chi2_contingency([[tgsRefCount, tgsAltCount], [(tgsRefCount+tgsAltCount)/2, (tgsRefCount+tgsAltCount)/2]])[1]
            ngsPvalue = chi2_contingency([[ngsRefCount, ngsAltCount], [(ngsRefCount+ngsAltCount)/2, (ngsRefCount+ngsAltCount)/2]])[1]
            if tgsPvalue > 0.05 or ngsPvalue > 0.05:
                continue
        print >> out, "\t".join([gene, refGene, tgsRefCount, tgsAltCount, ngsRefTPM, ngsRefCount, ngsAltTPM, ngsAltCount])
    out.close()


def runPhaser(targetDir):
    curDir = os.getcwd()
    os.chdir(targetDir)
    configFile = "config"
    with open(configFile) as f:
        config_string = StringIO.StringIO("[dummy_section]\n" + f.read())
    config_parser = RawConfigParser()
    config_parser.readfp(config_string)
    strand = config_parser.get("dummy_section", "ref_strand")

    cmd = "minimap2 -ax splice fake.fasta ccs.fastq >ccs.sam 2>/dev/null"
    subprocess.call(cmd, shell=True)

    cmd = "samtools sort ccs.sam > ccs.sorted.bam 2>/dev/null; samtools mpileup --min-BQ {} -f fake.fasta -s ccs.sorted.bam > ccs.mpileup 2>/dev/null".format(min_bq)
    subprocess.call(cmd, shell=True)
    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} " \
          "-o phased.partial -n {} --bhFDR 0.01 -p 0.01 1>partial.stdout.txt 2>&1".format(strand, ploidy)
    subprocess.call(cmd, shell=True)

    os.chdir(curDir)


def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    resultDict = {"partial": {}, "nopartial": {}}
    for asEvent in isoPairs:
        isoPair = isoPairs[asEvent]
        if partial:
            partialHaplotype = "phased.partial.cleaned.human_readable.txt"
            if os.path.exists(partialHaplotype):
                partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
                newIncIsos = list(set(partialDF.columns) & set(isoPair[0]))
                newExcIsos = list(set(partialDF.columns) & set(isoPair[1]))
                if len(newIncIsos) == 0 or len(newExcIsos) == 0: continue
                incIsoDf = partialDF.loc[:, newIncIsos].sum(axis=1)
                excIsoDf = partialDF.loc[:, newExcIsos].sum(axis=1)
                mergedDf = pd.concat([incIsoDf, excIsoDf], axis=1)
                if (pd.DataFrame.max(mergedDf) >= 5).all() and (pd.DataFrame.min(mergedDf) <= 3).all():
                    if len(mergedDf[(mergedDf.T != 0).any()]) == 1: continue
                    chi2Result = chi2_contingency(mergedDf)
                    chi2Pvalue = chi2Result[1]
                    if chi2Pvalue <= 0.001:
                        haploInc = incIsoDf.idxmax(axis=0)
                        haploExc = excIsoDf.idxmax(axis=0)
                        if haploInc == haploExc: continue
                        resultDict["partial"].update(
                            {asEvent: dict(zip([haploInc, haploExc], ["_".join(isoPair[0]), "_".join(isoPair[1])]))})

        if nopartial:
            nopartialHaplotype = "phased.nopartial.cleaned.human_readable.txt"
            if os.path.exists(nopartialHaplotype):
                nopartialDF = pd.read_csv(nopartialHaplotype, skiprows=[0], index_col=0, sep="\t")
                # if not set(isoPair[0] + isoPair[1]).issubset(set(nopartialDF.columns)): continue
                newIncIsos = list(set(nopartialDF.columns) & set(isoPair[0]))
                newExcIsos = list(set(nopartialDF.columns) & set(isoPair[1]))
                if len(newIncIsos) == 0 or len(newExcIsos) == 0: continue
                incIsoDf = nopartialDF.loc[:, newIncIsos].sum(axis=1)
                excIsoDf = nopartialDF.loc[:, newExcIsos].sum(axis=1)
                mergedDf = pd.concat([incIsoDf, excIsoDf], axis=1)
                if (pd.DataFrame.max(mergedDf) >= 5).all() and (pd.DataFrame.min(mergedDf) <= 3).all():
                    if len(mergedDf[(mergedDf.T != 0).any()]) == 1: continue
                    chi2Result = chi2_contingency(mergedDf)
                    chi2Pvalue = chi2Result[1]
                    if chi2Pvalue <= 0.001:
                        haploInc = incIsoDf.idxmax(axis=0)
                        haploExc = excIsoDf.idxmax(axis=0)
                        if haploInc == haploExc: continue
                        resultDict["nopartial"].update({asEvent: dict(
                            zip([haploInc, haploExc], [["_".join(isoPair[0])], ["_".join(isoPair[1])]]))})

        if not partial and not nopartial:
            raise Exception("You must specify either partial or nopartial file")
    return resultDict

def asas(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify allele-specific AS events for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "asas"))

    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    processedFq = os.path.join(baseDir, "preprocess", dataObj.tgs_plat.lower(), "rawFlnc.fq")
    if not validateFile(processedFq):
        processedFq = os.path.join(baseDir, "preprocess", "fmlrc", "rawFlnc.fq")
    collapsedGff = os.path.join(baseDir, "collapse", "tofu.collapsed.good.gff")
    # isoformFile = os.path.join(baseDir, "refine", "isoformGrouped.bed12+")
    isoformFile = os.path.join(baseDir, "hqIsoforms", "hq.collapsed.bed12+")

    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "SE.confident.bed12+")
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "IR.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A3SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=2, mergeByJunc=False)
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, isoformFile, asType="IR", filterByCount=2, mergeByJunc=False)
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, isoformFile, asType="A5SS", filterByCount=2, mergeByJunc=False)
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, isoformFile, asType="A3SS", filterByCount=2, mergeByJunc=False)

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    cmd = "rm -R by_loci"
    subprocess.call(cmd, shell=True)
    cmd = "select_loci_to_phase.py {} {} {} {} -c 10 1>/dev/null 2>&1"
    cmd = cmd.format(refParams.ref_genome, processedFq, collapsedGff, readStatFile)
    subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    pool = Pool(processes=dataObj.single_run_threads)
    for i in lociDir:
        pool.apply_async(runPhaser, (i,))
    pool.close()
    pool.join()

    resultDict = {}
    for i in lociDir:
        prevDir = os.getcwd()
        os.chdir(i)
        fakeGene = os.path.basename(i).split("_")[0]
        resultDict[fakeGene] = {}
        for asType in asPairs:
            if fakeGene not in asPairs[asType]:
                continue
            isoPairs = asPairs[asType][fakeGene]
            sigRelatedDict = relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=False, isoPairs=isoPairs)
            resultDict[fakeGene].update({asType: sigRelatedDict})
        os.chdir(prevDir)

    '''the structure of resultDict is resultDict = {fakeGene: {asType: {partialCategory: {asEvent: {haplotype1: PB.X.1, haplotype2: PB.X.2, ...}}}}}'''
    partialOut = open("partialAsRelatedHaplotype.txt", "w")
    for gene in resultDict:
        for asType in resultDict[gene]:
            for partial in resultDict[gene][asType]:
                for asEvent in resultDict[gene][asType][partial]:
                    haplos = resultDict[gene][asType][partial][asEvent]
                    haplotype1, haplotype2 = haplos.keys()
                    refGene = "_".join(asEvent.split(":")[0:-1])
                    print >>partialOut, "\t".join([gene, asType, refGene, asEvent, haplotype1, haplos[haplotype1], haplotype2, haplos[haplotype2]])
    partialOut.close()
    os.chdir(prevDir)
    print getCurrentTime() + " Identify allele-specific AS events for project {} sample {} done!".format(projectName, sampleName)

def ase(dataObj=None, refParams=None, dirSpec=None, refFa=None, altFa=None, useFreeBayes=False):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify allelic-specific expression genes for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "ase")
    logDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "log")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(baseDir)

    lociDir = glob.glob(os.path.join(dirSpec.out_dir, projectName, sampleName, "asas", "by_loci/*size*"))
    freeBayesAlleleSNP = {}
    hybridBam = os.path.join(dirSpec.out_dir, projectName, sampleName, "mapping", "rna-seq", "reassembly", "tmp.bam")
    if useFreeBayes:
        generate_regions(refParams.ref_genome, 1000000, outFile="split.region.txt")
        cmd = "freebayes-parallel <(cat split.region.txt) 40 -f {} {} -C 10 -F 0.2 --min-base-quality 20 > freebayes.SNP.vcf".format(
            refParams.ref_genome, hybridBam)
        subprocess.call(cmd, shell=True, executable="/bin/bash")


        vcf_reader = vcf.Reader(open("freebayes.SNP.vcf", "r"))
        for record in vcf_reader:
            snpName = "{}_{}".format(record.CHROM, record.POS)
            if record.INFO["AF"] == 0.5 and len(record.INFO["TYPE"]) == 1 and record.INFO["TYPE"][0] == "snp":
                freeBayesAlleleSNP[snpName] = record

    snp_position = open("snp_position.txt", "w")
    snpDict = {}
    for i in lociDir:
        partialHaplotype = os.path.join(i, "phased.partial.cleaned.human_readable.txt")
        partialVcf = os.path.join(i, "phased.partial.cleaned.vcf")
        if os.path.exists(partialHaplotype):
            haplo_vcf = vcf.Reader(open(partialVcf, "r"))
            for record in haplo_vcf:
                snpName = "{}_{}".format(record.CHROM, record.POS)
                if snpName in snpDict: continue
                snpDict[snpName] = ""
                if len(record.REF) != 1 or len(record.ALT) != 1: continue
                if len(freeBayesAlleleSNP) != 0 and snpName not in freeBayesAlleleSNP: continue
                print >> snp_position, "\t".join(map(str, [record.CHROM, record.POS, record.REF, "REF"]))
                print >> snp_position, "\t".join(map(str, [record.CHROM, record.POS, record.ALT[0], "ALT"]))
    snp_position.close()

    cmd = "bamtools split -in {} -reference".format(hybridBam)
    subprocess.call(cmd, shell=True)
    splitBams = glob.glob(os.path.join(os.path.dirname(hybridBam), "tmp.REF_*.bam"))
    taggedBams = []
    resultList = []
    pool = Pool(processes=dataObj.single_run_threads)
    for i in splitBams:
        newBam = os.path.join(os.path.dirname(i), "tagged." + os.path.basename(i))
        taggedBams.append(newBam)
        tmpRes = pool.apply_async(tagBamAndFilter, ("snp_position.txt", i, newBam, utilDir))
        resultList.append(tmpRes)
    pool.close()
    pool.join()

    # refReads = []
    # altReads = []
    # for res in resultList:
    #     tmp = res.get()
    #     refReads.extend(tmp[0])
    #     altReads.extend(tmp[1])
    # with open("ref_reads.lst", "w") as f:
    #     for i in refReads:
    #         print >> f, i
    # with open("alt_reads.lst", "w") as f:
    #     for i in altReads:
    #         print >> f, i
    # cmd = "samtools cat {} > tagged.bam".format(" ".join(taggedBams))
    # subprocess.call(cmd, shell=True)
    refReads = open("ref_reads.lst", "w")
    altReads = open("alt_reads.lst", "w")
    for i in taggedBams:
        samfile = pysam.AlignmentFile(i, "rb")
        for record in samfile:
            if record.get_tag("RG") == "REF" and record.get_tag("NH") == 1:
                print >> refReads, record.query_name
            if record.get_tag("RG") == "ALT" and record.get_tag("NH") == 1:
                print >> altReads, record.query_name
        samfile.close()
    refReads.close()
    altReads.close()
    removeFiles(fileList=splitBams)
    removeFiles(fileList=taggedBams)

    ##################
    isoformFile = os.path.join(dirSpec.out_dir, projectName, sampleName, "refine", "isoformGrouped.bed12+")
    refFa, altFa, haplo2flncCount = getHaploIsoformSeq(lociDir, isoformFile, refParams.ref_genome, refFa=refFa, altFa=altFa)
    cmd = "salmon index -t {} -i ref_salmon_index --keepDuplicates -p {} 2>{}/ase.ref_index.log".format(refFa, dataObj.single_run_threads, logDir)
    subprocess.call(cmd, shell=True)
    cmd = "salmon index -t {} -i alt_salmon_index --keepDuplicates -p {} 2>{}/ase.alt_index.log".format(altFa, dataObj.single_run_threads, logDir)
    subprocess.call(cmd, shell=True)

    if dataObj.ngs_reads_paired == "paired":
        leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
        rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                refLeftRepeatFq = "{}.ref.left_reads.fq.gz"
                refRightRepeatFq = "{}.ref.right_reads.fq.gz"
                altLeftRepeatFq = "{}.alt.left_reads.fq.gz"
                altRightRepeatFq = "{}.alt.right_reads.fq.gz"
                cmd = "seqkit grep -f ref_reads.lst {} > {}".format(" ".join(leftReads), refLeftRepeatFq)
                subprocess.call(cmd, shell=True)
                cmd = "seqkit grep -f ref_reads.lst {} > {}".format(" ".join(rightReads), refRightRepeatFq)
                subprocess.call(cmd, shell=True)
                cmd = "seqkit grep -f alt_reads.lst {} > {}".format(" ".join(leftReads), altLeftRepeatFq)
                subprocess.call(cmd, shell=True)
                cmd = "seqkit grep -f alt_reads.lst {} > {}".format(" ".join(rightReads), altRightRepeatFq)
                subprocess.call(cmd, shell=True)

                refRepeatOut = "repeat{}.ref_salmon_quant".format(i)
                cmd = "salmon quant -l A -i ref_salmon_index -1 {} -2 {} -p {} -o {}"
                cmd = cmd.format(refLeftRepeatFq, refRightRepeatFq, dataObj.single_run_threads, refRepeatOut)
                subprocess.call(cmd, shell=True)

                altRepeatOut = "repeat{}.alt_salmon_quant".format(i)
                cmd = "salmon quant -l A -i alt_salmon_index -1 {} -2 {} -p {} -o {}"
                cmd = cmd.format(altLeftRepeatFq, altRightRepeatFq, dataObj.single_run_threads, altRepeatOut)
                subprocess.call(cmd, shell=True)

                ref_quant_sf = os.path.join(refRepeatOut, "quant_sf")
                alt_quant_sf = os.path.join(altRepeatOut, "quant_sf")
                mergedASEfile = "mergedASEfile.txt"
                mergeASEfromSalmonOut(ref_quant_sf, alt_quant_sf, haplo2flncCount, mergedASEfile)
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                refLeftRepeatFq = "{}.ref.left_reads.fq.gz"
                altLeftRepeatFq = "{}.alt.left_reads.fq.gz"
                cmd = "seqkit grep -f ref_reads.lst {} > {}".format(" ".join(leftReads), refLeftRepeatFq)
                subprocess.call(cmd, shell=True)
                cmd = "seqkit grep -f alt_reads.lst {} > {}".format(" ".join(leftReads), altLeftRepeatFq)
                subprocess.call(cmd, shell=True)

                refRepeatOut = "repeat{}.ref_salmon_quant".format(i)
                cmd = "salmon quant -l A -i ref_salmon_index -r {} -p {} -o {}"
                cmd = cmd.format(refLeftRepeatFq, dataObj.single_run_threads, refRepeatOut)
                subprocess.call(cmd, shell=True)

                altRepeatOut = "repeat{}.alt_salmon_quant".format(i)
                cmd = "salmon quant -l A -i alt_salmon_index -r {} -p {} -o {}"
                cmd = cmd.format(altLeftRepeatFq, dataObj.single_run_threads, altRepeatOut)
                subprocess.call(cmd, shell=True)

                ref_quant_sf = os.path.join(refRepeatOut, "quant_sf")
                alt_quant_sf = os.path.join(altRepeatOut, "quant_sf")
                mergedASEfile = "mergedASEfile.txt"
                mergeASEfromSalmonOut(ref_quant_sf, alt_quant_sf, haplo2flncCount, mergedASEfile)
        elif dataObj.ngs_right_reads and dataObj.ngs_left_reads == None:
            rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                refLeftRepeatFq = "{}.ref.right_reads.fq.gz"
                altLeftRepeatFq = "{}.alt.right_reads.fq.gz"
                cmd = "seqkit grep -f ref_reads.lst {} > {}".format(" ".join(rightReads), refLeftRepeatFq)
                subprocess.call(cmd, shell=True)
                cmd = "seqkit grep -f alt_reads.lst {} > {}".format(" ".join(rightReads), altLeftRepeatFq)
                subprocess.call(cmd, shell=True)

                refRepeatOut = "repeat{}.ref_salmon_quant".format(i)
                cmd = "salmon quant -l A -i ref_salmon_index -r {} -p {} -o {}"
                cmd = cmd.format(refLeftRepeatFq, dataObj.single_run_threads, refRepeatOut)
                subprocess.call(cmd, shell=True)

                altRepeatOut = "repeat{}.alt_salmon_quant".format(i)
                cmd = "salmon quant -l A -i alt_salmon_index -r {} -p {} -o {}"
                cmd = cmd.format(altLeftRepeatFq, dataObj.single_run_threads, altRepeatOut)
                subprocess.call(cmd, shell=True)

                ref_quant_sf = os.path.join(refRepeatOut, "quant_sf")
                alt_quant_sf = os.path.join(altRepeatOut, "quant_sf")
                mergedASEfile = "mergedASEfile.txt"
                mergeASEfromSalmonOut(ref_quant_sf, alt_quant_sf, haplo2flncCount, mergedASEfile)
        else:
            raise Exception("The NGS data seem not to be single, please check it")

    print getCurrentTime() + " Identify allelic-specific expression genes for project {} sample {} done!".format(projectName, sampleName)
    os.chdir(prevDir)
