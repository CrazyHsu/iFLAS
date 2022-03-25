#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: tissue_spec_iso.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-22
Last modified: 2022-01-22
'''

import itertools
import numpy as np
from Bio import SeqIO
from heapq import nlargest
from scipy import stats
from collections import Counter
from commonFuncs import *
from commonObjs import *

def calVC(myList):
    return np.std(myList)/np.mean(myList)


def fastxid2sample(fastx, sample):
    myDict = {}
    if fastx.endswith(("fq", "fastq")):
        suffix = "fastq"
    elif fastx.endswith(("fa", "fasta")):
        suffix = "fasta"
    else:
        raise Exception("")
    for record in SeqIO.parse(fastx, suffix):
        myDict[record.id] = sample
    return myDict


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


def tissue_spec_iso(dataObj, rawDataObjs, dirSpec, isoformBed, isoDict):
    sampleNames = [x.sample_name for x in rawDataObjs]
    print getCurrentTime() + " Start finding tissue specific isoforms among".format(",".join(sampleNames))
    baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
    readsGroup = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")

    iso2reads = {}
    with open(readsGroup) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            iso = infoList[0]
            reads = infoList[1].split(",")
            iso2reads[iso] = reads

    sample2dataDict = {}
    mergeDict = {}
    for i in range(len(sampleNames)):
        dataLocation = rawDataObjs[i].data_processed_location
        sample2dataDict[sampleNames[i]] = dataLocation
        mergeDict.update(fastxid2sample(dataLocation, sampleNames[i]))

    juncDict = {}
    for iso in isoDict:
        if isoformBed[iso].juncChain not in juncDict:
            juncDict[isoformBed[iso].juncChain] = {"gene": isoformBed[iso].otherList[0], "iso": [iso]}
        else:
            juncDict[isoformBed[iso].juncChain]["iso"].append(iso)

    out = open("tissue_spec.txt", "w")
    print >> out, "\t".join(["gene", "junction", "\t".join(sampleNames), "chi2_pvalue", "vc", "category"])
    for junc in juncDict:
        reads = list(itertools.chain.from_iterable([iso2reads[x] for x in juncDict[junc]["iso"]]))
        if len(reads) < 5: continue
        readsAssigned = [(i, mergeDict[i]) for i in reads]
        assignedCount = Counter([x[1] for x in readsAssigned])

        allSampleCounts = []
        for sample in sampleNames:
            if sample in assignedCount:
                allSampleCounts.append(assignedCount[sample])
            else:
                allSampleCounts.append(0)
        meanCount = sum(allSampleCounts) / float(len(allSampleCounts))
        g, p, dof, expctd = stats.chi2_contingency([allSampleCounts, [meanCount] * len(allSampleCounts)])
        vc = calVC(allSampleCounts)
        if p < 0.05 and min(allSampleCounts) >= 5:
            top1, top2 = nlargest(2, allSampleCounts)
            if top1 >= top2 * 1.5 or vc >= 0.25:
                category = sampleNames[allSampleCounts.index(top1)]
                print >> out, "\t".join(map(str, [juncDict[junc]["gene"], junc, "\t".join(map(str, allSampleCounts)), p, vc, category]))
    out.close()

    print getCurrentTime() + " End finding tissue specific isoforms among".format(",".join(sampleNames))