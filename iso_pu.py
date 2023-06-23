#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iso_pu.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2022-04-23 14:22:38
Last modified: 2022-04-23 14:22:38
'''

import subprocess, os, pybedtools, pysam, itertools, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.patches as patches

# from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import Counter
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

from baggingPU import BaggingClassifierPU
from commonFuncs import *
from commonObjs import *

warnings.filterwarnings("ignore")


class SchemeModelEval(object):
    def __init__(self, trainingData=None, model=None, scheme=None, featureData=None):
        self.trainingData = trainingData
        self.model = model
        self.scheme = scheme

        self.finalEstimator = None
        self.predResults = None
        self.f1_score = None
        self.precision = None
        self.recall = None

        self.featureData = featureData

    def eval(self):
        X_train = self.trainingData.iloc[:, :-1]
        y_train = self.trainingData.iloc[:, -1]

        classifier = self.getClassifier()
        if self.scheme == "bagging":
            X = X_train
            y = y_train
            y = pd.Series(y)

            labels = Counter(y)
            sorted_labels = sorted(labels.items(), key=lambda x: x[1], reverse=True)
            max_samples = sorted_labels[0][1] if sorted_labels[0][0] == 0 else sorted_labels[1][1]

            # from baggingPU import BaggingClassifierPU
            cf = BaggingClassifierPU(
                classifier,
                n_estimators=100,
                max_samples=max_samples,
                n_jobs=-1
            )
            cf.fit(X, y)

            pu_score = cf.oob_decision_function_[:, 1]
            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": pu_score
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf
        elif self.scheme == "2step":
            X = X_train
            y = y_train
            y = pd.Series(y)

            cf = classifier
            cf.fit(X, y)

            ys = 2*y -1
            pred = cf.predict_proba(X)[:, 1]

            # Find the range of scores given to positive data points
            range_P = [min(pred * (ys > 0)), max(pred * (ys > 0))]

            # STEP 1
            # If any unlabeled point has a score above all known positives,
            # or below all known positives, label it accordingly
            iP_new = ys[(ys < 0) & (pred >= range_P[1])].index
            iN_new = ys[(ys < 0) & (pred <= range_P[0])].index
            ys.loc[iP_new] = 1
            ys.loc[iN_new] = 0

            # Classifier to be used for step 2
            cf2 = classifier

            # Limit to 10 iterations (this is arbitrary, but
            # otherwise this approach can take a very long time)
            for i in range(10):
                # If step 1 didn't find new labels, we're done
                if len(iP_new) + len(iN_new) == 0 and i > 0:
                    break

                # STEP 2
                # Retrain on new labels and get new scores
                cf2.fit(X, ys)
                pred = cf2.predict_proba(X)[:, -1]

                # Find the range of scores given to positive data points
                range_P = [min(pred * (ys > 0)), max(pred * (ys > 0))]

                # Repeat step 1
                iP_new = ys[(ys < 0) & (pred >= range_P[1])].index
                iN_new = ys[(ys < 0) & (pred <= range_P[0])].index
                ys.loc[iP_new] = 1
                ys.loc[iN_new] = 0

            # Lastly, get the scores assigned by this approach
            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": pred
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf2
        else:
            X = X_train
            y = y_train
            y = pd.Series(y)

            cf = classifier
            cf.fit(X, y)

            self.predResults = pd.DataFrame({
                "label": y,
                "pu_score": cf.predict_proba(X)[:, 1]
            }, columns=["label", "pu_score"])
            self.finalEstimator = cf

        # self.cfTest(self.finalEstimator, X_test, y_test)
        # self.filterIsoformsByScore(filterScore=0.8, outFile="validIsoforms.lst")

    def getClassifier(self):
        if self.model == "SVM":
            from sklearn.svm import SVC
            return SVC(kernel='rbf', gamma='auto', random_state=0)
        elif self.model == "DT":
            from sklearn.tree import DecisionTreeClassifier
            return DecisionTreeClassifier()
        elif self.model == "RF":
            from sklearn.ensemble import RandomForestClassifier
            return RandomForestClassifier(n_estimators=10, random_state=0)
        elif self.model == "GB":
            from sklearn.ensemble import GradientBoostingClassifier
            return GradientBoostingClassifier(n_estimators=100, random_state=0)
        elif self.model == "NB":
            from sklearn.naive_bayes import GaussianNB
            return GaussianNB()

    def estimatorTest(self, testData=None):
        X_test = testData.iloc[:, :-1]
        y_test = testData.iloc[:, -1]

        from sklearn.metrics import precision_recall_fscore_support
        y_pred = self.finalEstimator.predict(X_test)
        y_pred_1 = y_pred.copy()
        y_pred_1[np.where(y_pred_1 == 0)[0]] = -1
        precision, recall, f1_score, _ = precision_recall_fscore_support(y_test, y_pred_1)

        self.f1_score = f1_score[1]
        self.precision = precision[1]
        self.recall = recall[1]
        print("F1 score:", f1_score[1])
        print("Precision:", precision[1])
        print("Recall:", recall[1])


    def filterIsoformsByScore(self, filterScore=0.0, outFile="validIsoform.lst", lqAnnoIso=None):
        hqNovelIsoRes = self.predResults.loc[(self.predResults.label==0) & (self.predResults.pu_score>=filterScore), ]
        annoIsoRes = self.predResults.loc[self.predResults.label==1,]
        if lqAnnoIso is not None:
            lqAnnoIsoRes = lqAnnoIso.loc[:, "label"]
            lqAnnoIsoRes["pu_score"] = 1
            outResults = pd.concat([annoIsoRes, lqAnnoIsoRes, hqNovelIsoRes])
        else:
            outResults = pd.concat([annoIsoRes, hqNovelIsoRes])
        outResults.to_csv(outFile, sep="\t", index=True, index_label="isoform")


class GetFeatures(object):
    def __init__(self, inputBed=None, genomeFa=None, inputFa=None):
        self.inputBed = inputBed
        self.genomeFa = genomeFa
        self.inputFa = inputFa
        self.features = None
        self.inputBedObj = None

    def addIsoSupport(self):
        myDict = {}
        with open(self.inputBed) as f:
            for line in f.readlines():
                infoList = line.strip("\n").split("\t")
                isoName = infoList[3]
                flCount = int(infoList[-4])
                ratioIsoToGene = float(infoList[-2])
                myDict[isoName] = {"flCount": flCount, "ratioIsoToGene": ratioIsoToGene}

        newFeatures = pd.DataFrame.from_dict(myDict, orient="index")
        return pd.concat([self.features, newFeatures], axis=1)

    # def addIsoNgsSupp(self, isoNgsExp=None):
    #     self.features.update()
    #     return self.features

    def addJuncIndelInfo(self, samFile=None, junctionFile=None, refBed=None):
        juncDict = BedFile(junctionFile, type="bed12").getAllJuncDict()
        annoJuncDict = BedFile(refBed, type="bed12+").getAllJuncDict()
        sam = pysam.AlignmentFile(samFile, "r")

        indelsTotal = {}
        juncWithIndel = {}
        readCount = 0
        for read in sam.fetch():
            if read.is_unmapped: continue
            readCount += 1
            cigar = read.cigar
            pos = read.pos + 1
            spliceSites = []

            for cigarType, cigarLength in cigar:
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1
                    if cigarType == 3:
                        # juncName = "{}:{}-{}".format(read.reference_id, pos, posEnd)
                        # if juncName in juncDict:
                        #     spliceSites.append([pos, posEnd])
                        spliceSites.append([pos, posEnd])
                    pos = posEnd + 1

            pos = read.pos + 1
            posEnd = pos
            for cigarType, cigarLength in cigar:
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1

                if cigarType == 1:
                    indelStartPos = pos
                    indelEndPos = pos
                    indelLength = cigarLength
                    indelName = "{}:{}:{}+{}".format(read.reference_name, "insertion", indelStartPos, indelLength)

                    if indelName not in indelsTotal:
                        indelsTotal[indelName] = 1
                    else:
                        indelsTotal[indelName] += 1

                    for i in spliceSites:
                        juncName = "{}:{}-{}".format(read.reference_name, i[0], i[1])
                        if any(abs(indelStartPos-e) < 10 for e in i) or any(abs(indelEndPos-e) < 10 for e in i):
                            if juncName not in juncWithIndel:
                                juncWithIndel[juncName] = [indelName]
                            else:
                                juncWithIndel[juncName].append(indelName)

                if cigarType == 2:
                    indelStartPos = pos
                    indelEndPos = posEnd
                    indelLength = cigarLength
                    indelName = "{}:{}:{}-{}".format(read.reference_name, "deletion", indelStartPos, indelLength)

                    if indelName not in indelsTotal:
                        indelsTotal[indelName] = 1
                    else:
                        indelsTotal[indelName] += 1

                    for i in spliceSites:
                        juncName = "{}:{}-{}".format(read.reference_name, i[0], i[1])
                        if any(abs(indelStartPos-e) < 10 for e in i) or any(abs(indelEndPos-e) < 10 for e in i):
                            if juncName not in juncWithIndel:
                                juncWithIndel[juncName] = [indelName]
                            else:
                                juncWithIndel[juncName].append(indelName)

                if cigarType not in [1, 4, 5, 7, 8]:
                    pos = posEnd + 1
        sam.close()

        myDict = {}
        for isoName in self.inputBedObj.reads:
            if len(self.inputBedObj.reads[isoName].exons) > 1:
                myDict[isoName] = {
                    "nIndelsAroundJunc": 0,
                    "nJuncsWithIndels": 0,
                    "indelNearJunc": 0,
                    "ratioMinJuncCovToAllCov": 0,
                    "sdJuncCov": 0,
                    "minJuncRPKM": 0,
                    "withNovelJunc": False,
                    "minNovelJuncRPKM": 0
                }
                isoObj = self.inputBedObj.reads[isoName]
                isoLength = sum(isoObj.blockSizes)
                tmpDict = {}
                for junc in isoObj.introns:
                    juncName = "{}:{}-{}".format(isoObj.chrom, junc[0]+1, junc[1])
                    if juncName in juncDict:
                        juncCov = int(juncDict[juncName].split("\t")[4])
                        indelsNearJunc = list(set(juncWithIndel[juncName])) if juncName in juncWithIndel else []
                        validIndelsNearJunc = [indel for indel in indelsNearJunc if float(indelsTotal[indel]/juncCov) >= 0.5]
                        juncRPKM = juncCov / (isoLength/1000.0 * readCount/1000000.0)
                        juncAnno = True if juncName in annoJuncDict else False
                        tmpDict.update({juncName: {"juncCov": juncCov, "validIndels": validIndelsNearJunc,
                                                   "juncRPKM": juncRPKM, "juncAnno": juncAnno}})

                isoJuncCovs = [tmpDict[i]["juncCov"] for i in tmpDict]
                novelJuncRPKM = [tmpDict[i]["juncRPKM"] for i in tmpDict if tmpDict[i]["juncAnno"] == False]
                isoValidIndels = list(set(itertools.chain.from_iterable([tmpDict[i]["validIndels"] for i in tmpDict])))
                if isoValidIndels:
                    myDict[isoName]["indelNearJunc"] = 1
                    myDict[isoName]["nIndelsAroundJunc"] = len(isoValidIndels)
                    myDict[isoName]["nJuncsWithIndels"] = len([i for i in tmpDict if len(tmpDict[i]["validIndels"])])
                if isoJuncCovs:
                    myDict[isoName]["ratioMinJuncCovToAllCov"] = float(min(isoJuncCovs)) / sum(isoJuncCovs)
                    myDict[isoName]["sdJuncCov"] = np.std(isoJuncCovs)
                    myDict[isoName]["minJuncRPKM"] = min([tmpDict[i]["juncRPKM"] for i in tmpDict])
                if novelJuncRPKM:
                    myDict[isoName]["withNovelJunc"] = True
                    myDict[isoName]["minNovelJuncRPKM"] = min(novelJuncRPKM)


        newFeatures = pd.DataFrame.from_dict(myDict, orient="index")
        return pd.concat([self.features, newFeatures], axis=1)


    def getInputFa(self):
        tmpInputFa = os.path.join(os.getcwd(), "tmpInput.fa")
        cmd = '''cut -f 1-12 {} | 
            bedtools getfasta -fi {} -bed - -name -split -s | 
            seqkit replace -w 0 -p "(.*?):(.*)" -r '$1' > {}
        '''.format(self.inputBed, self.genomeFa, tmpInputFa)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        return tmpInputFa


    def getFeatures(self, refBed):
        from CPC2 import FindCDS
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        inputFaDict = SeqIO.to_dict(SeqIO.parse(self.inputFa, "fasta"))

        self.inputBedObj = BedFile(self.inputBed, type="bed12+")
        refBedObj = BedFile(refBed, type="bed12+")

        allAnnoJuncDict = refBedObj.getAllJuncDict()
        allAnnoExonDict = refBedObj.getAllExonDict()
        isoJuncChainDict = self.inputBedObj.getJuncChainDict()
        isos = self.inputBedObj.reads
        canonicalJuncFreq = self.getCanonicalJuncFreq(isoJuncChainDict, self.genomeFa)

        myDict = {}
        novelJuncList = []
        for isoName in isos:
            myDict[isoName] = {}
            juncChain = "{}:{}".format(isos[isoName].chrom, isos[isoName].juncChain)
            strand = isos[isoName].strand

            isoLength = sum(isos[isoName].blockSizes)
            exonNum = len(isos[isoName].blockSizes)
            annotation = isos[isoName].otherList[-1]
            orfSeq, startPos, orfStrand, orfFullness = FindCDS(inputFaDict[isoName].seq.upper()).longest_orf(strand)
            orfLen = len(orfSeq)
            canJuncRatio = canonicalJuncFreq[juncChain]
            gc = GC(inputFaDict[isoName].seq)
            myDict[isoName].update({
                "isoLength": isoLength,
                "exonNum": exonNum,
                "annotation": annotation,
                "orfLength": orfLen,
                "canJuncRatio": canJuncRatio,
                "GC": gc,
                "bite": False
            })

            for junc in isos[isoName].introns:
                juncStart, juncEnd = junc
                juncName = "{}:{}-{}".format(isos[isoName].chrom, juncStart+1, juncEnd)
                if juncName not in allAnnoJuncDict:
                    juncInfo = "\t".join(map(str, [isos[isoName].chrom, juncStart+1, juncEnd, isoName, ".", strand]))
                    novelJuncList.append(juncInfo)

        novelJuncObj = pybedtools.BedTool("\n".join(novelJuncList), from_string=True)
        annoExonObj = pybedtools.BedTool("\n".join(allAnnoExonDict.values()), from_string=True)
        intersectRes = novelJuncObj.intersect(annoExonObj, wa=True, wb=True, s=True)
        biteIsoforms = {}
        for i in intersectRes:
            infoList = str(i).strip("\n").split("\t")
            if infoList[3] not in biteIsoforms:
                biteIsoforms[infoList[3]] = ""

        for i in biteIsoforms:
            myDict[i]["bite"] = True

        return pd.DataFrame.from_dict(myDict, orient="index")


    def getCanonicalJuncFreq(self, isoJuncChainDict, genomeFasta):
        dinucleotideBedList = []
        for juncChain in isoJuncChainDict:
            for j in juncChain.split(":")[1].split(";"):
                chrom = isoJuncChainDict[juncChain].chrom
                strand = isoJuncChainDict[juncChain].strand
                juncStart, juncEnd = map(int, j.split("-"))
                juncName = "{}:{}-{}".format(chrom, juncStart+1, juncEnd)
                leftDinucleotidePos = "\t".join(map(str, [chrom, juncStart, juncStart + 2, ":".join([juncName, "left"]), ".", strand]))
                rightDinucleotidePos = "\t".join(map(str, [chrom, juncEnd - 2, juncEnd, ":".join([juncName, "right"]), ".", strand]))
                dinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])

        dinucleotideBedObj = pybedtools.BedTool("\n".join(dinucleotideBedList), from_string=True)
        dinucleotideBedRes = dinucleotideBedObj.sequence(genomeFasta, name=True, tab=True, s=True)
        juncCanonicalDict = {}
        for i in str(open(dinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            juncName, dinucleotideType = ":".join(infoList[0].split(":")[:2]), infoList[0].split(":")[2]

            if juncName not in juncCanonicalDict:
                juncCanonicalDict.update({juncName: {dinucleotideType: infoList[1]}})
            else:
                juncCanonicalDict[juncName].update({dinucleotideType: infoList[1]})

        canonicalJuncFreq = {}
        for juncChain in isoJuncChainDict:
            juncCount = len(juncChain.split(":")[1].split(";"))
            strand = isoJuncChainDict[juncChain].strand
            canonicalCount = 0
            for j in juncChain.split(":")[1].split(";"):
                juncStart = int(j.split("-")[0])+1
                juncEnd = int(j.split("-")[1])
                juncName = "{}:{}-{}".format(isoJuncChainDict[juncChain].chrom, juncStart, juncEnd)
                if strand == "+":
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["left"],
                                                 juncCanonicalDict[juncName]["right"])
                else:
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["right"],
                                                 juncCanonicalDict[juncName]["left"])
                if spliceMotif in ["GT-AG", "GC-AG", "AT-AC"]:
                    canonicalCount += 1
            canonicalJuncFreq[juncChain] = float(canonicalCount) / juncCount
        return canonicalJuncFreq


    def predCodingP(self, refBed):
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        self.features = self.features if not self.features.empty else self.getFeatures(refBed)

        from CPC2 import calculate_potential
        calculate_potential(self.inputFa, "+", 0, "cpc2out")
        codingPotential = pd.read_csv("cpc2out.txt", sep="\t", skiprows=1, header=None,
                                      usecols=[0, 2, 3, 4, 5, 6, 7], index_col=0)
        codingPotential.set_axis(["pepLength", "FickettScore", "pI", "orfIntegrity", "codingP", "codingLabel"], axis=1, inplace=True)

        return pd.concat([self.features, codingPotential], axis=1)


def getInvolvedIsos(asFile, asType="IR"):
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


def getInputBed(baseDir=None, refBed=None, useAsFile=False):
    isoformFile = os.path.join(baseDir, "refine", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    isoformBed = BedFile(isoformFile, type="bed12+").reads

    isoDict = {}
    if useAsFile:
        irFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "IR.confident.bed6+")
        seFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "SE.confident.bed12+")
        a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A3SS.confident.bed6+")
        a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "LR", "A5SS.confident.bed6+")

        isoDict.update(getInvolvedIsos(irFile, asType="IR"))
        isoDict.update(getInvolvedIsos(seFile, asType="SE"))
        isoDict.update(getInvolvedIsos(a5ssFile, asType="A5SS"))
        isoDict.update(getInvolvedIsos(a3ssFile, asType="A3SS"))
    else:
        isoDict = dict.fromkeys(isoformBed.keys(), "")

    # annoJuncDict, annoSingleExonList = BedFile(refBed, type="bed12+").getJuncChainDict()
    annoJuncDict = BedFile(refBed, type="bed12+").getJuncChainDict()
    isoform2reads = getDictFromFile(collapsedTrans2reads, sep="\t", inlineSep=",", valueCol=2)

    juncDict = {}
    gene2isoDict = {}
    isoSingleExonList = []
    for iso in isoDict:
        gene = isoformBed[iso].otherList[0]
        if gene not in gene2isoDict:
            gene2isoDict[gene] = {"isos": [iso], "count": len(isoform2reads[iso])}
        else:
            gene2isoDict[gene]["isos"].append(iso)
            gene2isoDict[gene]["count"] += len(isoform2reads[iso])

        if isoformBed[iso].juncChain == "":
            isoSingleExonList.append("{}\t{}\t{}\t{}\t{}\t{}".format(isoformBed[iso].chrom, isoformBed[iso].chromStart,
                                                                     isoformBed[iso].chromEnd, iso, ".",
                                                                     isoformBed[iso].strand))
            continue
        juncInfo = "{}:{}".format(isoformBed[iso].chrom, isoformBed[iso].juncChain)
        if juncInfo not in juncDict:
            isoLength = getBlockLength(isoformBed[iso].exons)
            juncDict[juncInfo] = {"iso": [iso], "gene": [isoformBed[iso].otherList[0]], "longest": [iso, isoLength]}
        else:
            juncDict[juncInfo]["iso"].append(iso)
            juncDict[juncInfo]["gene"].append(isoformBed[iso].otherList[0])
            if getBlockLength(isoformBed[iso].exons) > juncDict[juncInfo]["longest"][1]:
                juncDict[juncInfo]["longest"][0] = iso
                juncDict[juncInfo]["longest"][1] = getBlockLength(isoformBed[iso].exons)

    tmpOut = open(os.path.join(os.getcwd(), "inputBed.tmp"), "w")
    for junc in juncDict:
        longestIso = juncDict[junc]["longest"][0]
        gene = isoformBed[longestIso].otherList[0]
        isoSupport = sum([len(isoform2reads[x]) for x in juncDict[junc]["iso"]])
        geneSupport = gene2isoDict[gene]["count"]
        ratio = float(isoSupport) / geneSupport
        annotation = "annotated" if junc in annoJuncDict else "novel"
        print >> tmpOut, "\t".join(map(str, [str(isoformBed[longestIso]), ",".join(juncDict[junc]["iso"]), isoSupport, geneSupport, ratio, annotation]))
    tmpOut.close()
    return os.path.join(os.getcwd(), "inputBed.tmp")


def processInputFile(inputBed, genomeFa, refBed, junctionFile=None, samFile=None):
    gf = GetFeatures(inputBed=inputBed, genomeFa=genomeFa)
    gf.features = gf.getFeatures(refBed)
    gf.features = gf.predCodingP(refBed)
    gf.features = gf.addIsoSupport()

    # if isoNgsExpFile:
    #     features = gf.addIsoNgsSupp(isoNgsExpFile)

    if samFile:
        gf.features = gf.addJuncIndelInfo(samFile=samFile, junctionFile=junctionFile, refBed=refBed)

    gf.features = gf.features.round(3)
    colOrder = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "annotation", "GC", "orfLength", "orfIntegrity",
                "pepLength", "FickettScore", "pI", "codingP", "codingLabel", "canJuncRatio", "sdJuncCov",
                "withNovelJunc", "minNovelJuncRPKM", "bite", "minJuncRPKM", "nIndelsAroundJunc",
                "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc"]
    gf.features = gf.features.loc[:, colOrder]
    gf.features.to_csv("isoFeatures.txt", sep="\t", index=True, index_label="isoform")

    return gf.features


def selectBestModel(featureData, drawAUC=False, selectPUscore=False):
    from sklearn.model_selection import KFold

    usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
                    "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
                    "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]

    annoIsoData = featureData.loc[featureData.annotation == "annotated",]
    novelIsoData = featureData.loc[featureData.annotation == "novel",]

    leastExpIsos = novelIsoData.loc[(novelIsoData.ratioIsoToGene < 0.05),]
    unreliableJuncIsos = novelIsoData.loc[(novelIsoData.withNovelJunc == True) &
                                          (novelIsoData.minNovelJuncRPKM < 0.05),]
    validNegIndex = list(set(leastExpIsos.index) | set(unreliableJuncIsos.index))

    posIsoData = annoIsoData.loc[(annoIsoData.flCount >= 2) & (annoIsoData.minJuncRPKM >= 0.05),]
    posIsoDataInner, posIsoDataOuter = train_test_split(posIsoData, test_size=0.2, random_state=0)
    validNegIndexInner, validNegOuter = train_test_split(novelIsoData.loc[validNegIndex,], test_size=0.2,
                                                         random_state=0)

    outerPosData = posIsoDataInner.copy()
    outerPosData["label"] = 1
    posIsoDataOuterCopy = posIsoDataOuter.copy()
    posIsoDataOuterCopy["label"] = 1
    validNegOuterCopy = validNegOuter.copy()
    validNegOuterCopy["label"] = 0

    outerUnlabeledData = pd.concat([posIsoDataOuterCopy, validNegOuterCopy])
    outerUnlabeledData["label"] = 0

    outerPosData = outerPosData.loc[:, usedFeatures]
    outerPosData = outerPosData.sample(frac=1)
    outerUnlabeledData = outerUnlabeledData.loc[:, usedFeatures]
    outerUnlabeledData = outerUnlabeledData.sample(frac=1)

    outerPosData.replace({False: 0, True: 1}, inplace=True)
    outerUnlabeledData.replace({False: 0, True: 1}, inplace=True)
    outerTrainingData = pd.concat([outerPosData, outerUnlabeledData])

    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    colors = ["#1f497d", "#f79646", "#9bbb59", "#7f7f7f", "#8064a2"]
    models = ["RF", "GB", "DT", "SVM", "NB"]
    names = ["RF", "GB", "DT", "SVM", "NB"]

    aucScore = dict.fromkeys(models, {})
    puScoreTholdDict = {}
    for scheme in ["bagging"]:
        if drawAUC:
            fig1 = plt.figure(figsize=[6, 6])
            ax1 = fig1.add_subplot(111, aspect='equal')
        for model, color, name in zip(models, colors, names):
            tprs = []
            # aucs = []
            mean_fpr = np.linspace(0, 1, 100)
            # i = 1
            for train_i, test_i in kf.split(posIsoDataInner):
                train_index = posIsoDataInner.index[train_i]
                test_index = posIsoDataInner.index[test_i]
                posIsoDataTraining = posIsoDataInner.loc[train_index,]
                posIsoDataTest = posIsoDataInner.loc[test_index,]
                negIsoDataTest = novelIsoData.loc[validNegIndexInner.index,]

                posIsoDataTraining["label"] = 1
                posIsoDataTest["label"] = 0
                negIsoDataTest["label"] = 0

                trainingData = pd.concat([posIsoDataTraining, posIsoDataTest, negIsoDataTest])
                trainingData = trainingData.loc[:, usedFeatures]
                trainingData = trainingData.sample(frac=1)

                posIsoDataTest["label"] = 1
                testData = pd.concat([posIsoDataTest, negIsoDataTest])
                testData = testData.loc[:, usedFeatures]
                testData = testData.sample(frac=1)

                trainingData.replace({False: 0, True: 1}, inplace=True)
                testData.replace({False: 0, True: 1}, inplace=True)

                sme = SchemeModelEval(trainingData=trainingData, model=model, scheme=scheme)
                sme.eval()
                prediction = sme.finalEstimator.predict_proba(testData.iloc[:, :-1])
                pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=testData.index)
                results = pd.DataFrame({
                    "true_label": testData.label,
                    "train_label": trainingData.loc[trainingData.label == 0,].label,
                    "pu_score": pu_score.pu_score
                }, columns=["true_label", "train_label", "pu_score"])

                if selectPUscore:
                    thold = np.quantile(results.loc[results.true_label==0].pu_score, 0.95)
                    if model not in puScoreTholdDict:
                        puScoreTholdDict[model] = [thold]
                    else:
                        puScoreTholdDict[model].append(thold)

                fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
                tprs.append(interp(mean_fpr, fpr, tpr))
                # roc_auc = auc(fpr, tpr)
                # aucs.append(roc_auc)
                # i = i + 1

            mean_tpr = np.mean(tprs, axis=0)
            mean_auc = auc(mean_fpr, mean_tpr)
            aucScore[model].update({"trainingAUC": mean_auc})
            if drawAUC:
                plt.plot(mean_fpr, mean_tpr, lw=2, linestyle='--', color=color,
                         label='%s (CV AUC = %0.3f )' % (name, mean_auc))

        ######################

        for model, color, name in zip(models, colors, names):
            sme = SchemeModelEval(trainingData=outerTrainingData, model=model, scheme=scheme)
            sme.eval()
            prediction = sme.finalEstimator.predict_proba(outerUnlabeledData.iloc[:, :-1])
            pu_score = pd.DataFrame({"pu_score": prediction[:, 1]}, index=outerUnlabeledData.index)
            results = pd.DataFrame({
                "true_label": pd.concat([posIsoDataOuterCopy, validNegOuterCopy]).label,
                "train_label": outerUnlabeledData.label,
                "pu_score": pu_score.pu_score
            }, columns=["true_label", "train_label", "pu_score"])

            fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
            roc_auc = auc(fpr, tpr)
            aucScore[model].update({"testAUC": roc_auc})
            if drawAUC:
                plt.plot(fpr, tpr, lw=2, color=color, label='%s (Test AUC = %0.3f )' % (name, roc_auc))

        # plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey')
        if drawAUC:
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('ROC_AUC')
            plt.legend(loc="lower right")

            plt.savefig('PU.{}.ROC_AUC.pdf'.format(scheme))
            plt.close()

    tmpScore = 0
    tmpModel = None
    for i in aucScore:
        if np.mean(aucScore[i].values()) > tmpScore:
            tmpModel = i
            tmpScore = np.mean(aucScore[i].values())
    return tmpModel, puScoreTholdDict


def iso_pu1(dataObj=None, dirSpec=None, refParams=None, hqIsoParams=None, samFile=None, junctionFile=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Start finding high quality novel isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    prevDir = os.getcwd()
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    hqIsoformDir = os.path.join(baseDir, "hqIsoforms")
    resolveDir(hqIsoformDir)

    genomeFa = refParams.ref_genome
    refBed = refParams.ref_bed

    inputBed = None

    if samFile == None and os.path.exists(os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")):
        samFile = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")
    if junctionFile == None and os.path.exists(os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")):
        junctionFile = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")

    if hqIsoParams.feature_file and validateFile(hqIsoParams.feature_file):
        featureData = hqIsoParams.feature_file
    else:
        inputBed = inputBed if inputBed else getInputBed(baseDir=baseDir, refBed=refBed)
        featureData = processInputFile(inputBed, genomeFa, refBed, junctionFile=junctionFile, samFile=samFile)
        # featureData = pd.read_csv("isoFeatures.txt", sep="\t", index_col=0)

    annoIsoData = featureData.loc[featureData.annotation == "annotated",]
    novelIsoData = featureData.loc[featureData.annotation == "novel",]
    posIsoData = annoIsoData.loc[(annoIsoData.flCount >= int(hqIsoParams.pos_fl_coverage)) &
                                 (annoIsoData.minJuncRPKM >= float(hqIsoParams.pos_min_junc_rpkm)),]
    lqAnnoIsoData = annoIsoData.loc[~annoIsoData.index.isin(posIsoData.index),]

    posIsoDataCopy = posIsoData.copy()
    posIsoDataCopy["label"] = 1
    lqAnnoIsoDataCopy = lqAnnoIsoData.copy()
    lqAnnoIsoDataCopy["label"] = 1
    novelIsoDataCopy = novelIsoData.copy()
    novelIsoDataCopy["label"] = 0

    usedFeatures = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "GC", "orfLength", "orfIntegrity",
                    "pepLength", "FickettScore", "pI", "codingP", "canJuncRatio", "sdJuncCov", "minJuncRPKM",
                    "nIndelsAroundJunc", "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc", "label"]
    dataToClassify = pd.concat([posIsoDataCopy, novelIsoDataCopy])
    dataToClassify = dataToClassify.loc[:, usedFeatures]
    dataToClassify = dataToClassify.sample(frac=1)
    dataToClassify.replace({False: 0, True: 1}, inplace=True)

    if hqIsoParams.select_best_model:
        if hqIsoParams.auto_filter_score:
            bestModel, puScoreTholdDict = selectBestModel(featureData, drawAUC=hqIsoParams.draw_auc, selectPUscore=True)
            filter_score = np.mean(puScoreTholdDict[bestModel])
        else:
            bestModel, puScoreTholdDict = selectBestModel(featureData, drawAUC=hqIsoParams.draw_auc, selectPUscore=False)
            filter_score = float(hqIsoParams.filter_score)
    else:
        bestModel = "GB"
        filter_score = float(hqIsoParams.filter_score)
    sme = SchemeModelEval(trainingData=dataToClassify, model=bestModel, scheme="bagging")
    # sme = SchemeModelEval(trainingData=dataToClassify, model="GB", scheme="bagging")
    sme.eval()
    sme.predResults.to_csv("pu_score.txt", sep="\t")
    sme.filterIsoformsByScore(filterScore=filter_score, outFile="validIsoforms.lst", lqAnnoIso=lqAnnoIsoDataCopy)

    cmd = '''(grep 'novel' isoFeatures.txt | awk '{if($3>=15 && $4>=0.4){print $1}}'; cut -f 1 validIsoforms.lst) | %s/filter.pl -o - %s -2 4 -m i > hq.collapsed.bed12+''' % (utilDir, inputBed)
    subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    print getCurrentTime() + " Finding high quality novel isoforms for project {} sample {} done!".format(projectName, sampleName)

