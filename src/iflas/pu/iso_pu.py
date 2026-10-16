#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iso_pu.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2022-04-23 14:22:38
Last modified: 2022-04-23 14:22:38
'''

import os, pybedtools, pysam, itertools, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from collections import Counter
from numpy import interp  # scipy.interp removed in SciPy 1.x (was an alias of numpy.interp)
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

from .bagging import BaggingClassifierPU
from ..io.records import BedFile
from ..io.files import getDictFromFile
from ..intervalops import getBlockLength
from ..runtime import ToolRunner, quote as q

warnings.filterwarnings("ignore")


class SchemeModelEval(object):
    def __init__(self, trainingData=None, model=None, scheme=None, featureData=None,
                 bagging_n_jobs=-1, random_state=0):
        self.trainingData = trainingData
        self.random_state = random_state
        self.model = model
        self.scheme = scheme

        self.finalEstimator = None
        self.predResults = None
        self.f1_score = None
        self.precision = None
        self.recall = None

        self.featureData = featureData
        self.bagging_n_jobs = bagging_n_jobs

    def eval(self):
        X_train = self.trainingData.iloc[:, :-1]
        y_train = self.trainingData.iloc[:, -1]

        classifier = self.getClassifier()
        if self.scheme == "bagging":
            X = X_train
            y = y_train
            y = pd.Series(y)

            labels = Counter(y)
            sorted_labels = sorted(list(labels.items()), key=lambda x: x[1], reverse=True)
            max_samples = sorted_labels[0][1] if sorted_labels[0][0] == 0 else sorted_labels[1][1]

            cf = BaggingClassifierPU(
                classifier,
                n_estimators=100,
                max_samples=max_samples,
                n_jobs=self.bagging_n_jobs,
                random_state=self.random_state
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


    def filterIsoformsByScore(self, filterScore=0.0, outFile="validIsoform.lst", lqAnnoIso=None):
        hqNovelIsoRes = self.predResults.loc[(self.predResults.label==0) & (self.predResults.pu_score>=filterScore), ]
        annoIsoRes = self.predResults.loc[self.predResults.label==1,]
        if lqAnnoIso is not None:
            # a one-column frame plus a score column (legacy selected the Series and so added a junk row "pu_score": BUG-179)
            lqAnnoIsoRes = lqAnnoIso.loc[:, ["label"]].assign(pu_score=1)
            outResults = pd.concat([annoIsoRes, lqAnnoIsoRes, hqNovelIsoRes])
        else:
            outResults = pd.concat([annoIsoRes, hqNovelIsoRes])
        outResults.to_csv(outFile, sep="\t", index=True, index_label="isoform")


class GetFeatures(object):
    def __init__(self, inputBed=None, genomeFa=None, inputFa=None, runner=None):
        self.inputBed = inputBed
        self.genomeFa = genomeFa
        self.inputFa = inputFa
        self.runner = runner
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

    @staticmethod
    def scanIndelsNearJunctions(samFile):
        """NGS reads -> (indel counts by name, indels within 10 bp of each read's splice sites, mapped read count)."""
        indelsTotal, juncWithIndel, readCount = {}, {}, 0

        def record(read, label, start, end, length, spliceSites):
            indelName = "{}:{}:{}{}{}".format(read.reference_name, label, start, "+" if label == "insertion" else "-", length)
            indelsTotal[indelName] = indelsTotal.get(indelName, 0) + 1
            for site in spliceSites:
                if any(abs(start - e) < 10 for e in site) or any(abs(end - e) < 10 for e in site):
                    juncWithIndel.setdefault("{}:{}-{}".format(read.reference_name, site[0], site[1]), []).append(indelName)

        sam = pysam.AlignmentFile(samFile, "r")
        for read in sam.fetch():
            if read.is_unmapped:
                continue
            readCount += 1
            cigar = read.cigar
            pos, spliceSites = read.pos + 1, []
            for cigarType, cigarLength in cigar:                      # reference-consuming operations only
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1
                    if cigarType == 3:
                        spliceSites.append([pos, posEnd])
                    pos = posEnd + 1
            pos = read.pos + 1
            posEnd = pos
            for cigarType, cigarLength in cigar:
                if cigarType not in [1, 4, 5, 7, 8]:
                    posEnd = pos + cigarLength - 1
                if cigarType == 1:
                    record(read, "insertion", pos, pos, cigarLength, spliceSites)
                if cigarType == 2:
                    record(read, "deletion", pos, posEnd, cigarLength, spliceSites)
                if cigarType not in [1, 4, 5, 7, 8]:
                    pos = posEnd + 1
        sam.close()
        return indelsTotal, juncWithIndel, readCount

    def addJuncIndelInfo(self, samFile=None, junctionFile=None, refBed=None):
        juncDict = BedFile(junctionFile, type="bed12").getAllJuncDict()
        annoJuncDict = BedFile(refBed, type="bed12+").getAllJuncDict()
        indelsTotal, juncWithIndel, readCount = self.scanIndelsNearJunctions(samFile)

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
                        validIndelsNearJunc = [indel for indel in indelsNearJunc if float(indelsTotal[indel]) / juncCov >= 0.5]
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
        '''.format(q(self.inputBed), q(self.genomeFa), q(tmpInputFa))
        runner = self.runner or ToolRunner()
        runner.run(["bash", "-c", cmd])
        return tmpInputFa


    def getFeatures(self, refBed):
        from .cpc2 import FindCDS
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        inputFaDict = SeqIO.to_dict(SeqIO.parse(self.inputFa, "fasta"))

        self.inputBedObj = BedFile(self.inputBed, type="bed12+")
        refBedObj = BedFile(refBed, type="bed12+")

        allAnnoJuncDict = refBedObj.getAllJuncDict()
        allAnnoExonDict = refBedObj.getAllExonDict()
        isoJuncChainDict = self.inputBedObj.getJuncChainDict(withStrand=True)     # BUG-196: chains are per strand
        isos = self.inputBedObj.reads
        canonicalJuncFreq = self.getCanonicalJuncFreq(isoJuncChainDict, self.genomeFa)

        myDict = {}
        novelJuncList = []
        for isoName in isos:
            myDict[isoName] = {}
            strand = isos[isoName].strand
            juncChain = "{}:{}:{}".format(isos[isoName].chrom, isos[isoName].juncChain, strand)

            isoLength = sum(isos[isoName].blockSizes)
            exonNum = len(isos[isoName].blockSizes)
            annotation = isos[isoName].otherList[-1]
            orfSeq, startPos, orfStrand, orfFullness = FindCDS(inputFaDict[isoName].seq.upper()).longest_orf(strand)
            orfLen = len(orfSeq)
            canJuncRatio = canonicalJuncFreq[juncChain]
            gc = gc_fraction(inputFaDict[isoName].seq) * 100  # gc_fraction is 0-1; *100 to match old Bio.SeqUtils.GC (BUG-106)
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
        annoExonObj = pybedtools.BedTool("\n".join(list(allAnnoExonDict.values())), from_string=True)
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
                juncStart, juncEnd = list(map(int, j.split("-")))
                juncName = "{}:{}-{}{}".format(chrom, juncStart+1, juncEnd, strand)
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
                juncName = "{}:{}-{}{}".format(isoJuncChainDict[juncChain].chrom, juncStart, juncEnd, strand)
                if strand == "+":
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["left"],
                                                 juncCanonicalDict[juncName]["right"]).upper()
                else:
                    spliceMotif = "{}-{}".format(juncCanonicalDict[juncName]["right"],
                                                 juncCanonicalDict[juncName]["left"]).upper()
                if spliceMotif in ["GT-AG", "GC-AG", "AT-AC"]:
                    canonicalCount += 1
            canonicalJuncFreq[juncChain] = float(canonicalCount) / juncCount
        return canonicalJuncFreq


    def predCodingP(self, refBed):
        self.inputFa = self.inputFa if self.inputFa else self.getInputFa()
        self.features = self.features if not self.features.empty else self.getFeatures(refBed)

        from .cpc2 import calculate_potential
        calculate_potential(self.inputFa, "+", 0, "cpc2out", runner=self.runner)
        codingPotential = pd.read_csv("cpc2out.txt", sep="\t", skiprows=1, header=None,
                                      usecols=[0, 2, 3, 4, 5, 6, 7], index_col=0)
        # BUG-113: pandas>=2 removed set_axis(inplace=); rebind instead (semantically identical)
        codingPotential = codingPotential.set_axis(
            ["pepLength", "FickettScore", "pI", "orfIntegrity", "codingP", "codingLabel"], axis=1)

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
        isoDict = dict.fromkeys(list(isoformBed.keys()), "")

    # annoJuncDict, annoSingleExonList = BedFile(refBed, type="bed12+").getJuncChainDict()
    annoJuncDict = BedFile(refBed, type="bed12+").getJuncChainDict(withStrand=True)
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
        # BUG-196: the same intron chain on the other strand is another transcript (legacy merged them and labelled
        # them by the annotation of either strand)
        juncInfo = "{}:{}:{}".format(isoformBed[iso].chrom, isoformBed[iso].juncChain, isoformBed[iso].strand)
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
        print("\t".join(map(str, [str(isoformBed[longestIso]), ",".join(juncDict[junc]["iso"]), isoSupport, geneSupport, ratio, annotation])), file=tmpOut)
    tmpOut.close()
    return os.path.join(os.getcwd(), "inputBed.tmp")


def processInputFile(inputBed, genomeFa, refBed, junctionFile=None, samFile=None, runner=None):
    if not (samFile and junctionFile):       # BUG-191: was a KeyError on the missing NGS feature columns below
        raise ValueError("pu_filter needs the sample's NGS alignments and junctions (mapping/rna-seq/reassembly/tmp.bam "
                         "and junctions.bed): the PU features include NGS junction support, as in the published iFLAS. "
                         "Give the sample ngs_left_reads/ngs_right_reads, or a feature table (pu_filter.feature_file)")
    gf = GetFeatures(inputBed=inputBed, genomeFa=genomeFa, runner=runner)
    gf.features = gf.getFeatures(refBed)
    gf.features = gf.predCodingP(refBed)
    gf.features = gf.addIsoSupport()

    gf.features = gf.addJuncIndelInfo(samFile=samFile, junctionFile=junctionFile, refBed=refBed)

    gf.features = gf.features.round(3)
    colOrder = ["isoLength", "flCount", "ratioIsoToGene", "exonNum", "annotation", "GC", "orfLength", "orfIntegrity",
                "pepLength", "FickettScore", "pI", "codingP", "codingLabel", "canJuncRatio", "sdJuncCov",
                "withNovelJunc", "minNovelJuncRPKM", "bite", "minJuncRPKM", "nIndelsAroundJunc",
                "ratioMinJuncCovToAllCov", "nJuncsWithIndels", "indelNearJunc"]
    gf.features = gf.features.loc[:, colOrder]
    gf.features.to_csv("isoFeatures.txt", sep="\t", index=True, index_label="isoform")

    return gf.features


def selectBestModel(featureData, drawAUC=False, selectPUscore=False, random_state=0):
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
    posIsoDataInner, posIsoDataOuter = train_test_split(posIsoData, test_size=0.2, random_state=random_state)
    validNegIndexInner, validNegOuter = train_test_split(novelIsoData.loc[validNegIndex,], test_size=0.2,
                                                         random_state=random_state)

    outerPosData = posIsoDataInner.copy()
    outerPosData["label"] = 1
    posIsoDataOuterCopy = posIsoDataOuter.copy()
    posIsoDataOuterCopy["label"] = 1
    validNegOuterCopy = validNegOuter.copy()
    validNegOuterCopy["label"] = 0

    outerUnlabeledData = pd.concat([posIsoDataOuterCopy, validNegOuterCopy])
    outerUnlabeledData["label"] = 0

    outerPosData = outerPosData.loc[:, usedFeatures]
    outerPosData = outerPosData.sample(frac=1, random_state=random_state)
    outerUnlabeledData = outerUnlabeledData.loc[:, usedFeatures]
    outerUnlabeledData = outerUnlabeledData.sample(frac=1, random_state=random_state)

    outerPosData.replace({False: 0, True: 1}, inplace=True)
    outerUnlabeledData.replace({False: 0, True: 1}, inplace=True)
    outerTrainingData = pd.concat([outerPosData, outerUnlabeledData])

    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    colors = ["#1f497d", "#f79646", "#9bbb59", "#7f7f7f", "#8064a2"]
    models = ["RF", "GB", "DT", "SVM", "NB"]
    names = ["RF", "GB", "DT", "SVM", "NB"]

    # BUG-114 (fixed at CP3): one dict per model; dict.fromkeys(models, {}) shared a single dict
    aucScore = {m: {} for m in models}
    puScoreTholdDict = {}
    for scheme in ["bagging"]:
        if drawAUC:
            fig1 = plt.figure(figsize=[6, 6])
            fig1.add_subplot(111, aspect='equal')
        for model, color, name in zip(models, colors, names):
            tprs = []
            mean_fpr = np.linspace(0, 1, 100)
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
                trainingData = trainingData.sample(frac=1, random_state=random_state)

                posIsoDataTest["label"] = 1
                testData = pd.concat([posIsoDataTest, negIsoDataTest])
                testData = testData.loc[:, usedFeatures]
                testData = testData.sample(frac=1, random_state=random_state)

                trainingData.replace({False: 0, True: 1}, inplace=True)
                testData.replace({False: 0, True: 1}, inplace=True)

                sme = SchemeModelEval(trainingData=trainingData, model=model, scheme=scheme, random_state=random_state)
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
                    puScoreTholdDict.setdefault(model, []).append(thold)

                fpr, tpr, t = roc_curve(results.true_label, results.pu_score)
                tprs.append(interp(mean_fpr, fpr, tpr))

            mean_tpr = np.mean(tprs, axis=0)
            mean_auc = auc(mean_fpr, mean_tpr)
            aucScore[model].update({"trainingAUC": mean_auc})
            if drawAUC:
                plt.plot(mean_fpr, mean_tpr, lw=2, linestyle='--', color=color,
                         label='%s (CV AUC = %0.3f )' % (name, mean_auc))

        ######################

        for model, color, name in zip(models, colors, names):
            sme = SchemeModelEval(trainingData=outerTrainingData, model=model, scheme=scheme, random_state=random_state)
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
        if np.mean(list(aucScore[i].values())) > tmpScore:
            tmpModel = i
            tmpScore = np.mean(list(aucScore[i].values()))
    return tmpModel, puScoreTholdDict
