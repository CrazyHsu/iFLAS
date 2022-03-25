#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: generateHtml.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-20
Last modified: 2022-01-20
'''

import fitz
import pandas as pd
import numpy as np
import warnings
import os, argparse, glob, subprocess
from distutils.dir_util import copy_tree
from yattag import Doc, indent
from collections import Counter
from commonFuncs import *


################# Classes ####################
class MakeDict(object):
    def __init__(self, title, url):
        self.title = title
        self.url = url

    def toDict(self):
        return dict([("title", self.title), ("url", self.url)])


################# Functions ##################
def printJson(myDict, outfile=None, indent=4):
    import json
    out = open(outfile, "w")
    print >>out, json.dumps(myDict, indent=indent)
    out.close()


def getRelPath(myPath, targetPath=None, targetDir=None):
    if targetPath and targetDir: return
    if targetPath:
        targetDir = os.path.dirname(targetPath)
        return os.path.relpath(myPath, targetDir)
    if targetDir:
        return os.path.relpath(myPath, targetDir)

def column2breakPoint(colName, tableType="AS"):
    if tableType == "AS":
        myDict = {"ID": "all", "GeneID": "xs", "geneSymbol": "all", "chr": "all", "strand": "all",
                  "exonStart_0base": "all", "exonEnd": "all", "upstreamES": "all", "upstreamEE": "all",
                  "downstreamES": "all",
                  "downstreamEE": "all", "longExonStart_0base": "all", "longExonEnd": "all", "shortES": "all",
                  "shortEE": "all", "flankingES": "all", "flankingEE": "all",
                  "riExonStart_0base": "all", "riExonEnd": "all", "ID.1": "all", "IJC_SAMPLE_1": "xs sm md",
                  "SJC_SAMPLE_1": "xs sm md", "IJC_SAMPLE_2": "xs sm md", "SJC_SAMPLE_2": "xs sm md",
                  "IncFormLen": "all", "SkipFormLen": "all", "PValue": "xs", "FDR": "xs", "IncLevel1": "all",
                  "IncLevel2": "all", "IncLevelDifference": "xs",
                  }
    else:
        myDict = {"ID": "xs", "Ontology": "sm", "Description": "xs", "GeneRatio": "sm", "BgRatio": "sm", "pvalue": "xs",
                  "p.adjust": "xs", "qvalue": "all", "geneID": "all", "Count": "xs"}
    return myDict[colName]


def geneStrucBlock(isoformStrucFile, gene, doc=None, line=None, curDir=None):
    if validateFile(isoformStrucFile):
        line("h1", "The isoform structure in " + gene)
        isoformStrucFile = getRelPath(isoformStrucFile, targetDir=curDir)
        doc.stag("img", klass="img-responsive", src=isoformStrucFile)


def alleleAsBlock(alleleAsFile, doc=None, tag=None, line=None, curDir=None):
    if validateFile(alleleAsFile):
        line("h1", "Allele-specific alternative splicing")
        alleleAsFile = getRelPath(alleleAsFile, targetDir=curDir)
        doc.stag("img", klass="img-responsive", src=alleleAsFile)


def paTailLenAsBlock(paTailLenAsFile, doc=None, tag=None, line=None, curDir=None):
    if validateFile(paTailLenAsFile):
        line("h1", "AS-related poly(A) tail length differential")
        paTailLenAsFile = getRelPath(paTailLenAsFile, targetDir=curDir)
        doc.stag("img", klass="img-responsive", src=paTailLenAsFile)


def paTailLenApaBlock(paTailLenApaFile, doc=None, tag=None, line=None, curDir=None):
    if validateFile(paTailLenApaFile):
        line("h1", "APA-related poly(A) tail length differential")
        paTailLenApaFile = getRelPath(paTailLenApaFile, targetDir=curDir)
        doc.stag("img", klass="img-responsive", src=paTailLenApaFile)


def diffAsBlock(diffAsPlot, doc=None, tag=None, line=None, curDir=None):
    if validateFile(diffAsPlot):
        with tag("div", klass="col"):
            line("h1", "Differential alternative splicing pattern distribution")
            diffAsPlot = getRelPath(diffAsPlot, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=diffAsPlot)


def goEnrichmentBlock(goEnrichPlot, doc=None, tag=None, line=None, curDir=None):
    if validateFile(goEnrichPlot):
        with tag("div", klass="col"):
            line("h1", "GO enrichment of the differential alternative spliced genes")
            goEnrichPlot = getRelPath(goEnrichPlot, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=goEnrichPlot)


def readsCorrAndJuncBlock(basicStatisticsDict, doc=None, tag=None, line=None, curDir=None):
    if "readsCorrection" in basicStatisticsDict:
        readsCorrPlot = basicStatisticsDict["readsCorrection"][1]
        if validateFile(readsCorrPlot):
            with tag("div", klass="col-md-6"):
                line("h3", "The evaluation of reads before and after correction")
                readsCorrPlot = getRelPath(readsCorrPlot, targetDir=curDir)
                doc.stag("img", klass="img-responsive", src=readsCorrPlot)
    if "juncSupported" in basicStatisticsDict:
        juncSupportPlot = basicStatisticsDict["juncSupported"][1]
        if validateFile(juncSupportPlot):
            with tag("div", klass="col-md-6"):
                line("h3", "The junctions of full-length reads supported by NGS reads")
                juncSupportPlot = getRelPath(juncSupportPlot, targetDir=curDir)
                doc.stag("img", klass="img-responsive", src=juncSupportPlot)


def gcContentBlock(basicStatisticsDict, doc=None, tag=None, line=None, curDir=None):
    if "GC_of_raw_flnc" in basicStatisticsDict:
        gcInFlncPlot = basicStatisticsDict["GC_of_raw_flnc"][1]
        if validateFile(gcInFlncPlot):
            with tag("div", klass="col-md-6"):
                line("h3", "The GC content of flnc reads")
                gcInFlncPlot = getRelPath(gcInFlncPlot, targetDir=curDir)
                doc.stag("img", klass="img-responsive", src=gcInFlncPlot)
    if "GC_across_raw_flnc" in basicStatisticsDict:
        gcAcrossFlncPlot = basicStatisticsDict["GC_across_raw_flnc"][1]
        if validateFile(gcAcrossFlncPlot):
            with tag("div", klass="col-md-6"):
                line("h3", "The GC content across flnc reads")
                gcAcrossFlncPlot = getRelPath(gcAcrossFlncPlot, targetDir=curDir)
                doc.stag("img", klass="img-responsive", src=gcAcrossFlncPlot)


def asPatternBlock(basicStatisticsDict, doc=None, tag=None, line=None, curDir=None):
    annotationPlot = basicStatisticsDict["asPattern"]["asAnno"][1]
    spliceSitePlot = basicStatisticsDict["asPattern"]["asSpliceSite"][1]
    if validateFile(annotationPlot):
        with tag("div", klass="col-md-6"):
            line("h3", "The alternative splicing(AS) summary")
            annotationPlot = getRelPath(annotationPlot, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=annotationPlot)
    if validateFile(spliceSitePlot):
        with tag("div", klass="col-md-6"):
            line("h3", "The splice-site of alternative splicing(AS)")
            spliceSitePlot = getRelPath(spliceSitePlot, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=spliceSitePlot)

def lengthDistributionBlock(basicStatisticsDict, doc=None, tag=None, line=None, curDir=None):
    lenDistCurve = basicStatisticsDict["LengthDistribution"][1]
    lenDistBox = basicStatisticsDict["LengthDistribution"][2]
    if validateFile(lenDistBox):
        with tag("div", klass="col-md-6"):
            line("h3", "The length distribution (box) of full-length reads")
            lenDistBox = getRelPath(lenDistBox, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=lenDistBox)
    if validateFile(lenDistCurve):
        with tag("div", klass="col-md-6"):
            line("h3", "The length distribution (curve) of full-length reads")
            lenDistCurve = getRelPath(lenDistCurve, targetDir=curDir)
            doc.stag("img", klass="img-responsive", src=lenDistCurve)


def generateTable(fileIn, tableType="AS", doc=None, tag=None, text=None, line=None):
    df = pd.read_csv(fileIn, sep="\t")
    with tag("table", ("class", "table"), ("data-show-toggle", "true"), ("data-paging", "true"),
             ("data-filtering", "true"), ("data-sorting", "true")):
        with tag("thead"):
            with tag("tr", style="background-color: #c3e6cb"):
                for i in df.columns:
                    breakPoint = column2breakPoint(i, tableType=tableType)
                    line("th", i, ("data-breakpoints", breakPoint))
        with tag("tbody"):
            for idx, row in df.iterrows():
                with tag("tr"):
                    for i in row.tolist():
                        if isinstance(i, float):
                            formattedValue = ("%.25f" % i).rstrip('0').rstrip('.')
                            if "e" in str(i):
                                line("td", "{:.3e}".format(i), ("data-sort-value", formattedValue))
                            else:
                                line("td", "{:.6f}".format(i), ("data-sort-value", formattedValue))
                        else:
                            line("td", i)


def generateMainPage(reportDict):
    doc, tag, text, line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.stag('meta', charset='utf-8')
            doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
            line('title', 'iFLAS main page')
            doc.stag('link', rel='stylesheet', href='assets/css/bootstrap/css/bootstrap.min.css')
            doc.stag('link', rel='stylesheet', href='assets/css/font-awesome/css/font-awesome.min.css')
            doc.stag('link', rel='stylesheet', href="assets/css/iflas.css")
            doc.stag('link', rel='stylesheet', href="assets/css/footable.bootstrap.css")

        with tag("body"):
            with tag("div", id="sidebar-test"):
                with tag("div", klass="sidebar-header"):
                    with tag("h2"):
                        line("a", "iFLAS", href="iflas_report.html", klass="iflas")
                with tag("ul"):
                    getSideBar(reportDict, doc=doc, tag=tag, text=text, line=line, mainPage=True)

            with tag("div", klass="content"):
                with tag("div", klass="container"):
                    with tag("div", klass="row"):
                        with tag("div", klass="col"):
                            line("h1", "iFLAS: integrated Full Length Alternative Splicing analysis")
                            description = '''
                                iFLAS is an integrated pipeline for the identification, visualization and exploration of functional AS events by taking advantages of long-read RNA sequencing data. 
                                iFLAS is able to improve transcripts annotation, detect AS/polyadenylation events related differential poly(A) tails, identify haplotype-specific AS events, and detect differential AS events. 
                                Furthermore, iFLAS can generate rich graphical plots and comprehensive summary report, making it easy-to-use for researchers.
                                Below is the workflow of iFLAS.
                            '''
                            with tag("h3"):
                                text(description)
                            workflow_png = "assets/src/iFLAS_workflow.png"
                            doc.stag("img", klass="img-responsive", src=workflow_png)
                with tag("div", klass="footer"):
                    pass
            line("script", "", src="assets/js/jquery.min.js")
            line("script", "", src="assets/js/bootstrap.min.js")
            line("script", "", src="assets/js/iflas.js")
            line("script", "", src="assets/js/footable.js")

    mainPageOut = open("iflas_report.html", "w")
    res = indent(doc.getvalue(), indentation="    ")
    mainPageOut.write(res)
    mainPageOut.close()
    return MakeDict("mainPage", os.path.join(os.getcwd(), "iflas_report.html")).toDict()


def getSideBar(reportDict, sample=None, gene=None, basicStatistics=False, das=False, dasComp=None, doc=None, tag=None,
               text=None, line=None, mainPage=False):
    keys = reportDict.keys()
    newKeys = [x for x in keys if x != "das" and x != "allSampleMerged"]
    if "allSampleMerged" in keys:
        newKeys.insert(0, "allSampleMerged")
    if "das" in keys:
        newKeys.append("das")
    for tmpSample in newKeys:
        if sample and sample == tmpSample:
            expand = "true"
            faPlusOrMinus = "icon-minus"
            myClass = "list-unstyled collapse in"
            active = " active"
        else:
            expand = "false"
            faPlusOrMinus = "icon-plus"
            myClass = "list-unstyled collapse"
            active = ""
        if mainPage:
            relativeDir = os.path.join("", tmpSample)
        else:
            if sample == tmpSample:
                relativeDir = ""
            else:
                relativeDir = os.path.join("../", tmpSample)

        if tmpSample == "das":
            with tag("li"):
                with tag("div", ("class", "sample" + active), ("href", "#" + tmpSample), ("data-toggle", "collapse"),
                         ("aria-expanded", expand)):
                    with tag("a", ("class", "sectionName")):
                        text("Differential Alternative Splicing")
                    with tag("a", ("href", "#das"), ("data-toggle", "collapse"), ("aria-expanded", expand)):
                        line("i", "", klass=faPlusOrMinus)
                with tag("ul", ("class", myClass), ("id", "das"), ("aria-expanded", expand)):
                    for i in reportDict[tmpSample]:
                        with tag("li"):
                            href = os.path.join(relativeDir, i + ".html")
                            if das and i == dasComp:
                                with tag("a", ("href", href), ("class", "active"), ("aria-selected", "true")):
                                    text(" " + i)
                            else:
                                with tag("a", ("href", href), ("aria-selected", "false")):
                                    text(" " + i)
        elif tmpSample == "allSampleMerged":
            with tag("li"):
                with tag("div", ("class", "sample" + active), ("href", "#" + tmpSample), ("data-toggle", "collapse"),
                         ("aria-expanded", expand)):
                    with tag("a", ("class", "sectionName")):
                        text("All Sample Merged")
                    with tag("a", ("href", "#" + tmpSample), ("data-toggle", "collapse"), ("aria-expanded", expand)):
                        line("i", "", klass=faPlusOrMinus)
                with tag("ul", ("class", myClass), ("id", tmpSample), ("aria-expanded", expand)):
                    for i in reportDict[tmpSample]:
                        with tag("li"):
                            href = os.path.join(relativeDir, i + ".html")
                            if gene and gene == i and tmpSample == sample:
                                with tag("a", ("href", href), ("class", "active"), ("aria-selected", "true")):
                                    text(" " + i)
                            else:
                                with tag("a", ("href", href), ("aria-selected", "false")):
                                    text(" " + i)
        else:
            with tag("li"):
                with tag("div", ("class", "sample" + active), ("href", "#" + tmpSample), ("data-toggle", "collapse"),
                         ("aria-expanded", expand)):
                    with tag("a", ("class", "sectionName")):
                        text(tmpSample.capitalize())
                    with tag("a", ("href", "#" + tmpSample), ("data-toggle", "collapse"), ("aria-expanded", expand)):
                        line("i", "", klass=faPlusOrMinus)
                with tag("ul", ("class", myClass), ("id", tmpSample), ("aria-expanded", expand)):
                    if "basicStatistics" in reportDict[tmpSample]:
                        with tag("li"):
                            href = os.path.join(relativeDir, "basicStatistics.html")
                            if basicStatistics and tmpSample == sample:
                                with tag("a", ("href", href), ("class", "active"), ("aria-selected", "true")):
                                    text(" Basic Statistics")
                            else:
                                with tag("a", ("href", href), ("aria-selected", "false")):
                                    text(" Basic Statistics")
                    if "genes" in reportDict[tmpSample]:
                        for i in reportDict[tmpSample]["genes"]:
                            with tag("li"):
                                href = os.path.join(relativeDir, i + ".html")
                                if gene and gene == i and tmpSample == sample:
                                    with tag("a", ("href", href), ("class", "active"), ("aria-selected", "true")):
                                        text(" " + i)
                                else:
                                    with tag("a", ("href", href), ("aria-selected", "false")):
                                        text(" " + i)


def generateAllSampleMergedPage(reportDict):
    dataList = []
    for gene in reportDict["allSampleMerged"]:
        resolveDir("allSampleMerged", chdir=False)
        out = open(os.path.join("allSampleMerged", gene + ".html"), "w")
        doc, tag, text, line = Doc().ttl()
        doc.asis('<!DOCTYPE html>')
        with tag('html'):
            with tag('head'):
                doc.stag('meta', charset='utf-8')
                doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                line('title', 'The gene structure of all merged samples')
                doc.stag('link', rel='stylesheet', href='../assets/css/bootstrap/css/bootstrap.min.css')
                doc.stag('link', rel='stylesheet', href='../assets/css/font-awesome/css/font-awesome.min.css')
                doc.stag('link', rel='stylesheet', href="../assets/css/iflas.css")
                doc.stag('link', rel='stylesheet', href="../assets/css/footable.bootstrap.css")
            with tag("body"):
                with tag("div", id="sidebar-test"):
                    with tag("div", klass="sidebar-header"):
                        with tag("h2"):
                            line("a", "iFLAS", href="../iflas_report.html", klass="iflas")
                    with tag("ul"):
                        getSideBar(reportDict, sample="allSampleMerged", gene=gene, doc=doc, tag=tag, text=text,
                                   line=line)

                with tag("div", klass="content"):
                    with tag("div", klass="container"):
                        with tag("div", klass="row"):
                            with tag("div", klass="col"):
                                curDir = os.path.join(os.getcwd(), "allSampleMerged")
                                geneStrucBlock(reportDict["allSampleMerged"][gene], gene=gene, doc=doc, line=line,
                                               curDir=curDir)
                        with tag("div", klass="footer"):
                            pass

                line("script", "", src="../assets/js/jquery.min.js")
                line("script", "", src="../assets/js/bootstrap.min.js")
                line("script", "", src="../assets/js/iflas.js")
                line("script", "", src="../assets/js/footable.js")

                customJs = '''
                                <script>
                                    $('.table').footable();

                                    var offestFromTop = %d * 43 + 65;
                                    $('#sidebar-test').scrollTop(offestFromTop);

                                </script>
                            ''' % (0)
                doc.asis(customJs)

        res = indent(doc.getvalue(), indentation="    ")
        out.write(res)
        out.close()
        dataList.append(MakeDict("merged sample for {}".format(gene), os.path.join(os.getcwd(), "allSampleMerged", gene + ".html")).toDict())
    return dataList


def generateDasPage(reportDict):
    dataList = []
    for tmpSample in reportDict["das"]:
        dasDict = reportDict["das"][tmpSample]
        resolveDir("das", chdir=False)
        curDir = os.path.join(os.getcwd(), "das")
        out = open(os.path.join("das", tmpSample + ".html"), "w")
        doc, tag, text, line = Doc().ttl()
        doc.asis('<!DOCTYPE html>')
        with tag('html'):
            with tag('head'):
                doc.stag('meta', charset='utf-8')
                doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                line('title', 'The summary of differential alternative splicing events in ' + tmpSample)
                doc.stag('link', rel='stylesheet', href='../assets/css/bootstrap/css/bootstrap.min.css')
                doc.stag('link', rel='stylesheet', href='../assets/css/font-awesome/css/font-awesome.min.css')
                doc.stag('link', rel='stylesheet', href="../assets/css/iflas.css")
                doc.stag('link', rel='stylesheet', href="../assets/css/footable.bootstrap.css")
            with tag("body"):
                with tag("div", id="sidebar-test"):
                    with tag("div", klass="sidebar-header"):
                        with tag("h2"):
                            line("a", "iFLAS", href="../iflas_report.html", klass="iflas")
                    with tag("ul"):
                        getSideBar(reportDict, sample="das", das=True, dasComp=tmpSample, doc=doc, tag=tag, text=text,
                                   line=line)

                with tag("div", klass="content"):
                    with tag("div", klass="container"):
                        if "dasDistribution" in dasDict:
                            with tag("div", klass="row"):
                                diffAsBlock(dasDict["dasDistribution"], doc=doc, tag=tag, line=line, curDir=curDir)
                        line("h1", "Differential alternative spliced events by categories")
                        if "IR" in dasDict:
                            asType = "ir"
                            with tag("div", ("class", "row"), ("id", asType)):
                                with tag("div", ("class", "col asEvent")):
                                    with tag("button", ("class", "btn btn-primary"), ("type", "button"),
                                             ("data-toggle", "collapse"), ("data-target", "#{}Collapse".format(asType)),
                                             ("aria-expanded", "false"), ("aria-controls", "{}Collapse".format(asType))):
                                        text("Intron Retention")
                                with tag("div", ("class", "col collapse"), ("id", "{}Collapse".format(asType))):
                                    generateTable(dasDict["IR"], tableType="AS", doc=doc, tag=tag, text=text, line=line)
                        if "SE" in dasDict:
                            asType = "se"
                            with tag("div", ("class", "row"), ("id", asType)):
                                with tag("div", ("class", "col asEvent")):
                                    with tag("button", ("class", "btn btn-primary"), ("type", "button"),
                                             ("data-toggle", "collapse"), ("data-target", "#{}Collapse".format(asType)),
                                             ("aria-expanded", "false"), ("aria-controls", "{}Collapse".format(asType))):
                                        text("Exon Skipping")
                                with tag("div", ("class", "col collapse"), ("id", "{}Collapse".format(asType))):
                                    generateTable(dasDict["SE"], tableType="AS", doc=doc, tag=tag, text=text, line=line)
                        if "A5SS" in dasDict:
                            asType = "a5ss"
                            with tag("div", ("class", "row"), ("id", asType)):
                                with tag("div", ("class", "col asEvent")):
                                    with tag("button", ("class", "btn btn-primary"), ("type", "button"),
                                             ("data-toggle", "collapse"), ("data-target", "#{}Collapse".format(asType)),
                                             ("aria-expanded", "false"), ("aria-controls", "{}Collapse".format(asType))):
                                        text("Alternative Donor(A5SS)")
                                with tag("div", ("class", "col collapse"), ("id", "{}Collapse".format(asType))):
                                    generateTable(dasDict["A5SS"], tableType="AS", doc=doc, tag=tag, text=text, line=line)
                        if "A3SS" in dasDict:
                            asType = "a3ss"
                            with tag("div", ("class", "row"), ("id", asType)):
                                with tag("div", ("class", "col asEvent")):
                                    with tag("button", ("class", "btn btn-primary"), ("type", "button"),
                                             ("data-toggle", "collapse"), ("data-target", "#{}Collapse".format(asType)),
                                             ("aria-expanded", "false"), ("aria-controls", "{}Collapse".format(asType))):
                                        text("Alternative Acceptor(A3SS)")
                                with tag("div", ("class", "col collapse"), ("id", "{}Collapse".format(asType))):
                                    generateTable(dasDict["A3SS"], tableType="AS", doc=doc, tag=tag, text=text, line=line)
                        if "goEnrichPlot" in dasDict:
                            with tag("div", klass="row"):
                                goEnrichmentBlock(dasDict["goEnrichPlot"], doc=doc, tag=tag, line=line, curDir=curDir)
                        if "goEnrichResults" in dasDict:
                            line("h1", "GO enrichment results of the differential alternative spliced genes")
                            with tag("div", ("class", "row")):
                                with tag("div", ("class", "col enrichResults")):
                                    with tag("button", ("class", "btn btn-primary"), ("type", "button"),
                                             ("data-toggle", "collapse"), ("data-target", "#enrichResults"),
                                             ("aria-expanded", "false"), ("aria-controls", "enrichResults")):
                                        text("GO enrichment")
                                with tag("div", ("class", "col collapse"), ("id", "enrichResults")):
                                    generateTable(dasDict["goEnrichResults"], tableType="goEnrichment", doc=doc, tag=tag,
                                                  text=text, line=line)
                        doc.asis("<hr>")
                    with tag("div", klass="footer"):
                        pass

                line("script", "", src="../assets/js/jquery.min.js")
                line("script", "", src="../assets/js/bootstrap.min.js")
                line("script", "", src="../assets/js/iflas.js")
                line("script", "", src="../assets/js/footable.js")

                customJs = '''
                                <script>
                                    $('.table').footable();
    
                                    var offestFromTop = %d * 43 + 65;
                                    $('#sidebar-test').scrollTop(offestFromTop);
    
                                </script>
                            ''' % (len(reportDict.keys()) - 1)
                doc.asis(customJs)

        res = indent(doc.getvalue(), indentation="    ")
        out.write(res)
        out.close()
        dataList.append(MakeDict("differential alternavive spliced page for {}".format(tmpSample), os.path.join(os.getcwd(), "das", tmpSample + ".html")).toDict())
    return dataList


def generateSingleSampleStatisticsPage(reportDict):
    dataList = []
    keys = reportDict.keys()
    newKeys = [x for x in keys if x != "das" and x != "allSampleMerged"]
    count = 0
    for tmpSample in newKeys:
        if tmpSample == "das" or tmpSample == "allSampleMerged": continue
        if "basicStatistics" not in reportDict[tmpSample]: continue
        basicStatisticsDict = reportDict[tmpSample]["basicStatistics"]
        resolveDir(tmpSample, chdir=False)
        curDir = os.path.join(os.getcwd(), tmpSample)
        out = open(os.path.join(tmpSample, "basicStatistics.html"), "w")
        doc, tag, text, line = Doc().ttl()
        doc.asis('<!DOCTYPE html>')
        with tag('html'):
            with tag('head'):
                doc.stag('meta', charset='utf-8')
                doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                line('title', 'Basic statistics in sample ' + tmpSample)
                doc.stag('link', rel='stylesheet', href='../assets/css/bootstrap/css/bootstrap.min.css')
                doc.stag('link', rel='stylesheet', href='../assets/css/font-awesome/css/font-awesome.min.css')
                doc.stag('link', rel='stylesheet', href="../assets/css/iflas.css")
                doc.stag('link', rel='stylesheet', href="../assets/css/footable.bootstrap.css")
            with tag("body"):
                with tag("div", id="sidebar-test"):
                    with tag("div", klass="sidebar-header"):
                        with tag("h2"):
                            line("a", "iFLAS", href="../iflas_report.html", klass="iflas")
                    with tag("ul"):
                        getSideBar(reportDict, sample=tmpSample, basicStatistics=True, doc=doc, tag=tag, text=text,
                                   line=line)

                with tag("div", klass="content"):
                    with tag("div", klass="container"):
                        if "readsCorrection" in basicStatisticsDict or "juncSupported" in basicStatisticsDict:
                            with tag("div", klass="row"):
                                readsCorrAndJuncBlock(basicStatisticsDict, doc=doc, tag=tag, line=line, curDir=curDir)
                        if "GC_of_raw_flnc" in basicStatisticsDict or "GC_across_raw_flnc" in basicStatisticsDict:
                            with tag("div", klass="row"):
                                gcContentBlock(basicStatisticsDict, doc=doc, tag=tag, line=line, curDir=curDir)
                        if "LengthDistribution" in basicStatisticsDict:
                            with tag("div", klass="row"):
                                lengthDistributionBlock(basicStatisticsDict, doc=doc, tag=tag, line=line, curDir=curDir)
                        if "asPattern" in basicStatisticsDict:
                            with tag("div", klass="row"):
                                asPatternBlock(basicStatisticsDict, doc=doc, tag=tag, line=line, curDir=curDir)
                    with tag("div", klass="footer"):
                        pass

                line("script", "", src="../assets/js/jquery.min.js")
                line("script", "", src="../assets/js/bootstrap.min.js")
                line("script", "", src="../assets/js/iflas.js")
                line("script", "", src="../assets/js/footable.js")

                customJs = '''
                                <script>
                                    $('.table').footable();
            
                                    var offestFromTop = %d * 43 + 65;
                                    $('#sidebar-test').scrollTop(offestFromTop);
            
                                </script>
                            ''' % (count + 1 if "allSampleMerged" in reportDict else count)
                doc.asis(customJs)

        res = indent(doc.getvalue(), indentation="    ")
        out.write(res)
        out.close()

        count += 1
        dataList.append(MakeDict("basic statistics page for {}".format(tmpSample), os.path.join(os.getcwd(), tmpSample, "basicStatistics.html")).toDict())
    return dataList


def generateSingleSampleGenePage(reportDict):
    dataList = []
    keys = reportDict.keys()
    newKeys = [x for x in keys if x != "das" and x != "allSampleMerged"]
    count = 0
    for tmpSample in newKeys:
        if tmpSample == "das" or tmpSample == "allSampleMerged": continue
        if "genes" not in reportDict[tmpSample]: continue
        # allGenes = sorted(reportDict[tmpSample]["genes"])
        for gene in reportDict[tmpSample]["genes"]:
            resolveDir(tmpSample, chdir=False)
            curDir = os.path.join(os.getcwd(), tmpSample)
            out = open(os.path.join(tmpSample, gene + ".html"), "w")
            doc, tag, text, line = Doc().ttl()
            doc.asis('<!DOCTYPE html>')
            with tag('html'):
                with tag('head'):
                    doc.stag('meta', charset='utf-8')
                    doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                    line('title', 'Detailed information in gene ' + gene)
                    doc.stag('link', rel='stylesheet', href='../assets/css/bootstrap/css/bootstrap.min.css')
                    doc.stag('link', rel='stylesheet', href='../assets/css/font-awesome/css/font-awesome.min.css')
                    doc.stag('link', rel='stylesheet', href="../assets/css/iflas.css")
                    doc.stag('link', rel='stylesheet', href="../assets/css/footable.bootstrap.css")
                with tag("body"):
                    with tag("div", id="sidebar-test"):
                        with tag("div", klass="sidebar-header"):
                            with tag("h2"):
                                line("a", "iFLAS", href="../iflas_report.html", klass="iflas")
                        with tag("ul"):
                            getSideBar(reportDict, sample=tmpSample, gene=gene, doc=doc, tag=tag, text=text, line=line)

                    with tag("div", klass="content"):
                        with tag("div", klass="container"):
                            if "isoformStruc" in reportDict[tmpSample]["genes"][gene]:
                                with tag("div", klass="row"):
                                    with tag("div", klass="col"):
                                        geneStrucBlock(reportDict[tmpSample]["genes"][gene]["isoformStruc"], gene=gene,
                                                       doc=doc, line=line, curDir=curDir)
                            if "alleleAS" in reportDict[tmpSample]["genes"][gene]:
                                with tag("div", klass="row"):
                                    with tag("div", klass="col"):
                                        for i in reportDict[tmpSample]["genes"][gene]["alleleAS"]:
                                            alleleAsBlock(i, doc=doc, tag=tag, line=line, curDir=curDir)
                            if "palenAS" in reportDict[tmpSample]["genes"][gene]:
                                with tag("div", klass="row"):
                                    with tag("div", klass="col"):
                                        for i in reportDict[tmpSample]["genes"][gene]["palenAS"]:
                                            paTailLenAsBlock(i, doc=doc, tag=tag, line=line, curDir=curDir)
                            if "palenAPA" in reportDict[tmpSample]["genes"][gene]:
                                with tag("div", klass="row"):
                                    with tag("div", klass="col"):
                                        for i in reportDict[tmpSample]["genes"][gene]["palenAPA"]:
                                            paTailLenApaBlock(i, doc=doc, tag=tag, line=line, curDir=curDir)
                    with tag("div", klass="footer"):
                        pass

                    line("script", "", src="../assets/js/jquery.min.js")
                    line("script", "", src="../assets/js/bootstrap.min.js")
                    line("script", "", src="../assets/js/iflas.js")
                    line("script", "", src="../assets/js/footable.js")

                    customJs = '''
                                    <script>
                                        $('.table').footable();

                                        var offestFromTop = %d * 43 + 65;
                                        $('#sidebar-test').scrollTop(offestFromTop);

                                    </script>
                                ''' % (count + 1 if "allSampleMerged" in reportDict else count)
                    doc.asis(customJs)

            res = indent(doc.getvalue(), indentation="    ")
            out.write(res)
            out.close()
            dataList.append(MakeDict("single gene page of {} for {}".format(gene, tmpSample), os.path.join(os.getcwd(), tmpSample, gene + ".html")).toDict())
        count += 1
    return dataList


def convertPdf2png(inPdf=None, outDir=None, pageIndex=0):
    if outDir != None:
        outPNG = os.path.join(outDir, "{}.png".format(os.path.splitext(os.path.basename(inPdf))[0]))
    else:
        outPNG = os.path.join(os.path.dirname(inPdf), "{}.png".format(os.path.splitext(os.path.basename(inPdf))[0]))
    doc = fitz.open(inPdf)
    page = doc.loadPage(pageIndex)
    zoom_x = 5.0  # horizontal zoom
    zomm_y = 5.0  # vertical zoom
    mat = fitz.Matrix(zoom_x, zomm_y)  # zoom factor 2 in each dimension
    pix = page.getPixmap(matrix=mat)
    pix.writePNG(outPNG)
    return outPNG

# def convert(reportDict, iflasOut):
#     for sample in reportDict:
#         basicStatisticsDir = os.path.join(iflasOut, sample, "basicStatistics")
#         geneDir = os.path.join(iflasOut, sample, "genes")
#         resolveDir(basicStatisticsDir, chdir=False)
#         resolveDir(geneDir, chdir=False)
#         if sample == "mergedSample":
#             diffASDir = os.path.join(basicStatisticsDir, "diffAS")
#             goEnrichDir = os.path.join(basicStatisticsDir, "goEnrich")
#             resolveDir(diffASDir, chdir=False)
#             resolveDir(goEnrichDir, chdir=False)
#             diffASPdf = reportDict[sample]["basicStatistics"]["diffAS"][1]
#             goEnrichPdf = reportDict[sample]["basicStatistics"]["goEnrich"][1]
#             reportDict[sample]["basicStatistics"]["diffAS"][1] = convertPdf2png(diffASPdf, diffASDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["goEnrich"][1] = convertPdf2png(goEnrichPdf, goEnrichDir, pageIndex=0)
#         else:
#             readsCorrectionDir = os.path.join(basicStatisticsDir, "readsCorrection")
#             reportReadsContentEvalDir = os.path.join(basicStatisticsDir, "reportReadsContentEval")
#             asPatternDir = os.path.join(basicStatisticsDir, "asPattern")
#             isoformRankDir = os.path.join(basicStatisticsDir, "isoformRank")
#             resolveDir(readsCorrectionDir, chdir=False)
#             resolveDir(reportReadsContentEvalDir, chdir=False)
#             resolveDir(asPatternDir, chdir=False)
#             resolveDir(isoformRankDir, chdir=False)
#             readsCorrectionPdf = reportDict[sample]["basicStatistics"]["readsCorrection"][1]
#             gcFlncPdf = reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["GC_of_raw_flnc"][1]
#             gcAcrossFlncPdf = reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["GC_across_raw_flnc"][1]
#             LengthDistributionPdf1 = reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["LengthDistribution"][0]
#             LengthDistributionPdf2 = reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["LengthDistribution"][1]
#             asAnnoPdf = reportDict[sample]["basicStatistics"]["asPattern"]["asAnno"][1]
#             asSpliceSitePdf = reportDict[sample]["basicStatistics"]["asPattern"]["asSpliceSite"][1]
#             isoformRankPdf = reportDict[sample]["basicStatistics"]["isoformRank"][0]
#             reportDict[sample]["basicStatistics"]["readsCorrection"][1] = convertPdf2png(readsCorrectionPdf, readsCorrectionDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["GC_of_raw_flnc"][1] = convertPdf2png(gcFlncPdf, reportReadsContentEvalDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["GC_across_raw_flnc"][1] = convertPdf2png(gcAcrossFlncPdf, reportReadsContentEvalDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["LengthDistribution"][0] = convertPdf2png(LengthDistributionPdf1, reportReadsContentEvalDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["reportReadsContentEval"]["LengthDistribution"][1] = convertPdf2png(LengthDistributionPdf2, reportReadsContentEvalDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["asPattern"]["asAnno"][1] = convertPdf2png(asAnnoPdf, asPatternDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["asPattern"]["asSpliceSite"][1] = convertPdf2png(asSpliceSitePdf, asPatternDir, pageIndex=0)
#             reportDict[sample]["basicStatistics"]["isoformRank"][0] = convertPdf2png(isoformRankPdf, isoformRankDir, pageIndex=0)
#
#         for gene in reportDict[sample]:
#             targetGeneDir = os.path.join(geneDir, gene)
#             resolveDir(targetGeneDir, chdir=False)
#             if "isoformStruc" in reportDict[sample][gene]:
#                 isoformStrucPdf = reportDict[sample][gene]["isoformStruc"]
#                 reportDict[sample][gene]["isoformStruc"] = convertPdf2png(isoformStrucPdf, targetGeneDir, pageIndex=0)
#             if "alleleAS" in reportDict[sample][gene]:
#                 alleleASPdf = reportDict[sample][gene]["alleleAS"]
#                 reportDict[sample][gene]["isoformStruc"] = convertPdf2png(alleleASPdf, targetGeneDir, pageIndex=0)
#             if "paTailAS" in reportDict[sample][gene]:
#                 paTailASPdf = reportDict[sample][gene]["paTailAS"]
#                 reportDict[sample][gene]["isoformStruc"] = convertPdf2png(paTailASPdf, targetGeneDir, pageIndex=0)

def retrieveResults(dataToProcess, dirSpec, optionTools, args):
    reportDir = os.path.join(dirSpec.out_dir, "reports")
    resultDict = {}
    for dataObj in dataToProcess:
        projectName, sampleName = dataObj.project_name, dataObj.sample_name
        subReportDir = os.path.join(reportDir, "{}_{}".format(projectName, sampleName))
        resultDict.update({sampleName: {"basicStatistics": {}, "genes": {}}})

        if os.path.exists(os.path.join(subReportDir, "basicStatistics")):
            if os.path.exists(os.path.join(subReportDir, "basicStatistics", "readsCorrectResult.pdf")):
                readsCorrectResult = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "readsCorrectResult.pdf"))
                fileDiscription = "Reads correction evaluation in sample {}".format(sampleName)
                resultDict[sampleName]["basicStatistics"].update({"readsCorrection": [fileDiscription, readsCorrectResult]})

            gcOfFlncReads = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "GC_of_raw_flnc_reads.pdf"))
            gcAcrossFlncReads = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "GC_across_raw_flnc_read.pdf"))
            lengthDistributionCurve = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "LengthDistribution.curve.pdf"))
            lengthDistributionBox = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "LengthDistribution.box.pdf"))

            fileDiscription = ""
            resultDict[sampleName]["basicStatistics"].update({"GC_of_raw_flnc": [fileDiscription, gcOfFlncReads]})
            fileDiscription = ""
            resultDict[sampleName]["basicStatistics"].update({"GC_across_raw_flnc": [fileDiscription, gcAcrossFlncReads]})
            fileDiscription = ""
            resultDict[sampleName]["basicStatistics"].update({"LengthDistribution": [fileDiscription, lengthDistributionCurve, lengthDistributionBox]})
            if os.path.exists(os.path.join(subReportDir, "basicStatistics", "supportedByRNAseq.pdf")):
                juncSupportation = convertPdf2png(inPdf=os.path.join(subReportDir, "basicStatistics", "supportedByRNAseq.pdf"))
                fileDiscription = ""
                resultDict[sampleName]["basicStatistics"].update({"juncSupported": [fileDiscription, juncSupportation]})

        if os.path.exists(os.path.join(subReportDir, "asPattern")):
            asAnnoPdf = convertPdf2png(inPdf=os.path.join(subReportDir, "asPattern", "{}_{}.AS_annotation.pdf".format(projectName, sampleName)))
            asSpliceSitePdf = convertPdf2png(inPdf=os.path.join(subReportDir, "asPattern", "{}_{}.AS_spliceSite.pdf".format(projectName, sampleName)))
            fileDiscription = ""
            resultDict[sampleName]["basicStatistics"]["asPattern"] = {"asAnno": [fileDiscription, asAnnoPdf]}
            fileDiscription = ""
            resultDict[sampleName]["basicStatistics"]["asPattern"].update({"asSpliceSite": [fileDiscription, asSpliceSitePdf]})

        geneStrucDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "isoViewer")
        if os.path.exists(geneStrucDir):
            geneStrucReportDir = os.path.join(subReportDir, "geneStrucPlots")
            for i in glob.glob(os.path.join(geneStrucDir, "*")):
                gene = os.path.basename(i)
                plotPath = convertPdf2png(inPdf=os.path.join(i, "{}.pdf".format(gene)), outDir=geneStrucReportDir)
                if gene not in resultDict[sampleName]["genes"]:
                    resultDict[sampleName]["genes"].update({gene: {"isoformStruc": plotPath}})
                else:
                    resultDict[sampleName]["genes"][gene].update({"isoformStruc": plotPath})

        alleleAsDir = os.path.join(subReportDir, "alleleAsPlots")
        if os.path.exists(alleleAsDir):
            for i in glob.glob(os.path.join(alleleAsDir, "*.pdf")):
                if i == "alleleAS.pdf": continue
                gene = os.path.basename(i).split(".")[0]
                plotPath = convertPdf2png(inPdf=i)
                if gene not in resultDict[sampleName]["genes"]:
                    resultDict[sampleName]["genes"].update({gene: {"alleleAS": [plotPath]}})
                elif "alleleAS" not in resultDict[sampleName]["genes"][gene]:
                    resultDict[sampleName]["genes"][gene] = {"alleleAS": [plotPath]}
                else:
                    resultDict[sampleName]["genes"][gene]["alleleAS"].append(plotPath)

        palenAsDir = os.path.join(subReportDir, "palenAsPlots")
        if os.path.exists(palenAsDir):
            for i in glob.glob(os.path.join(palenAsDir, "*.pdf")):
                gene = os.path.basename(i).split(".palenAndAS.pdf")[0]
                plotPath = convertPdf2png(inPdf=i)
                if gene not in resultDict[sampleName]["genes"]:
                    resultDict[sampleName]["genes"].update({gene: {"palenAS": [plotPath]}})
                elif "palenAS" not in resultDict[sampleName]["genes"][gene]:
                    resultDict[sampleName]["genes"][gene].update({"palenAS": [plotPath]})
                else:
                    resultDict[sampleName]["genes"][gene]["palenAS"].append(plotPath)

        palenApaDir = os.path.join(subReportDir, "palenApaPlots")
        if os.path.exists(palenApaDir):
            for i in glob.glob(os.path.join(palenApaDir, "*.pdf")):
                gene = os.path.basename(i).split(".palenAndAPA.pdf")[0]
                plotPath = convertPdf2png(inPdf=i)
                if gene not in resultDict[sampleName]["genes"]:
                    resultDict[sampleName]["genes"].update({gene: {"palenAPA": [plotPath]}})
                elif "palenAPA" not in resultDict[sampleName]["genes"][gene]:
                    resultDict[sampleName]["genes"][gene].update({"palenAPA": [plotPath]})
                else:
                    resultDict[sampleName]["genes"][gene]["palenAPA"].append(plotPath)

    mergedGeneStrucDir = os.path.join(dirSpec.out_dir, "isoViewer_sample_merged")
    if os.path.exists(mergedGeneStrucDir):
        mergedGeneStrucReportDir = os.path.join(reportDir, "allSampleMergedGeneStrucPlots")
        resolveDir(mergedGeneStrucReportDir, chdir=False)
        resultDict.update({"allSampleMerged": {}})
        for i in glob.glob(os.path.join(mergedGeneStrucDir, "*")):
            gene = os.path.basename(i)
            plotPath = convertPdf2png(inPdf=os.path.join(i, "{}.pdf".format(gene)), outDir=mergedGeneStrucReportDir)
            resultDict["allSampleMerged"].update({gene: plotPath})

    dasReportDir = os.path.join(dirSpec.out_dir, "das")
    if os.path.exists(dasReportDir):
        dasDir = os.path.join(dirSpec.out_dir, "das")
        validAsType = ["IR", "SE", "A3SS", "A5SS"]
        resultDict.update({"das": {}})
        for i in glob.glob(os.path.join(dasDir, "*")):
            if i.endswith(".sigDiffAS"):
                compPair = os.path.basename(i).split(".sigDiffAS")[0]
                resolveDir(os.path.join(dasReportDir, compPair))
                resultDict["das"].update({compPair: {}})
                tmpSigAsFiles = []
                for j in glob.glob(os.path.join(i, "*")):
                    asType = os.path.basename(j).split(".")[0]
                    if asType not in validAsType: continue
                    makeLink(j, os.path.join(dasReportDir, compPair, os.path.basename(j)))
                    tmpSigAsFiles.append(os.path.join(dasReportDir, compPair, os.path.basename(j)))
                    resultDict["das"][compPair].update({asType: j})
                dasDistribution = convertPdf2png(inPdf=os.path.join(i, "sigDiff.AS_distribution.pdf"), outDir=i)
                resultDict["das"][compPair].update({"dasDistribution": dasDistribution})

                # cmd = "cat {} | grep -v 'ID' | cut -f 2 | sort -u > dasg.lst".format(" ".join(tmpSigAsFiles))
                # subprocess.call(cmd, shell=True)
                # from plotRscriptStrs import plotTargetGenesGoEnrichmentStr
                # outName = compPair
                # gene2goFile = args.gene2goFile if args.gene2goFile else optionTools.gene2go
                # if not gene2goFile:
                #     print "You don't provide gene2go file, the GO enrichment will not be carried out!"
                #     continue
                #
                # from rpy2 import robjects
                # from rpy2.rinterface import RRuntimeWarning
                # warnings.filterwarnings("ignore", category=RRuntimeWarning)
                # robjects.r(plotTargetGenesGoEnrichmentStr)
                # robjects.r.plotTargetGenesGoEnrichment("dasg.lst", outName, gene2goFile, outName, float(args.cutoff),
                #                                        args.filterBy, int(args.showCategory))
                enrichResult = os.path.abspath(os.path.join(i, "sigDiff.goEnrichResults.txt"))
                enrichPlot = convertPdf2png(inPdf=os.path.abspath(os.path.join(i, "sigDiff.goEnrichResults.pdf")))
                resultDict["das"][compPair].update({"goEnrichResults": enrichResult})
                resultDict["das"][compPair].update({"goEnrichPlot": enrichPlot})

    return resultDict

def generateHtml(dataToProcess, dirSpec, optionTools, args):
    # resolveDir(iflasOut)
    print getCurrentTime() + " Start generating html report..."
    resultDict = retrieveResults(dataToProcess, dirSpec, optionTools, args)
    htmlOut = os.path.join(dirSpec.out_dir, "reports", "html")
    # htmlOut = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/iflas_html_test"
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    resolveDir(htmlOut)
    copy_tree(os.path.join(scriptDir, "assets"), "assets")
    # resultDict = {"allSampleMerged": {"testGene1": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"},
    #               "das": {"sampl1_vs_sample2": {"IR": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/wangbo_embryo_vs_wangbo_endosperm.sigDiffAS/IR.sig.txt",
    #                                             "SE": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/wangbo_embryo_vs_wangbo_endosperm.sigDiffAS/SE.sig.txt",
    #                                             "A3SS": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/wangbo_embryo_vs_wangbo_endosperm.sigDiffAS/A3SS.sig.txt",
    #                                             "A5SS": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/wangbo_embryo_vs_wangbo_endosperm.sigDiffAS/A5SS.sig.txt",
    #                                             "asDistribution": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf",
    #                                             "goEnrichResults": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/wangbo_embryo_vs_wangbo_endosperm.sigDiffAS/A5SS.sig.txt",
    #                                             "goEnrichPlot": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"}},
    #               "sample1": {"basicStatistics": {"readsCorrection": ["readsCorrection", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                               "GC_of_raw_flnc": ["GC_of_raw_flnc", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                               "GC_across_raw_flnc": ["GC_across_raw_flnc", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                               "LengthDistribution": ["LengthDistribution", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                               "juncSupported": ["juncSupported", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                               "asPattern": {"asAnno": ["asAnnotation", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"],
    #                                                             "asSpliceSite": ["asSpliceSite", "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf"]}},
    #                           "gene": {"testGene1": {"isoformStruc": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf",
    #                                                  "alleleAS": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf",
    #                                                  "palenAsPlots": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf",
    #                                                  "palenAPA": "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/iflas_html/test.pdf",}}},
    #               }

    dataList = []
    dataList.append(generateMainPage(resultDict))
    if "allSampleMerged" in resultDict:
        dataList.extend(generateAllSampleMergedPage(resultDict))
    dataList.extend(generateSingleSampleStatisticsPage(resultDict))
    dataList.extend(generateSingleSampleGenePage(resultDict))
    if "das" in resultDict:
        dataList.extend(generateDasPage(resultDict))

    jsonDict = {"code": 0, "data": dataList}
    resolveDir(os.path.join(dirSpec.out_dir, "reports", "json"), chdir=False)
    jsonFile = os.path.join(dirSpec.out_dir, "reports", "json", "search.json")
    printJson(jsonDict, outfile=jsonFile, indent=4)

    print getCurrentTime() + " End generating html report!"

