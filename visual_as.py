#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: visual_as.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import glob
import cPickle as pickle
import numpy as np
from commonFuncs import *
from commonObjs import *

class MainSection(object):
    '''
    Set values of variables in main section.
    '''
    def __init__(self, **args):
        self.sectionName = getAttribute("sectionName", "[SinglePlotConfig]", **args)
        self.legend      = getAttribute("legend", True, **args)
        self.width       = getAttribute("width", 10.0, **args)
        self.height      = getAttribute("height", 11.0, **args)
        self.fontsize    = getAttribute("fontsize", 8, **args)
        self.output_file = getAttribute("fout", None, **args)
        self.shrink_introns = getAttribute("shrink_introns", False, **args)

    def printStr(self):
        return "%s\nlegend = %s\nwidth = %.1f\nheight = %.1f\nfontsize = %.1f\noutput_file = %s\nshrink_introns = %s\n" % \
               (self.sectionName, self.legend, self.width, self.height, self.fontsize, self.output_file, self.shrink_introns)


class PlotSection(object):
    '''
    Set values of variables in plot section.
    '''
    def __init__(self, **args):
        self.section_name = getAttribute("section_name", "[testGeneModel]", **args)
        self.plot_type = getAttribute("plot_type", "gene", **args)
        self.relative_size = getAttribute("relative_size", 10.0, **args)
        self.source_file = getAttribute("source_file", None, **args)
        self.title_string = getAttribute("title_string", "testTitle", **args)
        if self.plot_type == "gene":
            self.gene_name = getAttribute("gene_name", "testGeneName", **args)
            self.file_format = getAttribute("file_format", "gene_model", **args)
            self.hide = getAttribute("hide", False, **args)
        if self.plot_type == "junctions":
            self.labels = getAttribute("labels", True, **args)
            self.min_coverage = getAttribute("min_coverage", 5, **args)

    def printStr(self):
        returnStr = "%s\nplot_type = %s\nrelative_size = %.1f\nsource_file = %s\ntitle_string = %s\n" % \
                    (self.section_name, self.plot_type, self.relative_size, self.source_file, self.title_string)
        if self.plot_type == "gene":
            returnStr += "gene_name = %s\nfile_format = %s\nhide = %s\n" % (self.gene_name, self.file_format, self.hide)
        if self.plot_type == "junctions":
            returnStr += "labels = %s\nmin_coverage = %d\n" % (self.labels, self.min_coverage)
        return returnStr

#####################################
def getAttribute(key, default, **args):
    return default if key not in args else args[key]

def getFileRowCounts(inFile):
    return len(open(inFile).readlines())

def resizeTrackRatio(itemCounts):
    originRelativeSize = 10.0
    if itemCounts <= 10:
        return originRelativeSize
    elif 10 <= itemCounts and itemCounts <= 100:
        return 2 * originRelativeSize * np.log10(itemCounts) - originRelativeSize
    else:
        return originRelativeSize * 3


def parallelPlotter(gene, gpeTargetGenePickle, sampleTargetGenePickle, dataObjs, dirSpec, juncDict):
    projectName, sampleName = dataObjs.project_name, dataObjs.sample_name
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    isoViewDir = os.path.join(baseDir, "isoViewer")
    resolveDir(os.path.join(isoViewDir, gene))

    sampleTargetGeneObj = pickle.loads(sampleTargetGenePickle)
    targetGeneChrom = sampleTargetGeneObj.chrom
    targetGeneMinpos = sampleTargetGeneObj.minpos
    targetGeneMaxpos = sampleTargetGeneObj.maxpos
    targetGeneStrand = sampleTargetGeneObj.strand

    # Gene model section
    geneModelGPE = "{}.gpe".format(gene)
    geneModelGTF = "{}.gtf".format(gene)
    geneModelGFF = "{}.gff".format(gene)

    gpeTargetGene = pickle.loads(gpeTargetGenePickle)
    geneModelGPEOut = open(geneModelGPE, "w")
    gpeMinposList, gpeMaxposList = [targetGeneMinpos], [targetGeneMaxpos]
    for i in gpeTargetGene:
        gpeMinposList.append(i.txStart)
        gpeMaxposList.append(i.txEnd)
        print >> geneModelGPEOut, i
    geneModelGPEOut.close()
    plotMinpos = min(gpeMinposList)
    plotMaxpos = max(gpeMaxposList)
    targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
    os.system("genePredToGtf file {} {} -source=iFLAS".format(geneModelGPE, geneModelGTF))
    os.system("gene_model_to_splicegraph.py -m {} -o {} 2>/dev/null".format(geneModelGTF, geneModelGFF))
    geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                    gene_name=gene, relative_size=5.0,
                                    title_string="Gene model for %gene", hide=True)
    geneModelItemCount = getFileRowCounts(geneModelGPE)
    geneModelRelativeSize = resizeTrackRatio(geneModelItemCount)
    geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                       source_file=geneModelGFF, relative_size=geneModelRelativeSize,
                                       title_string="Gene model for %gene [{}({})]".format(
                                           targetGeneRegion, targetGeneStrand))

    # Corrected tgs reads from the pipeline
    postCorrIsoforms = sampleTargetGeneObj.reads
    postCorrIsoGPE = "{}.iFLAS.gpe".format(gene)
    postCorrIsoGTF = "{}.iFLAS.gtf".format(gene)
    postCorrIsoGFF = "{}.iFLAS.gff".format(gene)
    postCorrIsoGPEOut = open(postCorrIsoGPE, "w")
    for read in postCorrIsoforms:
        print >> postCorrIsoGPEOut, postCorrIsoforms[read].go_to_gpe()
    postCorrIsoGPEOut.close()
    os.system("genePredToGtf file {} {} -source=iFLAS".format(postCorrIsoGPE, postCorrIsoGTF))
    postCorrIsoItemCounts = len(postCorrIsoforms)
    postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
    os.system("gene_model_to_splicegraph.py -m {} -o {} -a 2>/dev/null".format(postCorrIsoGTF, postCorrIsoGFF))
    postCorrIsoPlotType = "splice_graph"
    postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse]", plot_type=postCorrIsoPlotType,
                                  source_file=postCorrIsoGFF, relative_size=postCorrIsoRelativeSize,
                                  title_string="Corrected isoforms and AS events in %s from long reads data" % gene)


    figOut = gene + ".pdf"
    cfgOut = open(gene + ".cfg", "w")
    majorItemCount = geneModelItemCount + postCorrIsoItemCounts
    figHeight = 10 if majorItemCount <= 50 else 15 if majorItemCount <= 150 else 20
    mainSec = MainSection(fout=figOut, height=figHeight)
    print >> cfgOut, mainSec.printStr()
    print >> cfgOut, geneModelHidePlot.printStr()
    print >> cfgOut, geneModelVisiblePlot.printStr()
    print >> cfgOut, postCorrIsoPlot.printStr()

    # Abundance from sam file in NGS pipeline
    if dataObjs.ngs_left_reads != None or dataObjs.ngs_right_reads != None:
        ngsSams = []
        for n in range(len(dataObjs.ngs_left_reads.split(";"))):
            repeatName = "repeat" + str(n)
            bamFile = os.path.join(baseDir, "mapping", "rna-seq", "alignment", "{}/{}.sorted.bam".format(repeatName, repeatName))
            targetSam = "{}.{}.sam".format(repeatName, gene)
            cmd = "samtools view -h {} {} > {}".format(bamFile, targetGeneRegion, targetSam)
            subprocess.call(cmd, shell=True)
            targetDepth = "{}.{}.depth".format(repeatName, gene)
            cmd = "sam_to_depths.py {} -o {} 2>/dev/null".format(targetSam, targetDepth)
            subprocess.call(cmd, shell=True)
            if juncDict:
                tmpOut = open("tmp.depth", "w")
                with open(targetDepth) as f:
                    for line in f.readlines():
                        if line.startswith("C") or line.startswith("D"):
                            print >> tmpOut, line.strip("\n")
                        else:
                            infoList = line.strip("\n").split("\t")
                            junc = "{}:{}-{}".format(infoList[1], int(infoList[3]), int(infoList[4])-1)
                            if junc in juncDict:
                                print >> tmpOut, line.strip("\n")
                tmpOut.close()
                os.rename("tmp.depth", targetDepth)
            ngsSams.append(targetSam)
            depthPlot = PlotSection(section_name="[Reads_junc_{}_{}]".format(repeatName, sampleName),
                                    plot_type="junctions", source_file=targetDepth, relative_size=5.0,
                                    min_coverage=10, title_string="%s junction depth in sample %s" % (gene, repeatName))
            readsPlot = PlotSection(section_name="[Reads_%s]" % (repeatName), plot_type="read_depth",
                                    source_file=targetSam, relative_size=5.0,
                                    title_string="%s read coverage in sample %s" % (gene, repeatName))
            print >> cfgOut, readsPlot.printStr()
            print >> cfgOut, depthPlot.printStr()

    cfgOut.close()
    os.system("plotter.py %s.cfg 2>/dev/null" % gene)
    # removeFiles(ngsSams)


def processTargetGenes(targetGenes):
    if os.path.isfile(targetGenes):
        return getDictFromFile(targetGenes)
    else:
        return targetGenes.split(",")


def visual_as_merge(dataToProcess=None, targetGenes=None, refParams=None, dirSpec=None):
    isoViewerDir = os.path.join(dirSpec.out_dir, "isoViewer_sample_merged")
    resolveDir(isoViewerDir)
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in all samples..."
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in all samples..."
        targetGenes = processTargetGenes(targetGenes)
    refGpeObj = GenePredObj(refParams.ref_gpe, False)
    for geneName in targetGenes:
        mergedData = {}
        for dataObj in dataToProcess:
            tmpDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "isoViewer", geneName)
            sampleName = "{}_{}".format(dataObj.project_name, dataObj.sample_name)
            if validateDir(tmpDir):
                tmpDict = {"tgs": [], "ngs": {}}
                isoformGff = os.path.abspath(os.path.join(tmpDir, "{}.iFLAS.gff".format(geneName)))
                isoformGpe = os.path.abspath(os.path.join(tmpDir, "{}.iFLAS.gpe".format(geneName)))
                ngsSamList = glob.glob(os.path.join(tmpDir, "repeat*.{}.sam".format(geneName)))
                if validateFile(isoformGff):
                    tmpDict["tgs"] = [isoformGff, isoformGpe]
                for i in range(len(ngsSamList)):
                    if validateFile(os.path.abspath(ngsSamList[i])):
                        repeat = os.path.basename(ngsSamList[i]).split(".")[0]
                        samFile = os.path.abspath(ngsSamList[i])
                        depthFile = os.path.abspath(os.path.join(tmpDir, "{}.{}.depth".format(repeat, geneName)))
                        tmpDict["ngs"].update({repeat: [samFile, depthFile]})
                mergedData.update({sampleName: tmpDict})
        if mergedData:
            resolveDir(os.path.join(isoViewerDir, geneName))

            isoformGpes = [mergedData[x]["tgs"][1] for x in mergedData]
            minPos = []
            maxPos = []
            for x in isoformGpes:
                gpeObjs = GenePredObj(x, False).geneName2gpeObj[geneName]
                minPos = [gpeObj.txStart for gpeObj in gpeObjs]
                maxPos = [gpeObj.txEnd for gpeObj in gpeObjs]

            refGeneObj = refGpeObj.geneName2gpeObj[geneName]

            targetGeneChrom = refGeneObj[0].chrom
            targetGeneStrand = refGeneObj[0].strand
            plotMinpos = min([x.txStart for x in refGeneObj] + minPos)
            plotMaxpos = max([x.txEnd for x in refGeneObj] + maxPos)

            # Gene model section
            geneModelGPE = "{}.gpe".format(geneName)
            geneModelGTF = "{}.gtf".format(geneName)
            geneModelGFF = "{}.gff".format(geneName)

            geneModelGPEOut = open(geneModelGPE, "w")
            for i in refGeneObj:
                print >> geneModelGPEOut, i
            geneModelGPEOut.close()
            targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
            os.system("genePredToGtf file {} {} -source=iFLAS".format(geneModelGPE, geneModelGTF))
            os.system("gene_model_to_splicegraph.py -m {} -o {} 2>/dev/null".format(geneModelGTF, geneModelGFF))
            sectionToPlot = []
            geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                            gene_name=geneName, relative_size=5.0,
                                            title_string="Gene model for %gene", hide=True)
            sectionToPlot.append(geneModelHidePlot)
            geneModelItemCount = getFileRowCounts(geneModelGPE)
            geneModelRelativeSize = resizeTrackRatio(geneModelItemCount)
            geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                               source_file=geneModelGFF, relative_size=geneModelRelativeSize,
                                               title_string="Gene model for %gene [{}({})]".format(
                                                   targetGeneRegion, targetGeneStrand))
            sectionToPlot.append(geneModelVisiblePlot)
            isoItemCounts = 0
            for sampleName in mergedData:
                tmpSample = mergedData[sampleName]
                postCorrIsoItemCounts = getFileRowCounts(tmpSample["tgs"][1])
                isoItemCounts += postCorrIsoItemCounts
                postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
                postCorrIsoPlotType = "isoforms"
                postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse_{}]".format(sampleName), plot_type=postCorrIsoPlotType,
                                              source_file=tmpSample["tgs"][0], relative_size=postCorrIsoRelativeSize,
                                              title_string="Corrected isoforms and AS events in {} from {} data".format(geneName, sampleName))
                sectionToPlot.append(postCorrIsoPlot)
                for repeatName in tmpSample["ngs"]:
                    ngsDepth = tmpSample["ngs"][repeatName][1]
                    depthPlot = PlotSection(section_name="[Reads_junc_{}_{}]".format(repeatName, sampleName),
                                            plot_type="junctions", source_file=ngsDepth, relative_size=5.0, min_coverage=10,
                                            title_string="{} junction depth in {} {}".format(geneName, sampleName, repeatName))
                    sectionToPlot.append(depthPlot)
                    ngsSam = tmpSample["ngs"][repeatName][0]
                    ngsPlot = PlotSection(section_name="[Reads_{}_{}]".format(repeatName, sampleName), plot_type="read_depth",
                                          source_file=ngsSam, relative_size=5.0,
                                          title_string="{} read coverage in {} {}".format(geneName, sampleName, repeatName))
                    sectionToPlot.append(ngsPlot)
            figOut = geneName + ".pdf"
            cfgOut = open(geneName + ".cfg", "w")
            majorItemCount = geneModelItemCount + isoItemCounts
            figHeight = 20 if majorItemCount <= 50 else 30 if majorItemCount <= 150 else 40
            mainSec = MainSection(fout=figOut, height=figHeight)
            print >> cfgOut, mainSec.printStr()
            for sec in sectionToPlot:
                print >> cfgOut, sec.printStr()
            cfgOut.close()
            os.system("plotter.py {}.cfg 2>/dev/null".format(geneName))
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in all samples done!"
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in all samples done!"


def visual_as(dataObj=None, targetGenes=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in project {} sample {}...".format(projectName, sampleName)
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in project {} sample {}...".format(projectName, sampleName)
        targetGenes = processTargetGenes(targetGenes)

    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    gpeObj = GenePredObj(refParams.ref_gpe, False)
    resolveDir(os.path.join(baseDir, "isoViewer"))
    # tgsIsoFile = os.path.join(baseDir, "refine", "isoformGrouped.bed12+")
    tgsIsoFile = os.path.join(baseDir, "hqIsoforms", "hq.collapsed.bed12+")
    readsDict = {}
    selectedGenes = {}
    with open(tgsIsoFile) as f:
        for line in f:
            readStruc = ReadLineStruc(line)
            if readStruc.chrom not in readsDict:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom] = {readStruc.strand: {readStruc.geneName: gene2reads}}
            elif readStruc.strand not in readsDict[readStruc.chrom]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand] = {readStruc.geneName: gene2reads}
            elif readStruc.geneName not in readsDict[readStruc.chrom][readStruc.strand]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand].update({readStruc.geneName: gene2reads})
            else:
                readsDict[readStruc.chrom][readStruc.strand][readStruc.geneName].update(readStruc)
            if targetGenes == None:
                selectedGenes.update({readStruc.geneName: ""})
            else:
                if readStruc.geneName in targetGenes:
                    selectedGenes.update({readStruc.geneName: ""})

    try:
        juncDict = {}
        dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")
        if os.path.exists(dataObj.ngs_junctions):
            with open(dataObj.ngs_junctions) as f:
                for line in f.readlines():
                    juncDict.update({Junction(line.strip("\n")).jPos: ""})
        pool = Pool(processes=dataObj.single_run_threads)
        for chrom in readsDict:
            for strand in readsDict[chrom]:
                for geneName in readsDict[chrom][strand]:
                    if geneName in selectedGenes:
                        gene2readsObj = readsDict[chrom][strand][geneName]
                        sampleTargetGenePickle = pickle.dumps(gene2readsObj)
                        gpeTargetGenePickle = pickle.dumps(gpeObj.geneName2gpeObj[gene2readsObj.geneName])
                        pool.apply_async(parallelPlotter, (gene2readsObj.geneName, gpeTargetGenePickle,
                                                           sampleTargetGenePickle, dataObj, dirSpec, juncDict))
        pool.close()
        pool.join()
    except Exception as e:
        print e
    os.chdir(prevDir)
    print getCurrentTime() + " Visualize the gene structure compared to the reference genome for project {} sample {} done!".format(projectName, sampleName)
