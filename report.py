#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: report.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-19
Last modified: 2022-01-19
'''
import os

import pandas as pd
import PyPDF2, glob, itertools, shutil
import warnings
from commonFuncs import *
from commonObjs import *
from multiprocessing import Pool
from rpy2 import robjects
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

def reportReadsCorrectedEval(dataObj=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Reads Corrected Evaluation for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    mappingDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "mapping")
    if not validateDir(mappingDir):
        print getCurrentTime() + " No Corrected Reads available used for evaluation for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return []
    rawMappedBed = os.path.join(mappingDir, "rawFlnc.addCVandID.bed12+")
    correctMappedBed = os.path.join(mappingDir, "flnc.addCVandID.bed12+")
    currDir = os.getcwd()
    basicStatisticsDir = os.path.join(currDir, "basicStatistics")
    resolveDir(basicStatisticsDir)
    from plotRscriptStrs import plotReadsCorrectedEvalStr
    robjects.r(plotReadsCorrectedEvalStr)
    robjects.r.plotReadsCorrectedEval(rawMappedBed, correctMappedBed, "readsCorrectResult.pdf")
    os.chdir(currDir)
    print getCurrentTime() + " Plotting Reads Corrected Evaluation for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)
    return ["readsCorrectResult.pdf"]

def gcAcrossRead(fxFile, outFile, interval=20):
    import math
    from Bio.SeqUtils import GC
    out = open(outFile, "w")
    fileType = validateFaAndFqFile(fxFile)
    for seq in SeqIO.parse(fxFile, fileType):
        if len(seq) < interval:
            chunkSize = 1
        else:
            chunkSize = int(math.ceil(len(seq)/interval))
        gcList = [seq.name]
        for i in xrange(0, len(seq), chunkSize):
            gcList.append(GC(seq[i:i+chunkSize].seq))
        print >>out, "\t".join(map(str, gcList))
    out.close()

def reportReadsContentEval(dataObj=None, refParams=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Reads Content Evaluation for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    if dataObj.data_processed_location and validateFile(dataObj.data_processed_location):
        flncFx = dataObj.data_processed_location
    else:
        flncFx = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "preprocess", dataObj.tgs_plat.lower(), "rawFlnc.fq")
        if not validateFile(flncFx):
            print getCurrentTime() + " No Reads Content can be evaluated for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
            return []

    currDir = os.getcwd()
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    basicStatisticsDir = os.path.join(currDir, "basicStatistics")
    resolveDir(basicStatisticsDir)
    cmd = "seqkit fx2tab -n -g {} | cut -f 1,2 > GC_of_raw_flnc_reads.log".format(flncFx)
    subprocess.call(cmd, shell=True)
    cmd = '''cut -f2 GC_of_raw_flnc_reads.log | {}/distrCurve.R -d -x='Binned GC%' -y='Fraction of reads' -v=50 -w=10 -mainS=28 -p=GC_of_raw_flnc_reads.pdf 2>/dev/null'''.format(utilDir)
    subprocess.call(cmd, shell=True)

    gcAcrossRead(flncFx, "GC_across_raw_flnc_read.log")
    cmd = '''cut -f2- GC_across_raw_flnc_read.log | {}/box.R -stack -nJ -ho=50 -x=Interval -y=GC% -oS=0.5 -w=10 -height=6 -mainS=28 -p=GC_across_raw_flnc_read.pdf 2>/dev/null'''.format(utilDir)
    subprocess.call(cmd, shell=True)

    cmd = "seqkit fx2tab -n -l {} | cut -f 2 > readsLength.lst".format(flncFx)
    subprocess.call(cmd, shell=True)
    cmd = "{}/gpe2bed.pl {} | {}/bedLength.pl | cut -f 13 > refGeneLength.lst".format(utilDir, refParams.ref_gpe, utilDir)
    subprocess.call(cmd, shell=True)
    cmd = '''{}/distrCurves.R -lgPosX=0.8 -lgPosY=0.8 -x1=0 -x2=10000 -d -x='Binned length (limited in 0-10000)' -w=10 *.lst -b=200 -mainS=28 -p=LengthDistribution.curve.pdf 2>/dev/null'''.format(utilDir)
    subprocess.call(cmd, shell=True)
    cmd = '''{}/boxes.R -ng -no -xlab='Category' -ylab='Length of reads or transcripts' *.lst -w=10 -p=LengthDistribution.box.pdf 2>/dev/null'''.format(utilDir)
    subprocess.call(cmd, shell=True)

    dataObj.ngs_junctions = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "mapping", "rna-seq", "reassembly", "junctions.bed")
    if os.path.exists(dataObj.ngs_junctions):
        isoformBed = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "refine", "isoformGrouped.bed12+")
        cmd = '''awk '$10>1' {} | {}/bed2gpe.pl | {}/transSupportByJunction.pl -j {} >supportedByRNAseq.tsv 2>supportedByRNAseq.summary'''.format(
            isoformBed, utilDir, utilDir, dataObj.ngs_junctions)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        cmd = '''awk 'BEGIN{OFS="\t"}{print $1,$2,$4/$3}' supportedByRNAseq.summary | %s/summary2d.R -binWidthX=0.9999999 -binWidthY=0.9999999 -x='Junction count of long reads' -y='Supported Junction Count' -fL=gray -fH=red -mainS=28 -w=10 -p=supportedByRNAseq.pdf 2>/dev/null''' % (utilDir)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(currDir)
    print getCurrentTime() + " Plotting Reads Content Evaluation for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)
    return ["GC_of_raw_flnc_reads.pdf", "GC_across_raw_flnc_read.pdf", "LengthDistribution.curve.pdf", "LengthDistribution.box.pdf"]

def reportASPattern(dataObj=None, dirSpec=None):
    print getCurrentTime() + " Start plotting AS Pattern for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    def _print(asType, annotation, count):
        return "{}\t{}\t{}".format(asType, annotation, count)

    def _getLineCount(myFile, sep="\t", grepStr=None, col=0, uniq=True, header=False):
        if header == False:
            df = pd.read_csv(myFile, sep=sep, header=None)
        else:
            df = pd.read_csv(myFile, sep=sep)
        if grepStr:
            df = df.loc[df.iloc[:, col] == grepStr]
        else:
            df = df.loc[df.iloc[:, col].index]
        return len(df)

    ######################### AS type pattern
    baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
    characAsDir = os.path.join(baseDir, "as_events", "characterization")
    asPatternDir = os.path.join(os.getcwd(), "asPattern")
    resolveDir(asPatternDir)
    if not validateDir(characAsDir):
        print getCurrentTime() + " No AS Pattern file available for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return []
    IR_anno = os.path.join(characAsDir, "IR.anno.lst")
    IR_novel = os.path.join(characAsDir, "IR.novel.lst")
    SE_anno = os.path.join(characAsDir, "SE.anno.lst")
    SE_novel = os.path.join(characAsDir, "SE.novel.lst")
    A3SS_anno = os.path.join(characAsDir, "A3SS.anno.lst")
    A3SS_novel = os.path.join(characAsDir, "A3SS.novel.lst")
    A5SS_anno = os.path.join(characAsDir, "A5SS.anno.lst")
    A5SS_novel = os.path.join(characAsDir, "A5SS.novel.lst")
    APA_anno = os.path.join(characAsDir, "annoPA.tsv")
    APA_novel = os.path.join(characAsDir, "novelPA.tsv")
    PA = os.path.join(characAsDir, "statistic.APA.tsv")

    outStr = "AS_type\tAnnotation\tCount"
    outStr += "\n{}".format(_print("IR", "Anno", _getLineCount(IR_anno)))
    outStr += "\n{}".format(_print("IR", "Novel", _getLineCount(IR_novel)))
    outStr += "\n{}".format(_print("SE", "Anno", _getLineCount(SE_anno)))
    outStr += "\n{}".format(_print("SE", "Novel", _getLineCount(SE_novel)))
    outStr += "\n{}".format(_print("A3SS", "Anno", _getLineCount(A3SS_anno)))
    outStr += "\n{}".format(_print("A3SS", "Novel", _getLineCount(A3SS_novel)))
    outStr += "\n{}".format(_print("A5SS", "Anno", _getLineCount(A5SS_anno)))
    outStr += "\n{}".format(_print("A5SS", "Novel", _getLineCount(A5SS_novel)))
    outStr += "\n{}".format(_print("APA", "Anno", _getLineCount(APA_anno)))
    outStr += "\n{}".format(_print("APA", "Novel", _getLineCount(APA_novel)))
    outStr += "\n{}".format(_print("PA", "Anno", _getLineCount(PA, grepStr="Known", col=2)))
    outStr += "\n{}".format(_print("PA", "Novel", _getLineCount(PA, grepStr="Novel", col=2)))
    asAnnoFile = "{}_{}.AS_annotation.txt".format(dataObj.project_name, dataObj.sample_name)
    out = open(asAnnoFile, "w")
    print >>out, outStr
    out.close()

    ######################### splice site
    irSpliceSite = pd.read_csv(os.path.join(characAsDir, "IR.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    seIncSpliceSite = pd.read_csv(os.path.join(characAsDir, "SE.inc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    seExcSpliceSite = pd.read_csv(os.path.join(characAsDir, "SE.exc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    a3ssIncSpliceSite = pd.read_csv(os.path.join(characAsDir, "A3SS.inc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    a3ssExcSpliceSite = pd.read_csv(os.path.join(characAsDir, "A3SS.exc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    a5ssIncSpliceSite = pd.read_csv(os.path.join(characAsDir, "A5SS.inc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    a5ssExcSpliceSite = pd.read_csv(os.path.join(characAsDir, "A5SS.exc.tmp.txt"), sep="\t", header=None, names=["Dinucleotide", "Frequency"])
    irSpliceSite["AS_type"], seIncSpliceSite["AS_type"], seExcSpliceSite["AS_type"], a3ssIncSpliceSite["AS_type"], \
    a3ssExcSpliceSite["AS_type"], a5ssIncSpliceSite["AS_type"], a5ssExcSpliceSite["AS_type"] = \
        "IR", "SE", "SE", "A3SS", "A3SS", "A5SS", "A5SS"
    irSpliceSite["Category"], seIncSpliceSite["Category"], seExcSpliceSite["Category"], a3ssIncSpliceSite["Category"], \
    a3ssExcSpliceSite["Category"], a5ssIncSpliceSite["Category"], a5ssExcSpliceSite["Category"] = \
        "Inc", "Inc", "Exc", "Inc", "Exc", "Inc", "Exc"
    spliceSite = pd.concat([irSpliceSite, seIncSpliceSite, seExcSpliceSite, a3ssIncSpliceSite, a3ssExcSpliceSite, a5ssIncSpliceSite, a5ssExcSpliceSite])
    asSpliceSiteFile = "{}_{}.AS_spliceSite.txt".format(dataObj.project_name, dataObj.sample_name)
    spliceSite.to_csv(asSpliceSiteFile, sep="\t", header=True, index=False)

    ######################### plot
    from plotRscriptStrs import plotAsCountStatisticsStr
    robjects.r(plotAsCountStatisticsStr)
    asAnnoPdf = "{}_{}.AS_annotation.pdf".format(dataObj.project_name, dataObj.sample_name)
    robjects.r.plotAsCountStatistics(asAnnoFile, asAnnoPdf)

    from plotRscriptStrs import plotAsDinucleotideStatisticsStr
    robjects.r(plotAsDinucleotideStatisticsStr)
    asSpliceSitePdf = "{}_{}.AS_spliceSite.pdf".format(dataObj.project_name, dataObj.sample_name)
    robjects.r.plotAsDinucleotideStatistics(asSpliceSiteFile, asSpliceSitePdf)
    print getCurrentTime() + " Plotting AS Pattern for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)
    return [asAnnoPdf, asSpliceSitePdf]

def reportTargetGeneStructure(dataObj=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Target Gene Structure for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    isoViewerDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "isoViewer")
    if not validateDir(isoViewerDir):
        print getCurrentTime() + " No Target Gene available use for visualization for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return
    allGenePlots = glob.glob("{}/*/*.pdf".format(isoViewerDir))
    writer = PyPDF2.PdfFileWriter()
    for i in allGenePlots:
        pdf = PyPDF2.PdfFileReader(open(i, "rb"))
        for page in range(pdf.getNumPages()):
            writer.addPage(pdf.getPage(page))
    geneStrucDir = os.path.join(os.getcwd(), "geneStrucPlots")
    resolveDir(geneStrucDir)
    output = open("allTargetGeneStructure.pdf", "wb")
    writer.write(output)
    output.close()
    print getCurrentTime() + " Plotting Target Gene Structure for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)

def reportNovelHqAS(dataObj=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Novel High-quality Isoform scores for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    isoformScoreFile = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "isoIdentity", "isoformScore.txt")
    if not validateFile(isoformScoreFile):
        print getCurrentTime() + " No Novel High-quality Isoform scores file available for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return []
    # isoformScore = pd.read_csv(isoformScoreFile, sep="\t", header=None, names=["gene", "isos", "count", "total_count", "freq", "annotation"])
    from plotRscriptStrs import plotNovelHqASStr
    robjects.r(plotNovelHqASStr)
    robjects.r.plotNovelHqAS(isoformScoreFile, "isoformScore.pdf")
    print getCurrentTime() + " Plotting Novel High-quality Isoform scores for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)
    return ["isoformScore.pdf"]

# def reportAlleleAS(dataObj=None, refParams=None, dirSpec=None):
#     print getCurrentTime() + " Start plotting Allele-Specific AS for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
#     baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
#     asHaploFile = os.path.join(baseDir, "alleleAS", "partialAsRelatedHaplotype.txt")
#     if not validateFile(asHaploFile):
#         print getCurrentTime() + " No Allele-Specific AS available for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
#         return
#     asHaplo = pd.read_csv(asHaploFile, header=None, sep="\t", names=["gene", "asType", "haplo1", "haplo1isos", "haplo2", "haplo2isos"])
#     asHaplo = asHaplo.loc[:, ["gene", "haplo1", "haplo1isos", "haplo2", "haplo2isos"]].drop_duplicates()
#     alleleAsDir = os.path.join(os.getcwd(), "alleleAsPlots")
#     scriptDir = os.path.dirname(os.path.abspath(__file__))
#     utilDir = os.path.join(scriptDir, "utils")
#     resolveDir(alleleAsDir)
#     isoformFile = os.path.join(baseDir, "collapse", "isoformGrouped.bed12+")
#     collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
#     flncBam = os.path.join(baseDir, "mapping", "flnc.mm2.sorted.bam")
#     isoBedObj = BedFile(isoformFile, type="bed12+")
#     collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
#     alleleAsPdfs = []
#     for i, row in asHaplo.iterrows():
#         outName = "{}.allele_as".format(row.gene)
#         resolveDir(outName)
#         haplo1isosOut = open("haplo1isosOut.bed", "w")
#         haplo2isosOut = open("haplo2isosOut.bed", "w")
#         chromStarts = []
#         chromEnds = []
#         chrom = set()
#         for x in row.haplo1isos.split("_"):
#             chromStarts.append(isoBedObj.reads[x].chromStart)
#             chromEnds.append(isoBedObj.reads[x].chromEnd)
#             chrom.add(isoBedObj.reads[x].chrom)
#             print >>haplo1isosOut, str(isoBedObj.reads[x])
#         for x in row.haplo2isos.split("_"):
#             chromStarts.append(isoBedObj.reads[x].chromStart)
#             chromEnds.append(isoBedObj.reads[x].chromEnd)
#             chrom.add(isoBedObj.reads[x].chrom)
#             print >>haplo2isosOut, str(isoBedObj.reads[x])
#         haplo1isosOut.close()
#         haplo2isosOut.close()
#         cmd = '''
#             {}/bed2gpe.pl -g 13 haplo1isosOut.bed > haplo1isosOut.gpe;
#             genePredToGtf file haplo1isosOut.gpe haplo1isosOut.gtf;
#             {}/bed2gpe.pl -g 13 haplo2isosOut.bed > haplo2isosOut.gpe;
#             genePredToGtf file haplo2isosOut.gpe haplo2isosOut.gtf;
#         '''.format(utilDir, utilDir)
#         subprocess.call(cmd, shell=True, executable="/bin/bash")
#
#         haplo1reads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in row.haplo1isos.split("_")])
#         haplo2reads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in row.haplo2isos.split("_")])
#         getSubSamByName(flncBam, nameList=haplo1reads, isBam=True, nameListIsFile=False, outPrefix="haplo1.flnc",
#                         sort=True, threads=dataObj.single_run_threads)
#         getSubSamByName(flncBam, nameList=haplo2reads, isBam=True, nameListIsFile=False, outPrefix="haplo2.flnc",
#                         sort=True, threads=dataObj.single_run_threads)
#         cmd = "samtools cat haplo1.flnc.sorted.bam haplo2.flnc.sorted.bam | samtools sort -@ {} > {}.flnc.sorted.bam"
#         cmd = cmd.format(dataObj.single_run_threads, outName)
#         subprocess.call(cmd, shell=True, executable="/bin/bash")
#         refGenome = refParams.ref_genome
#         gtfs = "{},{}".format(os.path.join(os.getcwd(), "haplo1.flnc.sorted.bam"), os.path.join(os.getcwd(), "haplo2.flnc.sorted.bam"))
#         mixedBam = os.path.join(os.getcwd(), "{}.flnc.sorted.bam".format(outName))
#         haploBams = "{},{}".format(os.path.join(os.getcwd(), "haplo1.flnc.sorted.bam"), os.path.join(os.getcwd(), "haplo2.flnc.sorted.bam"))
#
#         targetNgsBam = ""
#         if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
#             ngsBam = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")
#             cmd = "samtools view -h {} {} | samtools sort - > ngsReads.sorted.bam".format(ngsBam, "{}:{}-{}".format(list(chrom)[0], min(chromStarts), max(chromEnds)))
#             subprocess.call(cmd, shell=True, executable="/bin/bash")
#             targetNgsBam = os.path.join(os.getcwd(), "ngsReads.sorted.bam")
#
#         from plotRscriptStrs import plotAlleleAsStructureStr
#         robjects.r(plotAlleleAsStructureStr)
#         robjects.r.plotAlleleAsStructure(refGenome, gtfs, mixedBam, haploBams, targetNgsBam, list(chrom)[0], min(chromStarts), max(chromEnds), outName)
#         alleleAsPdfs.append(os.path.join("{}.pdf".format(outName)))
#         os.chdir(alleleAsDir)
#
#     writer = PyPDF2.PdfFileWriter()
#     for i in alleleAsPdfs:
#         pdf = PyPDF2.PdfFileReader(open(i, "rb"))
#         for page in range(pdf.getNumPages()):
#             writer.addPage(pdf.getPage(page))
#     output = open("alleleAS.pdf", "wb")
#     writer.write(output)
#     output.close()
#     print getCurrentTime() + " Plotting Allelic-Specific AS for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)

def reportAlleleAS1(dataObj=None, refParams=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Allele-Specific AS for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
    asHaploFile = os.path.join(baseDir, "alleleAS", "partialAsRelatedHaplotype.txt")
    if not validateFile(asHaploFile):
        print getCurrentTime() + " No Allele-Specific AS available for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return
    asHaplo = pd.read_csv(asHaploFile, header=None, sep="\t", names=["gene", "asType", "refGene", "asEvent", "haplo1", "haplo1isos", "haplo2", "haplo2isos"])
    alleleAsDir = os.path.join(os.getcwd(), "alleleAsPlots")
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    utilDir = os.path.join(scriptDir, "utils")
    resolveDir(alleleAsDir)
    isoformFile = os.path.join(baseDir, "hqIsoforms", "hq.collapsed.bed12+")
    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    flncBam = os.path.join(baseDir, "mapping", "flnc.mm2.sorted.bam")
    isoBedObj = BedFile(isoformFile, type="bed12+")
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    alleleAsPdfs = []
    highlightColorDict = {"IR": "#0000FF", "SE": "#00EE00", "A3SS": "#FFA500", "A5SS": "#DD77FF"}
    runParams = []
    for name, group in asHaplo.groupby(["gene", "refGene", "asType"]):
        outName = "{}.{}.allele_as".format(name[1], name[2])
        resolveDir(outName)
        haplo1isosOut = open("haplo1isosOut.bed", "w")
        haplo2isosOut = open("haplo2isosOut.bed", "w")
        chromStarts = []
        chromEnds = []
        chrom = set()
        highlightStarts = []
        highlightEnds = []
        haplo1isos = []
        haplo2isos = []
        for i, row in group.iterrows():
            asEvent = row.asEvent
            if name[2] == "SE":
                highlightStarts.append(int(asEvent.split("@")[1].split("-")[0]))
                highlightEnds.append(int(asEvent.split("@")[1].split("-")[-1]))
            else:
                highlightStarts.append(int(asEvent.split(":")[-1].split("-")[0]))
                highlightEnds.append(int(asEvent.split(":")[-1].split("-")[-1]))

            for x in row.haplo1isos.split("_"):
                haplo1isos.append(x)
                chromStarts.append(isoBedObj.reads[x].chromStart)
                chromEnds.append(isoBedObj.reads[x].chromEnd)
                chrom.add(isoBedObj.reads[x].chrom)
                print >> haplo1isosOut, str(isoBedObj.reads[x])
            for x in row.haplo2isos.split("_"):
                haplo2isos.append(x)
                chromStarts.append(isoBedObj.reads[x].chromStart)
                chromEnds.append(isoBedObj.reads[x].chromEnd)
                chrom.add(isoBedObj.reads[x].chrom)
                print >> haplo2isosOut, str(isoBedObj.reads[x])
        haplo1isosOut.close()
        haplo2isosOut.close()
        cmd = '''
            {}/bed2gpe.pl -g 13 haplo1isosOut.bed > haplo1isosOut.gpe;
            genePredToGtf file haplo1isosOut.gpe haplo1isosOut.gtf;
            {}/bed2gpe.pl -g 13 haplo2isosOut.bed > haplo2isosOut.gpe;
            genePredToGtf file haplo2isosOut.gpe haplo2isosOut.gtf;
        '''.format(utilDir, utilDir)
        subprocess.call(cmd, shell=True, executable="/bin/bash")

        haplo1reads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in haplo1isos])
        haplo2reads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in haplo2isos])
        getSubSamByName(flncBam, nameList=haplo1reads, isBam=True, nameListIsFile=False, outPrefix="haplo1.flnc",
                        sort=True, threads=dataObj.single_run_threads)
        getSubSamByName(flncBam, nameList=haplo2reads, isBam=True, nameListIsFile=False, outPrefix="haplo2.flnc",
                        sort=True, threads=dataObj.single_run_threads)
        cmd = "samtools cat haplo1.flnc.sorted.bam haplo2.flnc.sorted.bam | samtools sort -@ {} > {}.flnc.sorted.bam"
        cmd = cmd.format(dataObj.single_run_threads, outName)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        cmd = "samtools index {}.flnc.sorted.bam".format(outName)
        subprocess.call(cmd, shell=True)
        refGenome = refParams.ref_genome
        gtfs = "{},{}".format(os.path.join(os.getcwd(), "haplo1isosOut.gtf"), os.path.join(os.getcwd(), "haplo2isosOut.gtf"))
        mixedBam = os.path.join(os.getcwd(), "{}.flnc.sorted.bam".format(outName))
        haploBams = "{},{}".format(os.path.join(os.getcwd(), "haplo1.flnc.sorted.bam"), os.path.join(os.getcwd(), "haplo2.flnc.sorted.bam"))
        cmd = "samtools index haplo1.flnc.sorted.bam; samtools index haplo2.flnc.sorted.bam"
        subprocess.call(cmd, shell=True)

        targetNgsBam = ""
        if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
            ngsBam = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")
            cmd = "samtools view -h {} {} | samtools sort - > ngsReads.sorted.bam".format(ngsBam, "{}:{}-{}".format(list(chrom)[0], min(chromStarts), max(chromEnds)))
            subprocess.call(cmd, shell=True, executable="/bin/bash")
            targetNgsBam = os.path.join(os.getcwd(), "ngsReads.sorted.bam")
            cmd = "samtools index {}".format(targetNgsBam)
            subprocess.call(cmd, shell=True)

        params = [refGenome, gtfs, mixedBam, haploBams, targetNgsBam, list(chrom)[0], min(chromStarts), max(chromEnds),
                  ",".join(map(str, highlightStarts)), ",".join(map(str, highlightEnds)), highlightColorDict[name[2]],
                  outName, alleleAsDir]
        runParams.append(params)
        os.chdir(alleleAsDir)

    from plotRscriptStrs import plotAlleleAsStructureStr
    robjects.r(plotAlleleAsStructureStr)
    pool = Pool(processes=dataObj.single_run_threads)
    for params in runParams:
        pool.apply_async(robjects.r.plotAlleleAsStructure, (params[0], params[1], params[2], params[3], params[4],
                                                             params[5], params[6], params[7], params[8], params[9],
                                                             params[10], os.path.join(params[12], params[11])))
        alleleAsPdfs.append(os.path.join(params[12], "{}.pdf".format(params[11])))
    pool.close()
    pool.join()

    writer = PyPDF2.PdfFileWriter()
    for i in alleleAsPdfs:
        pdf = PyPDF2.PdfFileReader(open(i, "rb"))
        for page in range(pdf.getNumPages()):
            writer.addPage(pdf.getPage(page))
    output = open("alleleAS.pdf", "wb")
    writer.write(output)
    output.close()
    print getCurrentTime() + " Plotting Allele-Specific AS for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)

def reportPaTailAS(dataObj=None, refParams=None, dirSpec=None):
    print getCurrentTime() + " Start plotting Differential Poly(A) tail length related AS/APA for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
    palenASDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "palenAS")
    if not validateDir(palenASDir):
        print getCurrentTime() + " No Diffential Poly(A) tail length related AS available for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)
        return
    asType = ["IR", "SE", "A3SS", "A5SS"]
    palenAsFiles = []
    for i in asType:
        sigFile = os.path.join(palenASDir, "mergeByJunc", "{}.palenAndAS.sig.bed12+".format(i))
        palenAsFiles.append(sigFile)

    curDir = os.getcwd()
    palenAsPlotsDir = os.path.join(curDir, "palenAsPlots")
    if os.path.exists(palenAsPlotsDir):
        shutil.rmtree(palenAsPlotsDir)
    resolveDir(palenAsPlotsDir)
    cmd = "cat {} | sort -u > palenAndAS.sig.bed12+".format(" ".join(palenAsFiles))
    subprocess.call(cmd, shell=True)
    cmd = "rm *.pdf"
    subprocess.call(cmd, shell=True)
    from plotRscriptStrs import plotPaTailASStr
    robjects.r(plotPaTailASStr)
    robjects.r.plotPaTailAsStructure("palenAndAS.sig.bed12+")

    ##########################
    palenApaPlotsDir = os.path.join(curDir, "palenApaPlots")
    if os.path.exists(palenApaPlotsDir):
        shutil.rmtree(palenApaPlotsDir)
    resolveDir(palenApaPlotsDir)
    makeTxDbGtfStr = '''
        options(ucscChromosomeNames=FALSE)
        refGtfDb <- makeTxDbFromGFF("%s", format = "gtf")
        refGtfTrack <- GeneRegionTrack(refGtfDb, name="Gene model", transcriptAnnotation = "transcript", stackHeight = 0.5)
    ''' % (refParams.ref_gtf)
    robjects.r(makeTxDbGtfStr)
    readBedFile = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "refine", "reads.assigned.unambi.bed12+")
    palenApaFile = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "palenAS", "palenAPA", "apaRelatedPalen.txt")
    alignBam = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "mapping", "flnc.mm2.sorted.bam")
    cmd = "samtools index {}".format(alignBam)
    subprocess.call(cmd, shell=True)
    genePosInfo = BedFile(readBedFile, type="bed12+").getGenePos(bedType="bed12+", geneCol=15)

    runParams = []
    with open(palenApaFile) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            geneName = infoList[0]
            strand = infoList[1]
            paSite = infoList[2]
            apaSupport = infoList[3]
            paLenList = infoList[4]
            chrom, chromStart, chromEnd = genePosInfo[geneName]
            chromStart, chromEnd = int(chromStart), int(chromEnd)
            params = [geneName, paSite, paLenList, apaSupport, alignBam, chrom, chromStart, chromEnd]
            runParams.append(params)

    # from plotRscriptStrs import plotAlleleAsStructureStr
    # robjects.r(plotAlleleAsStructureStr)
    pool = Pool(processes=dataObj.single_run_threads)
    for params in runParams:
        pool.apply_async(robjects.r.plotPaTailApaStructure, (params[0], params[1], params[2], params[3], params[4],
                                                             params[5], params[6], params[7]))
    pool.close()
    pool.join()

    print getCurrentTime() + " Plotting Differential poly(A) tail length related AS/APA for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)

def reportDiffAS(dirSpec=None):
    print getCurrentTime() + " Start plotting Differential-related AS summary..."
    reportDir = os.path.join(dirSpec.out_dir, "das")
    resolveDir(reportDir)
    if len(glob.glob(os.path.join(dirSpec.out_dir, "das", "*.sigDiffAS"))) == 0:
        print getCurrentTime() + " No Differential-related AS summary can be plotted."
        return
    for i in glob.glob(os.path.join(dirSpec.out_dir, "das", "*.sigDiffAS")):
        os.chdir(i)
        irSig = pd.read_csv("{}/IR.sig.txt".format(i), sep="\t")
        seSig = pd.read_csv("{}/SE.sig.txt".format(i), sep="\t")
        a3ssSig = pd.read_csv("{}/A3SS.sig.txt".format(i), sep="\t")
        a5ssSig = pd.read_csv("{}/A5SS.sig.txt".format(i), sep="\t")
        sigDict = {"IR": ["IR", len(irSig)], "SE": ["SE", len(seSig)], "A3SS": ["A3SS", len(a3ssSig)], "A5SS": ["A5SS", len(a5ssSig)]}
        sigDf = pd.DataFrame.from_dict(sigDict, orient="index", columns=["AS_type", "Count"])
        sigDf.to_csv("sigDiff.AS.txt", sep="\t", index=False, header=True)
        from plotRscriptStrs import plotDiffASStr
        robjects.r(plotDiffASStr)
        robjects.r.plotDiffAS("sigDiff.AS.txt", "sigDiff.AS_distribution.pdf")
    print getCurrentTime() + " Plotting Differential-related AS summary done!"

def mergeAllPlots(plot2merge, outPdf):
    print getCurrentTime() + " Start merging plots..."
    if len(plot2merge) == 0:
        print getCurrentTime() + " No plot can be merged!"
        return
    writer = PyPDF2.PdfFileWriter()
    for i in plot2merge:
        pdf = PyPDF2.PdfFileReader(open(i, "rb"))
        for page in range(pdf.getNumPages()):
            writer.addPage(pdf.getPage(page))
    output = open(outPdf, "wb")
    writer.write(output)
    output.close()
    print getCurrentTime() + " Merging plots done!"

def report(dataToProcess=None, refInfoParams=None, dirSpec=None, optionTools=None, args=None):
    reportDir = os.path.join(dirSpec.out_dir, "reports")
    for dataObj in dataToProcess:
        refParams = refInfoParams[dataObj.ref_strain]
        projectName, sampleName = dataObj.project_name, dataObj.sample_name
        print getCurrentTime() + " Generate plot report for project {} sample {}...".format(projectName, sampleName)
        prevDir = os.getcwd()
        subDir = os.path.join(reportDir, "{}_{}".format(projectName, sampleName))
        resolveDir(subDir)

        plots2merge = []
        if args.basic or args.all:
            if dataObj.use_fmlrc2:
                plots2merge.extend(reportReadsCorrectedEval(dataObj=dataObj, dirSpec=dirSpec))
            plots2merge.extend(reportReadsContentEval(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec))
            os.chdir(subDir)
        if args.asp or args.all:
            plots2merge.extend(reportASPattern(dataObj=dataObj, dirSpec=dirSpec))
            os.chdir(subDir)
        if args.geneStruc or args.all:
            reportTargetGeneStructure(dataObj=dataObj, dirSpec=dirSpec)
            os.chdir(subDir)
        # if args.novelIso or args.all:
        #     plots2merge.extend(reportNovelHqAS(dataObj=dataObj, dirSpec=dirSpec))
        if args.asas or args.all:
            reportAlleleAS1(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            os.chdir(subDir)
        if args.palen or args.all:
            reportPaTailAS(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            os.chdir(subDir)

        os.chdir(prevDir)
        print getCurrentTime() + " Generate plot report for project {} sample {} done!".format(projectName, sampleName)
    if args.diff or args.all:
        reportDiffAS(dirSpec=dirSpec)
        os.chdir(reportDir)
    # reportTargetGenesGoEnrichment(dataObj=dataObj, dirSpec=dirSpec)
    if args.html or args.all:
        from generateHtml import generateHtml
        generateHtml(dataToProcess, dirSpec, optionTools, args)

