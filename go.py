#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: go.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-19
Last modified: 2022-01-19
'''


import warnings
from rpy2 import robjects
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)
from commonFuncs import *


def go(args, optionTools=None, dirSpec=None):
    print getCurrentTime() + " Perform GO enrichment for target genes..."
    if args.gene2goFile == None or optionTools.gene2go == None:
        raise Exception("You should provide gene2go file!")
    gene2goFile = args.gene2goFile if args.gene2goFile else optionTools.gene2go
    if not validateFile(gene2goFile):
        raise Exception("The gene2go file you provide doesn't exist, please check it!")
    targetGeneFile = args.targetGeneFile
    sampleName = args.sampleName
    # scriptDir = os.path.dirname(os.path.realpath(__file__))
    if not validateFile(targetGeneFile):
        validateIndex = []
        if "," in targetGeneFile:
            fileList = targetGeneFile.split(",")
            sampleList = sampleName.split(",")
            for x in range(len(fileList)):
                validateIndex.append(x)
            if len(validateIndex) != 0:
                validateTargetGeneFiles = [os.path.abspath(fileList[x]) for x in validateIndex]
                validateSampleNames = [sampleList[x] for x in validateIndex]
            else:
                raise Exception("The path of gene file you input seems not correct, please check it!")
        else:
            raise Exception("Please input the correct gene file you want to perform GO enrichment!")
    else:
        validateTargetGeneFiles = os.path.abspath(targetGeneFile)
        validateSampleNames = sampleName
    gene2goFile = os.path.abspath(gene2goFile)
    goDir = os.path.join(dirSpec.out_dir, "GO")
    resolveDir(goDir)
    from plotRscriptStrs import plotTargetGenesGoEnrichmentStr
    outName = args.outName
    robjects.r(plotTargetGenesGoEnrichmentStr)
    if isinstance(validateTargetGeneFiles, list) and isinstance(validateSampleNames, list):
        robjects.r.plotTargetGenesGoEnrichment(",".join(validateTargetGeneFiles), ",".join(validateSampleNames),
                                               gene2goFile, outName, float(args.cutoff), args.filterBy, int(args.showCategory))
    else:
        robjects.r.plotTargetGenesGoEnrichment(validateTargetGeneFiles, validateSampleNames, gene2goFile, outName,
                                               float(args.cutoff), args.filterBy, int(args.showCategory))
    print getCurrentTime() + " Perform GO enrichment for target genes done!"
