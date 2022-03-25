#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: collapse.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from commonFuncs import *

def collapse(dataObj=None, collapseParams=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Collapse with cDNA_cupcake for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)

    resolveDir(os.path.join(baseDir, "collapse"))
    logDir = os.path.join(baseDir, "log")
    resolveDir(logDir, chdir=False)
    processedFa = os.path.join(baseDir, "mapping", "flnc.processed.fa")
    flncSam = os.path.join(baseDir, "mapping", "flnc.mm2.sam")
    if collapseParams.dun_merge_5_shorter:
        cmd = "collapse_isoforms_by_sam.py --input {} -s {} --max_5_diff {} --max_3_diff {} " \
              "--flnc_coverage {} -c {} -i {} --max_fuzzy_junction {} --dun-merge-5-shorter -o tofu 1>{}/tofu.collapse.log 2>&1"
    else:
        cmd = "collapse_isoforms_by_sam.py --input {} -s {} --max_5_diff {} --max_3_diff {} " \
              "--flnc_coverage {} -c {} -i {} --max_fuzzy_junction {} -o tofu 1>{}/tofu.collapse.log 2>&1"
    cmd = cmd.format(processedFa, flncSam, collapseParams.max_5_diff, collapseParams.max_3_diff, collapseParams.fl_coverage,
                     collapseParams.min_coverage, collapseParams.min_identity, collapseParams.max_fuzzy_junction, logDir)
    subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    print getCurrentTime() + " Collapse with cDNA_cupcake for project {} sample {} done!".format(projectName, sampleName)
