#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: getPaSiteDist.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2022-03-12 16:06:59
Last modified: 2022-03-12 16:06:59
'''

in_file = "/data/CrazyHsu_data/Projects/isoseq/iFLAS_20220119/ont_output/B73_ont/B73_all/as_events/pa/paCluster.bed8+"

def generateTable(myFile):
    myDict = {}
    with open(myFile) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            reads = infoList[3]
            readsCount = int(infoList[4])
            paSite = int(infoList[6])
            paSiteName = "{}_{}".format(infoList[0], infoList[7])
            gene = infoList[8]
            if gene not in myDict:
                myDict[gene] = {paSiteName: [paSiteName, paSite, reads, readsCount]}
            else:
                myDict[gene].update({paSiteName: [paSiteName, paSite, reads, readsCount]})

    newDict = {}
    for g in myDict:
        for name in myDict[g]:
            if g not in newDict:
                newDict[g] = {"exc": myDict[g][name], "inc": myDict[g][name]}
            else:
                if myDict[g][name][1] < newDict[g]["exc"][1] and myDict[g][name][3] > newDict[g]["exc"][3]:
                    newDict[g]["exc"] = myDict[g][name]
                if myDict[g][name][1] > newDict[g]["inc"][1] and myDict[g][name][3] > newDict[g]["inc"][3]:
                    newDict[g]["inc"] = myDict[g][name]
        if newDict[g]["exc"][0] == newDict[g]["inc"][0]:
            newDict.pop(g, None)

    for g in newDict:
        for d in newDict[g]:
            for r in newDict[g][d][2].split(","):
                print "\t".join([g, newDict[g][d][0], d, r])

generateTable(in_file)