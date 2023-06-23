#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identify_as.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from commonFuncs import *

def identify_as(dataObj=None, refParams=None, dirSpec=None, hqIsoParams=None, args=None):
    refineDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "refine")
    hqIsoDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "hqIsoforms")
    if not validateDir(refineDir):
        from refine import refineJunc
        refineJunc(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, refine=False, adjust=False)

    if not validateDir(hqIsoDir):
        from iso_pu import iso_pu1
        iso_pu1(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, hqIsoParams=hqIsoParams)

    from find_as import find_all_as
    find_all_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
    from find_pa import find_pa
    find_pa(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, confidentPa=args.confidentPa,
            filterPaByRPKM=args.paRPKM, filterPaByCount=args.pa_support)
    from charaterize_as import charaterize_as
    charaterize_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
