#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:42:25 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    session, log, Parameters,folders,loadJSON)

from fileProcessing.functionCaller import HiMfunctionCaller



def test_processesPWDmatrices():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_processesPWDmatrices"]["rootFolder"]
    ilabel=testData["test_processesPWDmatrices"]["labels"]
    expectedOutputs = testData["test_processesPWDmatrices"]["expectedOutputs"]

    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)
                
    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    runParameters={}
    runParameters["rootFolder"]=rootFolder
    runParameters["parallel"]=False

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()  
    
    # sets parameters
    param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
    param.param['parallel']=HiM.parallel
    
    HiM.processesPWDmatrices(param, ilabel)

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
            
# test_processesPWDmatrices()