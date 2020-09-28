#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 20:55:18 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    session, log, Parameters,folders,loadJSON)

from fileProcessing.functionCaller import HiMfunctionCaller

def test_localDriftCorrection():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_localDriftCorrection"]["rootFolder"]
    expectedOutputs = testData["test_localDriftCorrection"]["expectedOutputs"]
    ilabel=testData["test_localDriftCorrection"]["labels"]

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
    runParameters["localAlignment"]=True

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()  
     
    # for ilabel in range(len(labels2Process)):
    label = labels2Process[ilabel]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    
    # sets parameters
    param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
    param.param['parallel']=HiM.parallel
        
    # [local drift correction]
    HiM.localDriftCorrection(param, ilabel)

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        print("{}--{}".format(os.path.getmtime(x),expectedOutputsTimeStamped[x]))
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)