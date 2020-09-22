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

from imageProcessing.localDriftCorrection import localDriftCorrection

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

    # session
    sessionName = "test_localDriftCorrection"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder=rootFolder,parallel=False)
  
    # for ilabel in range(len(labels2Process)):
    label = labels2Process[ilabel]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    
    # sets parameters
    param = Parameters(rootFolder, labelParameterFile)
    param.param['parallel']=False
        
    dataFolder = folders(param.param["rootFolder"])
    dataFolder.createsFolders(rootFolder, param)

    # [local drift correction]
    if label == "DAPI" and param.param["alignImages"]["localAlignment"]=='overwrite':
        errorCode, _, _ = localDriftCorrection(param, log1, session1)

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)