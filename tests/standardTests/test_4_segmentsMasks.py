#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:23:44 2020

@author: marcnol
"""


import os
import pytest

from fileProcessing.fileManagement import (
    session, log, Parameters,folders,loadJSON)

from fileProcessing.functionCaller import HiMfunctionCaller


def test_segmentsMasks():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_segmentsMasks"]["rootFolder"]
    fileName2Process = testData["test_segmentsMasks"]["fileName2Process"]
    expectedOutputs = testData["test_segmentsMasks"]["expectedOutputs"]
    labels=testData["test_segmentsMasks"]["labels"]

    runParameters={}
    runParameters["rootFolder"]=rootFolder
    runParameters["parallel"]=False

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()  
    
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

    for ilabel,fileName in zip(labels,fileName2Process):
        
        # sets parameters
        param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
        param.param['parallel']=HiM.parallel
            
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.createsFolders(rootFolder, param)
    
        HiM.segmentMasks(param, ilabel)      

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)        