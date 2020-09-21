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

from imageProcessing.segmentMasks import segmentMasks


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

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName = "makesProjections"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder=rootFolder,parallel=False)
     
    for ilabel,fileName in zip(labels,fileName2Process):
        
        # for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))
        
        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)
        param.param['parallel']=False
            
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.createsFolders(rootFolder, param)
    
        segmentMasks(param, log1, session1,fileName=fileName)
        

    if sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs):
        assert True
    else:
        assert False
