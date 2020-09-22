#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:08:51 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    session, log, Parameters,folders,loadJSON)

from imageProcessing.makeProjections import makeProjections, makes2DProjectionsFile


def test_makeProjections():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_makeProjections"]["rootFolder"]
    fileName2Process = testData["test_makeProjections"]["fileName2Process"]
    expectedOutputs = testData["test_makeProjections"]["expectedOutputs"]

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
     
    ilabel=2
    
    # for ilabel in range(len(labels2Process)):
    label = labels2Process[ilabel]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    log1.addSimpleText("**Analyzing label: {}**".format(label))
    
    # sets parameters
    param = Parameters(rootFolder, labelParameterFile)
    param.param['parallel']=False
        
    # [projects 3D images in 2d]
    dataFolder = folders(param.param["rootFolder"])
    dataFolder.createsFolders(rootFolder, param)

    # expectedOutputFileName =  dataFolder.outputFolders["zProject"] + os.sep + os.path.basename(fileName2Process).split(".")[0] + "_2d.npy"

    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)

    
    makes2DProjectionsFile(fileName2Process, param, log1, session1, dataFolder)                    

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
    
    
# test_makeProjections()