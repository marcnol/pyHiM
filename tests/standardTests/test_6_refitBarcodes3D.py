#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 20:28:01 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    session, log, Parameters,folders,loadJSON)

from imageProcessing.refitBarcodes3D import refitBarcodesClass

def test_refitBarcodes3D():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_refitBarcodes3D"]["rootFolder"]
    expectedOutputs = testData["test_refitBarcodes3D"]["expectedOutputs"]
    ilabel=testData["test_refitBarcodes3D"]["labels"]


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
    sessionName = "makesProjections"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder=rootFolder,parallel=False)
  
    # for ilabel in range(len(labels2Process)):
    label = labels2Process[ilabel]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    log1.addSimpleText("**Analyzing label: {}**".format(label))
    
    # sets parameters
    param = Parameters(rootFolder, labelParameterFile)
    param.param['parallel']=False
        
    dataFolder = folders(param.param["rootFolder"])
    dataFolder.createsFolders(rootFolder, param)

    # [refits spots in 3D]
    if label == "barcode":
        fittingSession = refitBarcodesClass(param, log1, session1)
        fittingSession.refitFolders()        

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        # print("{}--{}".format(os.path.getmtime(x),expectedOutputsTimeStamped[x]))
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
            
    
# test_refitBarcodes3D()