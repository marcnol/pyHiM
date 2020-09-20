#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:33:49 2020

@author: marcnol
"""


import os
import pytest

from fileProcessing.fileManagement import (
    session, log, loadJSON,Parameters,folders)

from imageProcessing.alignImages import alignImages

def test_alignImages():
    
    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"

    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_alignImages"]["rootFolder"]
    ilabel=testData["test_alignImages"]["label"]
   
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
        
    if label == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
        alignImages(param, log1, session1)

    dataFolder = folders(param.param["rootFolder"])
    dataFolder.createsFolders(rootFolder, param)
    dictShifts = loadJSON(dataFolder.outputFiles["dictShifts"] + ".json")

    if (dictShifts['ROI:001']['DAPI']==[-13.53, 46.0] and 
        dictShifts['ROI:001']['RT29']==[0.38, -2.98] and 
        dictShifts['ROI:001']['RT37']==[1.7, 11.37]):
        assert True
    else:
        assert False    
    
