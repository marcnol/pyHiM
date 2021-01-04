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

from fileProcessing.functionCaller import HiMfunctionCaller, HiM_parseArguments

def test_alignImages():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"

    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()

    rootFolder = testData["test_alignImages"]["rootFolder"]
    ilabel=testData["test_alignImages"]["label"]

    runParameters={}
    runParameters["rootFolder"]=rootFolder
    runParameters["parallel"]=False

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    HiM = HiMfunctionCaller(runParameters, sessionName="makesProjections")
    HiM.initialize()

    # sets parameters
    param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
    param.param['parallel']=HiM.parallel

    HiM.alignImages(param, ilabel)

    dataFolder = folders(param.param["rootFolder"])
    dataFolder.createsFolders(rootFolder, param)

    # dictShifts = loadJSON(dataFolder.outputFiles["dictShifts"] + ".json")
    dictFileName = os.path.splitext(dataFolder.outputFiles["dictShifts"])[0] + ".json"
    dictShifts = loadJSON(dictFileName)

    assert ("ROI:001" in dictShifts and
            "DAPI" in dictShifts["ROI:001"] and
            "RT29" in dictShifts["ROI:001"] and
            "RT37" in dictShifts["ROI:001"])

    # assert (dictShifts['ROI:001']['DAPI']==[-14.36, 45.29] and
    #     dictShifts['ROI:001']['RT29']==[0.28, -2.82] and
    #     dictShifts['ROI:001']['RT37']==[1.58, 11.59])

