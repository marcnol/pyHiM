#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:42:25 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    Session, Parameters,Folders,load_json)

from fileProcessing.functionCaller import HiMFunctionCaller



def test_processesPWDmatrices():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= load_json(testDataFileName)
    else:
        raise FileNotFoundError()
        
    root_folder = testData["test_processesPWDmatrices"]["rootFolder"]
    ilabel=testData["test_processesPWDmatrices"]["labels"]
    expectedOutputs = testData["test_processesPWDmatrices"]["expectedOutputs"]

    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)
                
    labels_to_process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    run_parameters={}
    run_parameters["rootFolder"]=root_folder
    run_parameters["parallel"]=False

    him = HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()  
    
    # sets parameters
    current_param = Parameters(run_parameters["rootFolder"], him.labels_to_process[ilabel]["parameterFile"])
    current_param.param_dict['parallel']=him.parallel
    
    him.process_pwd_matrices(current_param, ilabel)

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
            
# test_processesPWDmatrices()