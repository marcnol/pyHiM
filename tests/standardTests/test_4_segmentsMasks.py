#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:23:44 2020

@author: marcnol
"""


import os
import pytest

from fileProcessing.fileManagement import (
    Session, Parameters,Folders,load_json)

from fileProcessing.functionCaller import HiMFunctionCaller


def test_segmentsMasks():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= load_json(testDataFileName)
    else:
        raise FileNotFoundError()
        
    root_folder = testData["test_segmentsMasks"]["rootFolder"]
    filename_to_process = testData["test_segmentsMasks"]["filename_to_process"]
    expectedOutputs = testData["test_segmentsMasks"]["expectedOutputs"]
    labels=testData["test_segmentsMasks"]["labels"]

    run_parameters={}
    run_parameters["rootFolder"]=root_folder
    run_parameters["parallel"]=False

    him = HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()  
    
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

    for ilabel,file_name in zip(labels,filename_to_process):
        
        # sets parameters
        current_param = Parameters(run_parameters["rootFolder"], him.labels_to_process[ilabel]["parameterFile"])
        current_param.param_dict['parallel']=him.parallel
            
        data_folder = Folders(current_param.param_dict["rootFolder"])
        data_folder.create_folders(root_folder, current_param)
    
        him.segment_masks(current_param, ilabel)      

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)        