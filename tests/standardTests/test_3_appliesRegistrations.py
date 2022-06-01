#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:08:57 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    Session, Parameters,Folders,load_json)

from fileProcessing.functionCaller import HiMFunctionCaller


def test_appliesProjections():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= load_json(testDataFileName)
    else:
        raise FileNotFoundError()
        
    root_folder = testData["test_appliesRegistrations"]["root_folder"]
    filename_to_process = testData["test_appliesRegistrations"]["filename_to_process"]
    expectedOutputs = testData["test_appliesRegistrations"]["expectedOutput"]
    ilabel=testData["test_appliesRegistrations"]["label"]

    run_parameters={}
    run_parameters["root_folder"]=root_folder
    run_parameters["parallel"]=False
    
    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)
            
    labels_to_process = [
        {"label": "fiducial", "parameter_file": "infoList_fiducial.json"},
        {"label": "barcode", "parameter_file": "infoList_barcode.json"},
        {"label": "DAPI", "parameter_file": "infoList_DAPI.json"},
        {"label": "RNA", "parameter_file": "infoList_RNA.json"},
    ]

    him = HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()  

    # sets parameters
    current_param = Parameters(run_parameters["root_folder"], him.labels_to_process[ilabel]["parameter_file"])
    current_param.param_dict['parallel']=him.parallel
    
    him.apply_registrations(current_param, ilabel)

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
    
# test_appliesProjections()
