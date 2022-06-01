#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:08:51 2020

@author: marcnol
"""

import os
import pytest

from fileProcessing.fileManagement import (
    Session, Parameters,Folders,load_json)

from fileProcessing.functionCaller import HiMFunctionCaller


def test_makeProjections():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= load_json(testDataFileName)
    else:
        raise FileNotFoundError()
        
    root_folder = testData["test_makeProjections"]["rootFolder"]
    filename_to_process = testData["test_makeProjections"]["filename_to_process"]
    expectedOutputs = testData["test_makeProjections"]["expectedOutputs"]

    run_parameters={}
    run_parameters["rootFolder"]=root_folder
    run_parameters["parallel"]=False

    him = HiMFunctionCaller(run_parameters, session_name="HiM_analysis")
    him.initialize()  
     
    ilabel=2
    # sets parameters
    current_param = Parameters(run_parameters["rootFolder"], him.labels_to_process[ilabel]["parameterFile"])
    current_param.param_dict['parallel']=him.parallel
                
    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)

    # [projects 3D images in 2d]
    him.make_projections(current_param)    

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
    
    
test_makeProjections()