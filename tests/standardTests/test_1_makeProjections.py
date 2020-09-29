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

from fileProcessing.functionCaller import HiMfunctionCaller


def test_makeProjections():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"
    if os.path.exists(testDataFileName):
        testData= loadJSON(testDataFileName)
    else:
        raise FileNotFoundError()
        
    rootFolder = testData["test_makeProjections"]["rootFolder"]
    fileName2Process = testData["test_makeProjections"]["fileName2Process"]
    expectedOutputs = testData["test_makeProjections"]["expectedOutputs"]

    runParameters={}
    runParameters["rootFolder"]=rootFolder
    runParameters["parallel"]=False

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()  
     
    ilabel=2
    # sets parameters
    param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
    param.param['parallel']=HiM.parallel
                
    expectedOutputsTimeStamped={}
    for x in expectedOutputs:
        if os.path.exists(x):
            expectedOutputsTimeStamped[x]=os.path.getmtime(x)

    # [projects 3D images in 2d]
    HiM.makeProjections(param)    

    assert sum([os.path.exists(x) for x in expectedOutputs]) == len(expectedOutputs) 
    
    test=[]
    for key in expectedOutputsTimeStamped.keys():
        if os.path.getmtime(x)>expectedOutputsTimeStamped[x]:
            test.append(True)
            
    assert len(test)==len(expectedOutputs)
    
    
test_makeProjections()