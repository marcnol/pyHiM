#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:33:49 2020

@author: marcnol
"""


import os
import pytest

from fileProcessing.fileManagement import (
    Session, load_json,Parameters,Folders)

from fileProcessing.functionCaller import HiMFunctionCaller, him_parse_arguments

def test_alignImages():

    testDataFileName=os.getcwd()+os.sep+"tests"+os.sep+"standardTests"+os.sep+"testData.json"

    if os.path.exists(testDataFileName):
        testData= load_json(testDataFileName)
    else:
        raise FileNotFoundError()

    root_folder = testData["test_alignImages"]["rootFolder"]
    ilabel=testData["test_alignImages"]["label"]

    run_parameters={}
    run_parameters["rootFolder"]=root_folder
    run_parameters["parallel"]=False

    labels_to_process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    him = HiMFunctionCaller(run_parameters, session_name="makesProjections")
    him.initialize()

    # sets parameters
    current_param = Parameters(run_parameters["rootFolder"], him.labels_to_process[ilabel]["parameterFile"])
    current_param.param_dict['parallel']=him.parallel

    him.align_images(current_param, ilabel)

    data_folder = Folders(current_param.param_dict["rootFolder"])
    data_folder.create_folders(root_folder, current_param)

    # dict_shifts = load_json(data_folder.output_files["dictShifts"] + ".json")
    dict_filename = os.path.splitext(data_folder.output_files["dictShifts"])[0] + ".json"
    dict_shifts = load_json(dict_filename)

    assert ("ROI:001" in dict_shifts and
            "DAPI" in dict_shifts["ROI:001"] and
            "RT29" in dict_shifts["ROI:001"] and
            "RT37" in dict_shifts["ROI:001"])

    # assert (dict_shifts['ROI:001']['DAPI']==[-14.36, 45.29] and
    #     dict_shifts['ROI:001']['RT29']==[0.28, -2.82] and
    #     dict_shifts['ROI:001']['RT37']==[1.58, 11.59])

