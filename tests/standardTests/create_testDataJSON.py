#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:42:53 2020

@author: marcnol
"""

import os, argparse
from fileProcessing.fileManagement import (
    Session, Parameters,Folders,save_json)

#%% defines data folder

parser = argparse.ArgumentParser()
parser.add_argument("-F", "--root_folder", help="Folder with images")
args = parser.parse_args()

if args.root_folder:
    root_folder = args.root_folder
else:
    # root_folder = os.getcwd()
    root_folder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"

#%% constructs dictionary with settings

testData={}

testData['test_makeProjections']={}
testData['test_makeProjections']['root_folder']=root_folder
testData['test_makeProjections']['filename_to_process']= root_folder +os.sep+"scan_006_DAPI_001_ROI_converted_decon_ch00.tif"
testData['test_makeProjections']['expectedOutputs']= [root_folder+"/zProject/scan_006_DAPI_001_ROI_converted_decon_ch00_2d.npy"]

testData['test_alignImages']={}
testData['test_alignImages']['root_folder']=root_folder
testData['test_alignImages']['label']= 0

testData['test_appliesRegistrations']={}
testData['test_appliesRegistrations']['root_folder']=root_folder
testData['test_appliesRegistrations']['filename_to_process']= root_folder+os.sep+"scan_001_RT29_001_ROI_converted_decon_ch00.tif"
testData['test_appliesRegistrations']['expectedOutput']= [root_folder+"/alignImages/scan_001_RT29_001_ROI_converted_decon_ch01_2d_registered.npy"]
testData['test_appliesRegistrations']['label']= 1

testData['test_segmentsMasks']={}
testData['test_segmentsMasks']['root_folder']=root_folder
testData['test_segmentsMasks']['filename_to_process']= [root_folder+os.sep+"scan_001_RT29_001_ROI_converted_decon_ch01.tif",
                                                     root_folder+os.sep+"scan_006_DAPI_001_ROI_converted_decon_ch00.tif"]
testData['test_segmentsMasks']['expectedOutputs']= [root_folder+"/segmentedObjects/segmentedObjects_barcode.dat",
                                                    root_folder+"/segmentedObjects/scan_006_DAPI_001_ROI_converted_decon_ch00_Masks.npy"]
testData['test_segmentsMasks']['labels']= [1, 2]

testData['test_processesPWDmatrices']={}
testData['test_processesPWDmatrices']["root_folder"] = root_folder
testData['test_processesPWDmatrices']['expectedOutputs']= [root_folder+"/buildsPWDmatrix/buildsPWDmatrix_HiMscMatrix.npy"]
testData['test_processesPWDmatrices']['labels']= 2

#%% saves dictionary with settings

file_name="testData.json"
save_json(file_name, testData)

