#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:42:53 2020

@author: marcnol
"""

import os, argparse
from fileProcessing.fileManagement import (
    session, log, Parameters,folders,saveJSON)

#%% defines data folder

parser = argparse.ArgumentParser()
parser.add_argument("-F", "--rootFolder", help="Folder with images")
args = parser.parse_args()

if args.rootFolder:
    rootFolder = args.rootFolder
else:
    # rootFolder = os.getcwd()
    rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"

#%% constructs dictionary with settings

testData={}

testData['test_makeProjections']={}
testData['test_makeProjections']['rootFolder']=rootFolder
testData['test_makeProjections']['fileName2Process']= rootFolder +os.sep+"scan_006_DAPI_001_ROI_converted_decon_ch00.tif"
testData['test_makeProjections']['expectedOutputs']= [rootFolder+"/zProject/scan_006_DAPI_001_ROI_converted_decon_ch00_2d.npy"]

testData['test_alignImages']={}
testData['test_alignImages']['rootFolder']=rootFolder
testData['test_alignImages']['label']= 0

testData['test_appliesRegistrations']={}
testData['test_appliesRegistrations']['rootFolder']=rootFolder
testData['test_appliesRegistrations']['fileName2Process']= rootFolder+os.sep+"scan_001_RT29_001_ROI_converted_decon_ch00.tif"
testData['test_appliesRegistrations']['expectedOutput']= [rootFolder+"/alignImages/scan_001_RT29_001_ROI_converted_decon_ch01_2d_registered.npy"]
testData['test_appliesRegistrations']['label']= 1

testData['test_segmentsMasks']={}
testData['test_segmentsMasks']['rootFolder']=rootFolder
testData['test_segmentsMasks']['fileName2Process']= [rootFolder+os.sep+"scan_001_RT29_001_ROI_converted_decon_ch01.tif",
                                                     rootFolder+os.sep+"scan_006_DAPI_001_ROI_converted_decon_ch00.tif"]
testData['test_segmentsMasks']['expectedOutputs']= [rootFolder+"/segmentedObjects/segmentedObjects_barcode.dat",
                                                    rootFolder+"/segmentedObjects/scan_006_DAPI_001_ROI_converted_decon_ch00_Masks.npy"]
testData['test_segmentsMasks']['labels']= [1, 2]

testData['test_processesPWDmatrices']={}
testData['test_processesPWDmatrices']["rootFolder"] = rootFolder
testData['test_processesPWDmatrices']['expectedOutputs']= [rootFolder+"/buildsPWDmatrix/buildsPWDmatrix_HiMscMatrix.npy"]
testData['test_processesPWDmatrices']['labels']= 2


testData['test_refitBarcodes3D']={}
testData['test_refitBarcodes3D']["rootFolder"] = rootFolder
testData['test_refitBarcodes3D']['expectedOutputs']= [rootFolder+"/segmentedObjects/segmentedObjects_barcode2D.dat"]
testData['test_refitBarcodes3D']['labels']= 1

testData['test_localDriftCorrection']={}
testData['test_localDriftCorrection']["rootFolder"] = rootFolder
testData['test_localDriftCorrection']['expectedOutputs']= [rootFolder+"/alignImages/LocalShiftsViolinPlot_ROI:001.png"]
testData['test_localDriftCorrection']['labels']= 2

#%% saves dictionary with settings

fileName="testData.json"
saveJSON(fileName, testData)

