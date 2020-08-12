#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:24:59 2020

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse

from imageProcessing.projectsBarcodes import projectsBarcodes
from fileProcessing.fileManagement import (session, log, Parameters)

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "/home/marcnol/data/Experiment_20/Embryo_1"
        # rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        # rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'

    print("parameters> rootFolder: {}".format(rootFolder))

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
    ]

    # session
    session1 = session(rootFolder, "processingPipeline")

    # setup logs
    log1 = log(rootFolder)

    labelParameterFile = labels2Process[1]["parameterFile"]
    param = Parameters(rootFolder, labelParameterFile)

    projectsBarcodes(param, log1, session1)
