#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:31:24 2020

@author: marcnol
"""

import os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    folders,
    isnotebook,
    session,
    writeString2File,
    Parameters, 
    log)

from matrixOperations.alignBarcodesMasks import processesPWDmatrices

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = os.getcwd()
        # rootFolder = "."
        # rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug"

    print("parameters> rootFolder: {}".format(rootFolder))
    now = datetime.now()

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName="alignBarcodesMasks"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(sessionName))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M analysis {}".format(now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # [builds PWD matrix for all folders with images]
        if label == "DAPI":
            log1.addSimpleText("**Analyzing label: {}**".format(label))
            processesPWDmatrices(param, log1, session1)
            
