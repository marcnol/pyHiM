#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:50:55 2020

@author: marcnol
"""

import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    session, writeString2File, log, Parameters)


from imageProcessing.refitBarcodes3D import refitBarcodes3D

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
        rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0"

    print("parameters> rootFolder: {}".format(rootFolder))

    now = datetime.now()
    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName = "refitBarcodes3D"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format("processingPipeline"))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M {}: {}".format(sessionName, now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # [builds PWD matrix for all folders with images]
        if label == "barcode":
            refitBarcodes3D(param, log1, session1)
