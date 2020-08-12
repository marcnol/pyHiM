#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:35:41 2020

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================
import os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    writeString2File,
    session, log, Parameters
    )

from imageProcessing.alignImages import alignImages, appliesRegistrations

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-x", "--fileName", help="fileName to analyze")
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = os.getcwd()

    if args.fileName:
        fileName = args.fileName
    else:
        fileName = None
        
    print("parameters> rootFolder: {}".format(rootFolder))
    now = datetime.now()

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName = "registersImages"
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
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # [registers fiducials using a barcode as reference]
        if label == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
            log1.report(
                "Making image registrations, ilabel: {}, label: {}".format(ilabel, label), "info",
            )
            alignImages(param, log1, session1,fileName)

        # [applies registration to DAPI and barcodes]
        if label != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            log1.report(
                "Applying image registrations, ilabel: {}, label: {}".format(ilabel, label), "info",
            )
            appliesRegistrations(param, log1, session1,fileName)

        print("\n")
        del param
    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    del log1, session1
    print("Elapsed time: {}".format(datetime.now() - begin_time))


