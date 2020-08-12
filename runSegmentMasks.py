#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:31:14 2020

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================

# ---- stardist
from __future__ import print_function, unicode_literals, absolute_import, division

import os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    session, log, Parameters, writeString2File)
from imageProcessing.segmentMasks import segmentMasks

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")



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
    sessionName = "segmentMasks"
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

        # [applies registration to DAPI and barcodes]
        if label != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            # [segments DAPI and spot masks]
            if label != "RNA" and param.param["acquisition"]["label"] != "RNA":
                segmentMasks(param, log1, session1,fileName)

        print("\n")
        del param
    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    del log1, session1
    print("Elapsed time: {}".format(datetime.now() - begin_time))

# 