#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol

This file contains routines to process Hi-M datasets

The user needs to provide either a folder by argument or in the source code.
The main() will search for parameter files within the folder provided. All ope-export PATH="$PATH:/home/marcnol/Repositories/pyHiM/"
-ration of the code will be defined in the parameters file.

"""
# =============================================================================
# IMPORTS
# =============================================================================q

import os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import Parameters, log, writeString2File, session

from imageProcessing.alignImages import alignImages, appliesRegistrations
from imageProcessing.makeProjections import makeProjections
from imageProcessing.segmentMasks import segmentMasks
from imageProcessing.localDriftCorrection import localDriftCorrection
from imageProcessing.projectsBarcodes import projectsBarcodes
from matrixOperations.alignBarcodesMasks import processesPWDmatrices

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("--parallel", help="Runs in parallel mode", action="store_true")

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters={}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        runParameters["rootFolder"] = os.getcwd()
       
    if args.parallel:
        runParameters["parallel"] = args.parallel
    else:
        runParameters["parallel"] = False

    return runParameters

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=parseArguments()    

    print("\n--------------------------------------------------------------------------")

    print("parameters> rootFolder: {}".format(runParameters["rootFolder"]))
    now = datetime.now()

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    session1 = session(runParameters["rootFolder"], "processingPipeline")

    # setup logs
    log1 = log(runParameters["rootFolder"])
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format("processingPipeline"))
    if log1.fileNameMD=='.md':
        log1.fileNameMD='HiM_report.md'
        
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M analysis {}".format(now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(runParameters["rootFolder"], labelParameterFile)

        if runParameters["parallel"]:
            param.param['parallel']=True
        else:
            param.param['parallel']=False


        # [projects 3D images in 2d]
        makeProjections(param, log1, session1)

        # [registers fiducials using a barcode as reference]
        if label == "fiducial" and param.param["acquisition"]["label"] == "fiducial":
            log1.report(
                "Making image registrations, ilabel: {}, label: {}".format(ilabel, label), "info",
            )
            alignImages(param, log1, session1)

        # [applies registration to DAPI and barcodes]
        if label != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
            log1.report(
                "Applying image registrations, ilabel: {}, label: {}".format(ilabel, label), "info",
            )
            appliesRegistrations(param, log1, session1)

            # [segments DAPI and spot masks]
            if label != "RNA" and param.param["acquisition"]["label"] != "RNA":
                segmentMasks(param, log1, session1)

        # [2D projects all barcodes in an ROI]
        if label == "barcode":
            projectsBarcodes(param, log1, session1)

        # [refits spots in 3D]

        # [local drift correction]
        if label == "DAPI" and param.param["alignImages"]["localAlignment"]=='overwrite':
            errorCode, _, _ = localDriftCorrection(param, log1, session1)

        # [builds PWD matrix for all folders with images]
        if label == "DAPI":
            processesPWDmatrices(param, log1, session1)

        print("\n")
        del param

    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    del log1, session1
    print("Elapsed time: {}".format(datetime.now() - begin_time))
