#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:19:10 2020

@author: marcnol
"""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    session, writeString2File, log, Parameters)

from imageProcessing.makeProjections import makeProjections

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")


# =============================================================================
# Local functions
# =============================================================================

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-x", "--fileName", nargs='+', help="fileName to analyze")
    parser.add_argument("--parallel", help="Runs in parallel mode", action="store_true")

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters={}
    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        runParameters["rootFolder"] = os.getcwd()

    if args.fileName:
        runParameters["fileName"] = args.fileName
    else:
        runParameters["fileName"] = None
        
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
    
    print("parameters> rootFolder: {}".format(runParameters["rootFolder"]))
    now = datetime.now()

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName = "makesProjections"
    session1 = session(runParameters["rootFolder"], sessionName)

    # setup logs
    log1 = log(rootFolder=runParameters["rootFolder"],parallel=runParameters["parallel"])
    log1.addSimpleText("\n-------------------------{}-------------------------\n".format(sessionName))
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
        makeProjections(param, log1, session1, runParameters["fileName"] )

        print("\n")
        del param
        
    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    del log1, session1
    print("Elapsed time: {}".format(datetime.now() - begin_time))
