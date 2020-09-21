#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:50:55 2020

@author: marcnol
"""

import argparse
from datetime import datetime

from fileProcessing.fileManagement import (
    session, writeString2File, log, Parameters, daskCluster)

from dask.distributed import Client, LocalCluster, get_client, as_completed, fire_and_forget

from imageProcessing.refitBarcodes3D import refitBarcodesClass

# =============================================================================
# MAIN
# =============================================================================

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
        runParameters["rootFolder"] = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0"

    if args.parallel:
        runParameters["parallel"] = args.parallel
    else:
        runParameters["parallel"] = False
        
    return runParameters


if __name__ == "__main__":

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
    sessionName = "refitBarcodes3D"
    session1 = session(runParameters["rootFolder"], sessionName)

    # setup logs
    log1 = log(runParameters["rootFolder"],parallel=runParameters["parallel"])
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format("processingPipeline"))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M {}: {}".format(sessionName, now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

      
    ilabel=1
    label = labels2Process[1]["label"]
    labelParameterFile = labels2Process[ilabel]["parameterFile"]
    log1.addSimpleText("**Analyzing label: {}**".format(label))

    # sets parameters
    param = Parameters(runParameters["rootFolder"], labelParameterFile)
 
    if not runParameters['parallel']:
        # [builds PWD matrix for all folders with images]
        fittingSession = refitBarcodesClass(param, log1, session1,parallel=runParameters['parallel'])
        fittingSession.refitFolders()
    else:

        daskClusterInstance = daskCluster(20)
        print("Go to http://localhost:8787/status for information on progress...")
        fittingSession = refitBarcodesClass(param, log1, session1,parallel=runParameters['parallel'])
        with LocalCluster(n_workers=daskClusterInstance.nThreads,
                            # processes=True,
                            # threads_per_worker=1,
                            # memory_limit='2GB',
                            # ip='tcp://localhost:8787',
                            ) as cluster, Client(cluster) as client:
    
            result = client.submit(fittingSession.refitFolders)
        
            a = client.gather(result)
    


    # # [builds PWD matrix for all folders with images]
    # fittingSession = refitBarcodesClass(param, log1, session1,parallel=runParameters['parallel'])
    # fittingSession.refitFolders()

    print("Elapsed time: {}".format(datetime.now() - now))
    
        