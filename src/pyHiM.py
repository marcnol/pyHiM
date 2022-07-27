#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol

This file contains routines to process Hi-M datasets

"""
# =============================================================================
# IMPORTS
# =============================================================================q

from datetime import datetime

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

from fileProcessing.fileManagement import Parameters, printLog
from fileProcessing.functionCaller import HiMfunctionCaller, HiM_parseArguments

# =============================================================================
# MAIN
# =============================================================================

def main():
    begin_time = datetime.now()

    runParameters=HiM_parseArguments()

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()
    session1, log1=HiM.session1, HiM.log1

    HiM.lauchDaskScheduler(threadsRequested = runParameters["threads"],maximumLoad=0.8)
    param = Parameters(rootFolder = runParameters["rootFolder"], fileName = 'infoList.json')
    labels=param.param['labels']

    printLog('$ Started logging to: {}'.format(HiM.logFile))
    printLog("$ labels to process: {}\n".format(labels))

    for label in labels:#range(len(HiM.labels2Process)):

        # sets parameters
        param = Parameters(rootFolder = runParameters["rootFolder"], label = label, fileName = 'infoList.json')

        printLog("--------------------------------------------------------------------------")
        printLog(">                  Analyzing label: {}           ".format(param.param["acquisition"]["label"]))
        printLog("--------------------------------------------------------------------------")

        param.param['parallel']=HiM.parallel
        param.param['fileNameMD']=HiM.fileNameMD

        # [projects 3D images in 2d]
        if "makeProjections" in runParameters["cmd"]:
            HiM.makeProjections(param)

        # [registers fiducials using a barcode as reference]
        if "alignImages" in runParameters["cmd"]:
            HiM.alignImages(param, label)

        # [applies registration to DAPI and barcodes]
        if "appliesRegistrations" in runParameters["cmd"]:
            HiM.appliesRegistrations(param, label)

        # [aligns fiducials in 3D]
        if "alignImages3D" in runParameters["cmd"]:
            HiM.alignImages3D(param, label)

        # [segments DAPI and sources in 2D]
        if "segmentMasks" in runParameters["cmd"]:
            HiM.segmentMasks(param, label)

        # [segments masks in 3D]
        if "segmentMasks3D" in runParameters["cmd"]:
            HiM.segmentMasks3D(param, label)

        # [segments sources in 3D]
        if "segmentSources3D" in runParameters["cmd"]:
            HiM.segmentSources3D(param, label)

        # [2D projects all barcodes in an ROI]
        if "projectBarcodes" in runParameters["cmd"]:
            HiM.projectsBarcodes(param, label)

        # [refits spots in 3D]
        if "refitBarcodes3D" in runParameters["cmd"]:
            HiM.refitBarcodes(param, label)

        # [local drift correction]
        if "localDriftCorrection" in runParameters["cmd"]:
            HiM.localDriftCorrection(param, label)

        # [filters barcode localization table]
        if "filter_localizations" in runParameters["cmd"]:
            HiM.filter_localizations(param, label)

        # [registers barcode localization table]
        if "register_localizations" in runParameters["cmd"]:
            HiM.register_localizations(param, label)

        # [build traces]
        if "build_traces" in runParameters["cmd"]:
            HiM.build_traces(param, label)

        # [builds matrices]
        if "build_matrix" in runParameters["cmd"]:
            HiM.build_matrix(param, label)

        # [builds PWD matrix for all folders with images]
        if "buildHiMmatrix" in runParameters["cmd"]:
            HiM.processesPWDmatrices(param, label)

        print("\n")
        del param

    # exits
    HiM.session1.save(HiM.log1)
    printLog("\n===================={}====================\n".format("Normal termination"))

    if runParameters["parallel"]:
        HiM.cluster.close()
        HiM.client.close()

    del HiM

    printLog("Elapsed time: {}".format(datetime.now() - begin_time))


if __name__ == "__main__":
    main()