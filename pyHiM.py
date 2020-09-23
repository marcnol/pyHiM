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

from datetime import datetime

from fileProcessing.fileManagement import Parameters
from fileProcessing.functionCaller import HiMfunctionCaller, HiM_parseArguments

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")
               
# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    runParameters=HiM_parseArguments()    

    HiM = HiMfunctionCaller(runParameters, sessionName="HiM_analysis")
    HiM.initialize()
    session1, log1=HiM.session1, HiM.log1
    
    HiM.lauchDaskScheduler()

    for ilabel in range(len(HiM.labels2Process)):
        HiM.log1.addSimpleText("**Analyzing label: {}**".format(HiM.labels2Process[ilabel]["label"]))

        # sets parameters
        param = Parameters(runParameters["rootFolder"], HiM.labels2Process[ilabel]["parameterFile"])
        param.param['parallel']=HiM.parallel

        # [projects 3D images in 2d]
        HiM.makeProjections(param)

        # [registers fiducials using a barcode as reference]
        HiM.alignImages(param, ilabel)
        
        # [applies registration to DAPI and barcodes]
        HiM.appliesRegistrations(param, ilabel)
        
        # [segments DAPI and spot masks]
        HiM.segmentMasks(param, ilabel)

        # [2D projects all barcodes in an ROI]
        HiM.projectsBarcodes(param, ilabel)
            
        # [refits spots in 3D]
        HiM.refitBarcodes(param, ilabel)
           
        # [local drift correction]
        HiM.localDriftCorrection(param, ilabel)

        # [builds PWD matrix for all folders with images]
        HiM.processesPWDmatrices(param, ilabel)

        print("\n")
        del param

    # exits
    HiM.session1.save(HiM.log1)
    HiM.log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

    if runParameters["parallel"]:
        HiM.cluster.close()   
        HiM.client.close()

    del HiM
    
    print("Elapsed time: {}".format(datetime.now() - begin_time))
