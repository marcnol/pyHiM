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

from makeProjections import makeProjections
from alignImages import alignImages, appliesRegistrations
from segmentMasks import segmentMasks
from fileManagement import Parameters, log,writeString2File, session
import os
import argparse
from datetime import datetime
from alignBarcodesMasks import processesPWDmatrices
from projectsBarcodes import projectsBarcodes


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    begin_time = datetime.now()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-F","--rootFolder", help="Folder with images")
    args = parser.parse_args()
    
    print("\n--------------------------------------------------------------------------")
    
    if args.rootFolder:
        rootFolder=args.rootFolder
    else:
        #rootFolder='/home/marcnol/data/Embryo_debug_dataset'
        rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        
        #rootFolder='/home/marcnol/data/Experiment_20/Embryo_1'
    
        #rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001'
        #rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001_test'
        #rootFolder='/mnt/PALM_dataserv/DATA/merFISH_2019/Experiment_15/2019_05_15/deconvolved_RT_1/006_Embryo/rawData'

    print("parameters> rootFolder: {}".format(rootFolder))
    now=datetime.now()
    
    labels2Process = [{'label':'fiducial', 'parameterFile': 'infoList_fiducial.json'},
                      {'label': 'barcode', 'parameterFile': 'infoList_barcode.json'},
                      {'label':'DAPI', 'parameterFile': 'infoList_DAPI.json'},
                      {'label':'RNA', 'parameterFile': 'infoList_RNA.json'}]
                        
     # session
    session1=session(rootFolder,'processingPipeline')
 
    # setup logs
    log1=log(rootFolder)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format('processingPipeline'))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(log1.fileNameMD,"# Hi-M analysis {}".format(now.strftime("%d/%m/%Y %H:%M:%S")),'w') # initialises MD file
    
    for ilabel in range(len(labels2Process)):
        label=labels2Process[ilabel]['label']
        labelParameterFile=labels2Process[ilabel]['parameterFile']
        log1.addSimpleText("**Analyzing label: {}**".format(label))
       
        # sets parameters
        param = Parameters(rootFolder,labelParameterFile)
        
        # [projects 3D images in 2d]
        makeProjections(param,log1,session1)

        # [registers fiducials using a barcode as reference]
        if label=='fiducial' and param.param['acquisition']['label']=='fiducial':
            log1.report('Making image registrations, ilabel: {}, label: {}'.format(ilabel,label),'info')
            alignImages(param,log1,session1)

        # [applies registration to DAPI and barcodes]
        if label!='fiducial' and param.param['acquisition']['label']!='fiducial':
            log1.report('Applying image registrations, ilabel: {}, label: {}'.format(ilabel,label),'info')
            appliesRegistrations(param,log1,session1)

            # [segments DAPI and spot masks]
            if label!='RNA' and param.param['acquisition']['label']!='RNA':
                segmentMasks(param,log1,session1)

        # [2D projects all barcodes in an ROI]
        if label=='barcode':
            projectsBarcodes(param,log1,session1)
            
        # [refits spots in 3D]

        # [local drift correction]
        
        # [builds PWD matrix for all folders with images]
        if label=='DAPI':
            processesPWDmatrices(param,log1,session1)        
        
        print("\n")        
        del param    

   # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format('Normal termination'))
    
    del log1, session1
    print("Elapsed time: {}".format(datetime.now()-begin_time))