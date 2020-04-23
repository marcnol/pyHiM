#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol
"""
# =============================================================================
# IMPORTS
# =============================================================================

from makeProjections import makeProjections
from alignImages import alignImages, appliesRegistrations
from segmentMasks import segmentMasks
from fileManagement import Parameters, log,writeString2File, session
import os
import argparse
from datetime import datetime
from alignBarcodesMasks import processesPWDmatrices

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
        #rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
        #rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
        
        rootFolder='/home/marcnol/data/Experiment_20/Embryo_1'
    
        #rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001'
        #rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001_test'
        #rootFolder='/mnt/PALM_dataserv/DATA/merFISH_2019/Experiment_15/2019_05_15/deconvolved_RT_1/006_Embryo/rawData'

    print("parameters> rootFolder: {}".format(rootFolder))


    
    labels2Process = [{'label':'fiducial', 'parameterFile': 'infoList_fiducial.json'},
                      {'label':'DAPI', 'parameterFile': 'infoList_DAPI.json'},
                      {'label': 'barcode', 'parameterFile': 'infoList_barcode.json'}] 
                        
     # session
    now=datetime.now()
    sessionRootName=now.strftime("%d%m%Y_%H%M%S")
    sessionFileName=rootFolder+os.sep+'Session_'+sessionRootName+'.json'
    session1=session('processingPipeline',sessionFileName)
 
    # setup logs
    logFileName=rootFolder+os.sep+'HiM_analysis'+sessionRootName+'.log'
    log1=log(logFileName)
    log1.eraseFile()
    log1.report("Starting to log to: {}".format(logFileName))
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format('processingPipeline'))
    writeString2File(log1.fileNameMD,"# Hi-M analysis {}".format(now.strftime("%d/%m/%Y %H:%M:%S")),'w') # initialises MD file
    
    for ilabel in range(len(labels2Process)):
        label=labels2Process[ilabel]['label']
        labelParameterFile=labels2Process[ilabel]['parameterFile']

        log1.addSimpleText("**Analyzing label: {}**".format(label))
       
        # sets parameters
        param = Parameters(labelParameterFile)
        param.initializeStandardParameters()
        paramFile = rootFolder+os.sep+labelParameterFile
        param.loadParametersFile(paramFile)
        
        param.param['rootFolder']=rootFolder
        
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
            segmentMasks(param,log1,session1)
            
        # [refits spots in 3D]

        # [local drift correction]
        
        
        print("\n")        
        del param    

    # [builds PWD matrix for all folders with images]
    ilabel=1 # uses DAPI for parameters file
    label=labels2Process[ilabel]['label'] 
    labelParameterFile=labels2Process[1]['parameterFile']
   
    # sets parameters
    param = Parameters(labelParameterFile)
    param.loadParametersFile(rootFolder+os.sep+labelParameterFile)
    param.param['rootFolder']=rootFolder
  
    processesPWDmatrices(param,log1,session1)        
    
    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format('Normal termination'))
    
    del log1, session1
    print("Elapsed time: {}".format(datetime.now()-begin_time))