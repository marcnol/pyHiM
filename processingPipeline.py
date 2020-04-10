#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol
"""

from makeProjections import makeProjections
from alignImages import alignImages, appliesRegistrations
from fileManagement import Parameters, log, session
import os
from datetime import datetime

if __name__ == '__main__':
    begin_time = datetime.now()
    
    #rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
    #rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001'
    rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001_test'
    
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

        print("\n")        
        del param    
    
    # exits
    session1.save(log1)
    log1.addSimpleText("\n===================={}====================\n".format('Normal termination'))
    
    del log1, session1
    print("Elapsed time: {}".format(datetime.now()-begin_time))