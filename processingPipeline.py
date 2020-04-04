#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:11:01 2020

@author: marcnol
"""

from makeProjections import makeProjections
from fileManagement import Parameters, folders, log
import os
#from os import path

if __name__ == '__main__':
    
    rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
    
    parameterFiles=['infoList_DAPI.json',
                    'infoList_barcode.json',
                    'infoList_fiducial.json']

    # setup logs
    logFile='makeProjections.log'
    logFileName=rootFolder+os.sep+logFile
    log1=log(logFileName)
    log1.eraseFile()
    log1.report("Starting to log to: {}".format(logFileName))
 
    for label in parameterFiles:
        log1.report("*******Analyzing label: {}**********".format(label))
        log1.report("*******************************************************\n\n".format(label))
       
        # sets parameters
        param = Parameters()
        param.initializeStandardParameters()
        paramFile = rootFolder+os.sep+label
        param.loadParametersFile(paramFile)
        
        param.param['rootFolder']=rootFolder
        
        # [projects 3D images in 2d]
        makeProjections(param,log1)

        print("\n\n")        
        del param
    
    log1.report('** Normal exit **')
