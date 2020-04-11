#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 23:17:58 2020

@author: marcnol
"""


import glob,os,sys
from os import path
import matplotlib.pylab as plt
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib import cm
from skimage import io
from imageProcessing import Image
from fileManagement import folders
from fileManagement import session,writeString2File

def makes2DProjectionsFile(fileName,param,log1,session1,dataFolder):
    
    if fileName in session1.data and param.param['zProject']['operation']!='overwrite':
        if param.param['zProject']['display']:
            # creates image object
            Im = Image()
            Im.loadImage2D(fileName,log1,dataFolder)
            Im.imageShow()
        log1.report("File already projected: {}".format(os.path.basename(fileName)))        
    else:
        
        log1.report('Analysing file: {}'.format(os.path.basename(fileName)))
          
        # creates image object
        Im = Image()
        
        # loads image
        Im.loadImage(fileName)
        
        # makes actual 2d projection
        Im.zProjectionRange(param,log1)
        
        # outputs information from file
        if Im.fileName:
            Im.printImageProperties()
        
        # saves output 2d zProjection as png
        if param.param['zProject']['display']:
            pngFileName=dataFolder.zProjectFolder+os.sep+os.path.basename(fileName)+'_2d.png'
            Im.imageShow(save=param.param['zProject']['saveImage'],outputName=pngFileName)
            writeString2File(log1.fileNameMD,"{}\n ![]({})\n".format(os.path.basename(fileName),pngFileName),'a') # initialises MD file

        # saves output 2d zProjection as matrix
        Im.saveImage2D(log1,dataFolder.zProjectFolder)
        
        del Im
        
def makeProjections(param,log1,session1):
    sessionName='makesProjections'
 
    # processes folders and files 
    dataFolder=folders(param.param['rootFolder'])
    dataFolder.setsFolders()
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report('folders read: {}'.format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD,"## {}: {}\n".format(sessionName,param.param['acquisition']['label']),'a') # initialises MD file
    
    for currentFolder in dataFolder.listFolders:
        #currentFolder=dataFolder.listFolders[0]
        filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder,param)
    
        # generates lists of files to process    
        param.files2Process(filesFolder)
        log1.report('About to read {} files\n'.format(len(param.fileList2Process)))
        
        for fileName in param.fileList2Process:
            makes2DProjectionsFile(fileName,param,log1,session1,dataFolder)
            session1.add(fileName,sessionName)


    
    
