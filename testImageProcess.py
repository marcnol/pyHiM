#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:10 2020

@author: marcnol

"""


#Imports
import glob,os,sys
#from copy import deepcopy
import matplotlib.pylab as plt
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib import cm
from skimage import io
from imageProcessing import Image
import parameters as parameters
import logging

#import cPickle as pickle
#from scipy.spatial import ConvexHull
#from scipy.spatial.distance import cdist
#from tqdm import tqdm_notebook as tqdm

#Internal packages
#import Fitting_v4 as ft
#import workers_cells_v3 as wkc

def loggerInitialize(logFileName='tmp.log'):
    # create logger
    logger = logging.getLogger('Main logger')
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(filename='example.log', filemode='w')
    
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # add formatter to ch
    ch.setFormatter(formatter)
    
    # add ch to logger
    logger.addHandler(ch)
    
    return logger

if __name__ == '__main__':
    
    param = parameters.Parameters()

    # defines folder
    folder='/home/marcnol/Documents/Images/'
    tiffFile='scan_001_DAPI_001_ROI_converted_decon_ch00.tif'
    logFile='testImageProcess.log'
    fileName=folder+tiffFile
    logFileName=folder+logFile
    log=loggerInitialize()
    #log.basicConfig(filename=logFileName,level=log.INFO)
    
    # creates image object
    Im = Image()
    
    # loads image
    Im.loadImage(fileName)

    if Im.fileName:
        Im.printImageProperties()
        
    Im.zProjectionRange(param,log)
    print("zRange={}".format(Im.zRange))
    
    # calculates MIP
    #Im.maxIntensityProjection()
    
    # Show an image in each subplot
    '''
    fig, ax = plt.subplots(nrows=1, ncols=3)
    ax[0].imshow(Im.data[25])
    ax[1].imshow(Im.normalizeImage(Im.data)[25])
    ax[2].imshow(Im.data_2D,cmap=cm.plasma)
    '''
    
    Im.imageShow(save=False)
    log.info('Normal exit.')
    
    # exits

    
    
