#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 23:17:58 2020

@author: marcnol

This file contains functions to project 3D images to 2D

Operation will be defined in the parameters file. Options are:
    - user-defined range
    - all z range
    - optimal range based on detection of focal plane and use of user defined window around it
    

"""
# =============================================================================
# IMPORTS
# =============================================================================


import glob, os

# import matplotlib.pylab as plt
import numpy as np
import argparse
from datetime import datetime
from dask.distributed import Client, wait, LocalCluster
from multiprocessing.pool import ThreadPool
import threading

# import cv2
import matplotlib.pyplot as plt
from imageProcessing.imageProcessing import Image

from fileProcessing.fileManagement import (
    folders,session, writeString2File, folders, session, log, Parameters)


# =============================================================================
# FUNCTIONS
# =============================================================================


def makes2DProjectionsFile(fileName, param, log1, session1, dataFolder):

    if fileName in session1.data and param.param["zProject"]["operation"] != "overwrite":
        # creates image object
        Im = Image()
        Im.loadImage2D(fileName, log1, dataFolder.outputFolders["zProject"])
        if param.param["zProject"]["display"]:
            Im.imageShow()
        log1.report("File already projected: {}".format(os.path.basename(fileName)))
    else:

        log1.report("Analysing file: {}".format(os.path.basename(fileName)))

        # creates image object
        Im = Image()

        # loads image
        Im.loadImage(fileName)

        # makes actual 2d projection
        Im.zProjectionRange(param, log1)

        # outputs information from file
        if Im.fileName:
            Im.printImageProperties()

        # saves output 2d zProjection as png
        if param.param["zProject"]["display"]:
            pngFileName = dataFolder.outputFolders["zProject"] + os.sep + os.path.basename(fileName) + "_2d.png"
            Im.imageShow(save=param.param["zProject"]["saveImage"], outputName=pngFileName)
            writeString2File(
                log1.fileNameMD, "{}\n ![]({})\n".format(os.path.basename(fileName), pngFileName), "a",
            )  # initialises MD file

        # saves output 2d zProjection as matrix
        Im.saveImage2D(log1, dataFolder.outputFolders["zProject"])

        del Im


def makeProjections(param, log1, session1,fileName=None):
    sessionName = "makesProjections"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
    )  # initialises MD file

        
    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        if param.param['parallel']:
            threads=list()
            files2ProcessFiltered = [x for x in param.fileList2Process if \
                                     (fileName==None) \
                                     or (fileName!=None \
                                     and (os.path.basename(x) in [os.path.basename(x1) for x1 in fileName]))]

            if len(files2ProcessFiltered)>0:
                print("Cluster with {} workers started".format(len(files2ProcessFiltered)))

                # dask
                with LocalCluster(n_workers=len(files2ProcessFiltered)
                                ) as cluster, Client(cluster) as client:

                    # client = Client(processes=False)#,n_workers=len(files2ProcessFiltered))
                    threads=[client.submit(makes2DProjectionsFile,x, param, log1, session1, dataFolder) for x in files2ProcessFiltered]            
                    
                    for index, thread in enumerate(threads):
                        print("Waiting for thread: {}".format(index+1))
                        wait(threads)        
                
                # ThreadPool
                # pool = ThreadPool(processes=len(param.fileList2Process))
                # threads=[pool.apply_async(makes2DProjectionsFile,args=(x, param, log1, session1, dataFolder)) for x in files2ProcessFiltered]            
                # pool.close()
                # pool.join()

                # simple threads
                # for x in files2ProcessFiltered:
                #     x1=threading.Thread(target=makes2DProjectionsFile,args=(x, param, log1, session1, dataFolder))
                #     threads.append(x1)
                #     x1.start() 
                # for index, thread in enumerate(threads):
                #     print("Waiting for thread: {}".format(index+1))
                #     thread.join()
        else:
            
            for index, fileName2Process in enumerate(param.fileList2Process):
    
                if (fileName==None) or (fileName!=None and (os.path.basename(fileName2Process) in [os.path.basename(x) for x in fileName])):
                    makes2DProjectionsFile(fileName2Process, param, log1, session1, dataFolder)                    
                    session1.add(fileName2Process, sessionName)
                else:
                    pass

