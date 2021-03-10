#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:54:49 2021

@author: marcnol


- loads a list of images from a folder
- segments image in 3D
- saves output in TIF

Steps:
    - defines runParameters
1    - loads iteratively images and segments volumes 
    - saves output

"""


import os, argparse, sys, glob
from datetime import datetime
import numpy as np
from skimage import io
import numpy as np
from tifffile import imsave
from tqdm import tqdm, trange
from skimage import exposure
from imageProcessing import (
    _removesInhomogeneousBackground2D,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    )

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-O","--outputFile", help="Provide input file ")

    parser.add_argument("-M", "--mode", help="Mode: tif, fits, npy")
    parser.add_argument("-W", "--wildcard", help="For instance '*tif'. If none is give, it will get all the tif files in the folder.")
    parser.add_argument("-Z", "--zPlane", help="z-plane for 3D images")
    parser.add_argument("--threshold_over_std", help="threshold_over_std for thresholding. Default = 2 (ie. 2% of max intensity range)")
    parser.add_argument("--area_min", help="area_min of each object, in pixels. Default = 3")
    parser.add_argument("--area_max", help="area_max of each object, in pixels. Default = 1000")
    parser.add_argument("--nlevels", help="nlevels for deblending. Default = 32")
    parser.add_argument("--contrast", help="contrast for deblending. Default = 0.001")
    parser.add_argument("--blockSizeXY", help="blockSizeXY for breaking image into blocks. Default = 256")
    
    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    runParameters = dict()

    if args.rootFolder:
        runParameters["rootFolder"] = args.rootFolder
    else:
        print("\n> rootFolder NOT FOUND, using PWD")
        runParameters["rootFolder"] = os.getenv("PWD")  # os.getcwd()

    if args.outputFile:
        runParameters["outputFile"] = runParameters["rootFolder"] + os.sep + args.outputFile
    else:
        runParameters["outputFile"] = None

    if args.mode:
        runParameters["mode"] = args.mode
    else:
        runParameters["mode"] = "tif"

    if args.wildcard:
        runParameters["wildcard"] = args.wildcard
    else:
        runParameters["wildcard"] = "None"

    if args.zPlane:
        runParameters["zPlane"] = int(args.zPlane)
    else:
        runParameters["zPlane"] = np.nan

    if args.threshold_over_std:
        runParameters["threshold_over_std"] = int(args.threshold_over_std)
    else:
        runParameters["threshold_over_std"] = 2
        
    if args.area_min:
        runParameters["area_min"] = int(args.area_min)
    else:
        runParameters["area_min"] = 3
        
    if args.area_max:
        runParameters["area_max"] = int(args.area_max)
    else:
        runParameters["area_max"] = 1000

    if args.nlevels:
        runParameters["nlevels"] = int(args.nlevels)
    else:
        runParameters["nlevels"] = 32
        
    if args.contrast:
        runParameters["contrast"] = float(args.contrast)
    else:
        runParameters["contrast"] = 0.001
                
    if args.blockSizeXY:
        runParameters["blockSizeXY"] = float(args.blockSizeXY)
    else:
        runParameters["blockSizeXY"] = 256
        
    print("\n> Arguments parsed from command line list: {}\nNo problem found.".format(runParameters))

    return runParameters


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    # - defines runParameters
    runParameters = parseArguments()

    # - gets list of images in folder using wildcard, *tif by default
    if "None" in runParameters["wildcard"]:
        extension = runParameters["mode"]
        searchString = runParameters["rootFolder"] + os.sep + "*." + extension
    else:
        print("wildcard: {}".format(runParameters["wildcard"]))
        searchString = runParameters["rootFolder"] + os.sep + runParameters["wildcard"]
        extension = runParameters["wildcard"].split('.')[1]

    files2Process = glob.glob(searchString)
    files2Process.sort()
    
    if len(files2Process) < 1:
        print("No images found!")
        sys.exit()
    else:
        print("Number of images to process: {}".format(len(files2Process)))

    # sets expected image size
    img = io.imread(files2Process[0]).squeeze()
    imgSize = img.shape
    print("Expected image sizes: {}".format(imgSize))

    # creates output image
    outputImage = np.zeros(imgSize)

    # - loads iteratively images and sums
    i = 0
    threshold_over_std = runParameters["threshold_over_std"]
    area_min = runParameters["area_min"]
    area_max = runParameters["area_max"]
    nlevels = runParameters["nlevels"]
    contrast = runParameters["contrast"]
    blockSizeXY = runParameters["blockSizeXY"]
    
    for file in tqdm(files2Process):
        print("\nLoading and processing: {}".format(file))

        # loads file
        image3D = io.imread(file).squeeze()

        print("Segmenting images...")
        binaryMask, newLabeledImage = _segments3DvolumesByThresholding(image3D,
                                                           threshold_over_std=threshold_over_std, 
                                                           sigma = 3, 
                                                           boxSize=(32, 32),
                                                           filter_size=(3, 3),
                                                           area_min=area_min,
                                                           area_max=area_max,
                                                           nlevels=nlevels,
                                                           contrast=contrast,
                                                           deblend3D=True)
        # breaks and saves image in 3D blocks
        savesImageAsBlocks(newLabeledImage,file,blockSizeXY=blockSizeXY)

        # - save outputs
        fileName = file.split('.')[0]
        if runParameters["outputFile"] is None:
            outputFile = fileName + '_segmented'
        else:
            outputFile = runParameters["outputFile"] + '_' + str(i) + '_segmented'
            i+=1            
        
        outfile=outputFile + '.' + extension
        print("\n> Saving image : {}".format(outfile))
        imsave(outfile, newLabeledImage)

            
