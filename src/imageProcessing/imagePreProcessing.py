#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 13:44:47 2021

@author: marcnol

- loads a list of images from a folder
- re-scales exposures
- removes inhomogeneous background
- adjusts image levels
- saves output in TIF

Steps:
    - defines runParameters
    - gets list of images in folder using wildcard, *tif by default
    - loads iteratively images and applies operations described above
    - saves output

"""


import os, argparse, sys, glob
from datetime import datetime

from skimage import io
import numpy as np
from tifffile import imsave
from tqdm import tqdm, trange
from skimage import exposure
from imageProcessing import (
    _removesInhomogeneousBackground,
    imageAdjust,
    savesImageAsBlocks,
    )

from dask.distributed import Client, LocalCluster, get_client, as_completed, fire_and_forget
import multiprocessing

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-O","--outputFile", help="Provide input file ")
    parser.add_argument("-M", "--mode", help="Mode: tif, fits, npy")
    parser.add_argument("-W", "--wildcard", help="For instance '*tif'. If none is give, it will get all the tif files in the folder.")
    parser.add_argument("-Z", "--zPlane", help="z-plane for 3D images")
    parser.add_argument("--lower_threshold", help="lower_threshold for adjusting levels")
    parser.add_argument("--higher_threshold", help="higher_threshold for adjusting levels")
    parser.add_argument("--blockSizeXY", help="blockSizeXY for breaking image into blocks. Default = 256")
    parser.add_argument("--parallel", help="Runs in parallel mode", action="store_true")



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
        runParameters["zPlane"] = 0

    if args.lower_threshold:
        runParameters["lower_threshold"] = float(args.lower_threshold)
    else:
        runParameters["lower_threshold"] = 0.9

    if args.higher_threshold:
        runParameters["higher_threshold"] = float(args.higher_threshold)
    else:
        runParameters["higher_threshold"] = 0.9999

    if args.blockSizeXY:
        runParameters["blockSizeXY"] = float(args.blockSizeXY)
    else:
        runParameters["blockSizeXY"] = 256

    if args.parallel:
        runParameters["parallel"] = args.parallel
    else:
        runParameters["parallel"] = False

    print("\n> Arguments parsed from command line list: {}\nNo problem found.".format(runParameters))

    return runParameters

def lauchDaskScheduler(requestedNumberNodes, maximumLoad=0.6, memoryPerWorker=2000):
    numberCoresAvailable = multiprocessing.cpu_count()

    # we want at least 2 GB per worker
    _, _, free_m = map(int, os.popen("free -t -m").readlines()[-1].split()[1:])

    maxNumberThreads = int(np.min([numberCoresAvailable * maximumLoad, free_m / memoryPerWorker]))

    nThreads = int(np.min([maxNumberThreads, requestedNumberNodes]))

    print("Go to http://localhost:8787/status for information on progress...")

    cluster = LocalCluster(
        n_workers=nThreads,
        # processes=True,
        # threads_per_worker=1,
        # memory_limit='2GB',
        # ip='tcp://localhost:8787',
    )
    client = Client(cluster)

    return client, cluster

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
    blockSizeXY = runParameters["blockSizeXY"]

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

    if runParameters["parallel"]:
        client, cluster = lauchDaskScheduler(6)

    # - loads iteratively images and sums
    i = 0
    for file in tqdm(files2Process):
        # print("Loading and processing: {}".format(file))

        # loads file
        images = list()
        rawImage = io.imread(file).squeeze()
        images.append(rawImage)

        # breaks and saves original image as 3D blocks
        savesImageAsBlocks(rawImage,file,blockSizeXY=blockSizeXY)

        # autoscales exposures
        images = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in images]

        # removes inhomogeneous background
        print("\nRemoving inhomogeneous background...")
        images = [_removesInhomogeneousBackground(x) for x in images]

        # rescales grey levels
        print("\nRescaling grey levels...")
        images = [imageAdjust(x, lower_threshold=runParameters["lower_threshold"],
                              higher_threshold=runParameters["higher_threshold"])[0] for x in images]

        # saves output
        fileName = file.split('.')[0]
        if runParameters["outputFile"] is None:
            outputFile = fileName + '_preProcessed'
        else:
            outputFile = runParameters["outputFile"] + '_' + str(i) + '_preProcessed'
            i+=1

        for index, image in enumerate(images):
            outfile=outputFile + '_index:' + str(index) + '.' + extension
            imsave(outfile, image)
            print("\n> Saving image : {}".format(outfile))

    if runParameters["parallel"]:
        cluster.close()
        client.close()

    print("Elapsed time: {}".format(datetime.now() - begin_time))
