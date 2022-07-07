#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:21:01 2020

@author: marcnol

processSNDchannel


"""


# =============================================================================
# IMPORTS
# =============================================================================q

import os, glob
try:
    del os.environ["MPLBACKEND"] # unset before importing matplotlib
except:
    print("No environment variable MPLBACKEND found. Continuing anyway.")
import argparse
from datetime import datetime
from matplotlib import pyplot as plt
import numpy as np
from roipoly import MultiRoi
from astropy.table import Table, Column, vstack
from astropy.visualization import simple_norm

from imageProcessing.imageProcessing import Image

from fileProcessing.fileManagement import ( Parameters, log, session,
                                           folders, writeString2File)


# =============================================================================
# FUNCTIONS
# =============================================================================q

def imageShow(data_2D, normalization = "simple",size=(10, 10)):
    fig = plt.figure()
    fig.set_size_inches(size)

    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()

    if normalization == "simple":
        norm = simple_norm(data_2D, "sqrt", percent=99.9)

    fig.add_axes(ax)
    ax.imshow(data_2D, origin="lower", cmap="Greys_r", norm=norm)
    ax.set_title("2D Data")
    return fig,ax

def createsUserMask(param,log1,fileName, outputFileName):

    # loads image

    # displays image
    Im = Image(param,log1)
    Im.data_2D = np.load(fileName).squeeze()
    Im.imageShow(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple ROIs
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    numberROIs = len(multiroi_named.rois)
    print("Number of ROIs drawn: {}".format(numberROIs))

    masks = np.zeros((Im.data_2D.shape[0], Im.data_2D.shape[1], numberROIs))

    # Display image with all ROIs
    fig,ax = imageShow(Im.data_2D)

    roi_names = []
    i = 0
    for name, roi in multiroi_named.rois.items():
        roi.display_roi()
        roi.display_mean(Im.data_2D)
        roi_names.append(name)
        masks[:, :, i] = roi.get_mask(Im.data_2D)
        i += 1
    plt.legend(roi_names, bbox_to_anchor=(1.2, 1.05))
    #plt.show()

    print("Saving and closing image with ROIs...")
    fig.savefig(outputFileName + "_segmentedSources.png")
    plt.close(fig)

    # saves result
    np.save(outputFileName, masks)


def processesUserMasks(param, log1, processingList):

    if param.param["segmentedObjects"]["operation"] == "overwrite":

        # session
        sessionName = "processesUserMasks"

        # processes folders and files
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.setsFolders()
        log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        log1.report("folders read: {}".format(len(dataFolder.listFolders)))

        positionROIinformation = param.param["acquisition"]["positionROIinformation"]
        numberMaskedFiles = 0

        allresultsTable = Table()
        for currentFolder in dataFolder.listFolders:

            filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
            dataFolder.createsFolders(currentFolder, param)
            log1.report("-------> Processing Folder: {}".format(currentFolder))

            # generates lists of files to process
            param.files2Process(filesFolder)
            log1.report("About to process {} files\n".format(len(param.fileList2Process)))

            if len(param.fileList2Process) > 0:
                fileList2Process = [
                    file
                    for file in glob.glob(dataFolder.outputFolders["segmentedObjects"] + os.sep + "*")
                    if "SNDmask" in file.split("_")
                ]

                if processingList["cleanAllMasks"]:

                    # [clears all SND masks]
                    numberMaskedFiles = len(fileList2Process)
                    for fileName in fileList2Process:
                        os.remove(fileName)
                        log1.report(
                            "Removing SND mask: {}".format(os.path.basename(fileName)), "info",
                        )

                elif processingList["addMask"] == "" and not processingList["cleanAllMasks"]:

                    # [assigns cells to exsting masks]
                    resultsTable = assignsSNDmask2Cells(fileList2Process, positionROIinformation)
                    allresultsTable = vstack([allresultsTable, resultsTable])

                    log1.report("assigning masks", "info")

                elif not processingList["cleanAllMasks"] and len(processingList["addMask"]) > 0:

                    # [makes new set of masks]
                    for fileName in param.fileList2Process:

                        # gets filename information
                        ROI = os.path.basename(fileName).split("_")[positionROIinformation]

                        # checks that SND channel was projected and aligned
                        registeredFileName = (
                            dataFolder.outputFolders["alignImages"]
                            + os.sep
                            + os.path.basename(fileName).split(".")[0]
                            + "_2d_registered.npy"
                        )

                        if os.path.exists(registeredFileName):
                            outputFileName = (
                                dataFolder.outputFolders["segmentedObjects"]
                                + os.sep
                                + os.path.basename(registeredFileName).split(".")[0]
                                + "_SNDmask"
                                + "_"
                                + processingList["addMask"]
                                + ".npy"
                            )
                            createsUserMask(param,log1,registeredFileName, outputFileName)
                            numberMaskedFiles += 1
                            log1.report(
                                "Segmented image for SND channel ROI#{}: {}".format(ROI, outputFileName), "info",
                            )
                        else:
                            log1.report(
                                "Could not find registered image for SND channel:{}".format(registeredFileName),
                                "error",
                            )

        tableOutputFileName = dataFolder.outputFolders["segmentedObjects"] + os.sep + "SNDassignedCells.ecsv"
        print(f"\n>Number of cells found: {len(allresultsTable)}")
        print("\n>Saving output to: {}".format(tableOutputFileName))
        if os.path.exists(tableOutputFileName):
            print("\Results appended to {}".format(tableOutputFileName))
            previousresultsTable = Table.read(tableOutputFileName, format="ascii.ecsv")  # ascii.ecsv
            allresultsTable = vstack([previousresultsTable, allresultsTable])
        if len(allresultsTable) > 0:
            print("\nTable written")
            allresultsTable.write(tableOutputFileName, format="ascii.ecsv", overwrite=True)

        return numberMaskedFiles
    else:
        return 0


def assignsSNDmask2Cells(fileList2Process, positionROIinformation):
    resultsTable = Table()

    numberFilesProcessed = 0
    #print(f"\nfiles2Process: {fileList2Process}")

    for fileName in fileList2Process:
        print(f"\n-----> Processing: {fileName}")
        ROI = os.path.basename(fileName).split("_")[positionROIinformation]

        # [checks if DAPI mask exists for the file to process]
        fileNameDAPImask = (
            os.path.dirname(fileName)
            + os.sep
            + "_".join(os.path.basename(fileName).split("_")[0:7])
            + "_ch00_Masks.npy"
        )

        print(f"\nWill search masks in: {fileNameDAPImask}")

        if os.path.exists(fileNameDAPImask) and fileName.split('.')[-1] == 'npy':

            # load DAPI mask
            maskDAPI = Image()
            maskDAPI.data_2D = np.load(fileNameDAPImask).squeeze()

            # load SND mask
            print(f"\nWill attemp to match masks from: {fileName}")
            maskSND = Image()
            maskSND.data_2D = np.load(fileName, allow_pickle=False).squeeze()

            # matches cells and SND masks
            newMatrix = maskDAPI.data_2D * maskSND.data_2D
            cellsWithinMask = np.unique(newMatrix)

            if len(cellsWithinMask) > 0:

                newTable = Table()
                colMask = Column(
                    [fileName.split("_")[-1].split(".")[0]] * len(cellsWithinMask), name="MaskID #", dtype=str,
                )
                colROI = Column(int(ROI) * np.ones(len(cellsWithinMask)), name="ROI #", dtype=int)
                colCellID = Column(cellsWithinMask, name="CellID #", dtype=int)
                newTable.add_column(colROI, index=0)
                newTable.add_column(colCellID, index=1)
                newTable.add_column(colMask, index=2)

                resultsTable = vstack([resultsTable, newTable])
                del newTable

            # outputs list of cells within (1) and outside SND mask (0)
            log1.report(
                "Matched DAPI and SND masks: \n SND: {}\n DAPI: {}\n".format(fileName, fileNameDAPImask), "info",
            )
            numberFilesProcessed += 1

        else:
            print("\nERROR: Could not find expected DAPI mask: {}\n or file does not have ({}) the expected extension (.npy)".format(fileNameDAPImask,fileName.split('.')[-1]))

    return resultsTable

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    begin_time = datetime.now()

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    parser.add_argument("-A", "--addMask", help="Add manual segmentation")
    parser.add_argument("--cleanAllMasks", help="Clear all masks", action="store_true")

    args = parser.parse_args()

    print("\n--------------------------------------------------------------------------")
    processingList = {}

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder = "."
        # rootFolder='/home/marcnol/data/Experiment_4/0_Embryo/alignImages/'
        # fileNameRNA = rootFolder+'scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy'

    if args.addMask:
        processingList["addMask"] = args.addMask
    else:
        processingList["addMask"] = ""

    processingList["cleanAllMasks"] = args.cleanAllMasks

    print("parameters> rootFolder: {}".format(rootFolder))
    now = datetime.now()


    param = Parameters(rootFolder = rootFolder, fileName = 'infoList.json')
    labels=param.param['labels']

    sessionName = "processSNDchannel"

    # setup logs
    log1 = log(rootFolder, fileNameRoot=sessionName)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(sessionName))
    log1.report("Process SND channel MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Process SND channel {}".format(now.strftime("%d/%m/%Y %H:%M:%S")), "w",
    )  # initialises MD file

    for label in labels:

        # sets parameters
        param = Parameters(rootFolder = rootFolder, label = label, fileName = 'infoList.json')
        print("**Analyzing label: {}**".format(label))

        # processes Secondary masks
        if label == "RNA":
            numberMaskedFiles = processesUserMasks(param, log1, processingList)
            if numberMaskedFiles > 0:
                print("Files processed: {}".format(numberMaskedFiles), "info")
            else:
                print("Nothing processed...", "warning")
        print("\n")
        del param

    # exits
    print("\n===================={}====================\n".format("Normal termination"))

    print("Elapsed time: {}".format(datetime.now() - begin_time))
