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
import argparse
from datetime import datetime
from matplotlib import pyplot as plt
import logging
import numpy as np
from roipoly import MultiRoi
from astropy.table import Table, Column, vstack
from imageProcessing import Image
from alignBarcodesMasks import processesPWDmatrices
from projectsBarcodes import projectsBarcodes
from fileManagement import folders, writeString2File, saveJSON, loadJSON, RT2fileName
from makeProjections import makeProjections
from alignImages import alignImages, appliesRegistrations
from segmentMasks import segmentMasks
from fileManagement import Parameters, log, writeString2File, session

"""
logging.basicConfig(format='%(levelname)s ''%(processName)-10s : %(asctime)s '
                            '%(module)s.%(funcName)s:%(lineno)s %(message)s',
                    level=logging.INFO)
"""
# =============================================================================
# FUNCTIONS
# =============================================================================q


def createsUserMask(fileName, outputFileName):

    # loads image

    # displays image
    Im = Image()
    Im.data_2D = np.load(fileName).squeeze()
    Im.imageShow(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple ROIs
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    numberROIs = len(multiroi_named.rois)
    print("Number of ROIs drawn: {}".format(numberROIs))

    masks = np.zeros((Im.data_2D.shape[0], Im.data_2D.shape[1], numberROIs))

    # Display image with all ROIs
    Im.imageShow(show=True, normalization="simple")

    roi_names = []
    i = 0
    for name, roi in multiroi_named.rois.items():
        roi.display_roi()
        roi.display_mean(Im.data_2D)
        roi_names.append(name)
        masks[:, :, i] = roi.get_mask(Im.data_2D)
        i += 1
    plt.legend(roi_names, bbox_to_anchor=(1.2, 1.05))
    plt.show()

    # saves result
    np.save(outputFileName, masks)


def processesUserMasks(param, log1, session1, processingList):

    if param.param["segmentedObjects"]["operation"] == "overwrite":

        # session
        sessionName = "processesUserMasks"

        # processes folders and files
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.setsFolders()
        log1.addSimpleText(
            "\n===================={}====================\n".format(sessionName)
        )
        log1.report("folders read: {}".format(len(dataFolder.listFolders)))

        positionROIinformation = param.param["acquisition"]["positionROIinformation"]
        numberMaskedFiles = 0

        allresultsTable = Table()
        for currentFolder in dataFolder.listFolders:
            # currentFolder=dataFolder.listFolders[0] # only one folder processed so far...
            filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
            dataFolder.createsFolders(currentFolder, param)
            log1.report("-------> Processing Folder: {}".format(currentFolder))

            # generates lists of files to process
            param.files2Process(filesFolder)
            log1.report(
                "About to process {} files\n".format(len(param.fileList2Process))
            )

            if len(param.fileList2Process) > 0:
                fileList2Process = [
                    file
                    for file in glob.glob(
                        dataFolder.outputFolders["segmentedObjects"] + os.sep + "*"
                    )
                    if "SNDmask" in file.split("_")
                ]

                if processingList["cleanAllMasks"]:

                    # [clears all SND masks]
                    numberMaskedFiles = len(fileList2Process)
                    for fileName in fileList2Process:
                        os.remove(fileName)
                        log1.report(
                            "Removing SND mask: {}".format(os.path.basename(fileName)),
                            "info",
                        )

                elif (
                    processingList["addMask"] == ""
                    and not processingList["cleanAllMasks"]
                ):

                    # [assigns cells to exsting masks]
                    resultsTable = assignsSNDmask2Cells(
                        fileList2Process, positionROIinformation
                    )
                    allresultsTable = vstack([allresultsTable, resultsTable])
                    # allresultsTable=resultsTable

                    log1.report("assigning masks", "info")

                elif (
                    not processingList["cleanAllMasks"]
                    and len(processingList["addMask"]) > 0
                ):

                    # [makes new set of masks]
                    for fileName in param.fileList2Process:

                        # gets filename information
                        ROI = os.path.basename(fileName).split("_")[
                            positionROIinformation
                        ]

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
                            createsUserMask(registeredFileName, outputFileName)
                            numberMaskedFiles += 1
                            log1.report(
                                "Segmented image for SND channel ROI#{}: {}".format(
                                    ROI, outputFileName
                                ),
                                "info",
                            )
                        else:
                            log1.report(
                                "Could not find registered image for SND channel:{}".format(
                                    registeredFileName
                                ),
                                "error",
                            )

        tableOutputFileName = (
            dataFolder.outputFolders["segmentedObjects"]
            + os.sep
            + "SNDassignedCells.ecsv"
        )
        if os.path.exists(tableOutputFileName):
            previousresultsTable = Table.read(
                tableOutputFileName, format="ascii.ecsv"
            )  # ascii.ecsv
            allresultsTable = vstack([previousresultsTable, allresultsTable])
        if len(allresultsTable)>0:
            allresultsTable.write(tableOutputFileName, format="ascii.ecsv", overwrite=True)

        return numberMaskedFiles
    else:
        return 0


def assignsSNDmask2Cells(fileList2Process, positionROIinformation):
    resultsTable = Table()

    numberFilesProcessed = 0
    for fileName in fileList2Process:
        ROI = os.path.basename(fileName).split("_")[positionROIinformation]

        # [checks if DAPI mask exists for the file to process]
        fileNameDAPImask = (
            os.path.dirname(fileName)
            + os.sep
            + "_".join(os.path.basename(fileName).split("_")[0:7])
            + "_ch00_Masks.npy"
        )
        if os.path.exists(fileNameDAPImask):

            # load DAPI mask
            maskDAPI = Image()
            maskDAPI.data_2D = np.load(fileNameDAPImask).squeeze()

            # load SND mask
            maskSND = Image()
            maskSND.data_2D = np.load(fileName).squeeze()

            # matches cells and SND masks

            # matches cells and SND masks
            newMatrix = maskDAPI.data_2D * maskSND.data_2D
            cellsWithinMask = np.unique(newMatrix)

            if len(cellsWithinMask) > 0:

                newTable = Table()
                colMask = Column(
                    [fileName.split("_")[-1].split(".")[0]] * len(cellsWithinMask),
                    name="MaskID #",
                    dtype=str,
                )
                colROI = Column(
                    int(ROI) * np.ones(len(cellsWithinMask)), name="ROI #", dtype=int
                )
                colCellID = Column(cellsWithinMask, name="CellID #", dtype=int)
                newTable.add_column(colROI, index=0)
                newTable.add_column(colCellID, index=1)
                newTable.add_column(colMask, index=2)

                # newTable=Table()
                # newTable['CellID #']=cellsWithinMask
                # newTable['ROI #']=int(ROI)*np.ones(len(cellsWithinMask))
                # newTable['MaskID #']=[fileName.split('_')[-1].split('.')[0]]*len(cellsWithinMask)

                resultsTable = vstack([resultsTable, newTable])
                del newTable

            # outputs list of cells within (1) and outside SND mask (0)

            log1.report(
                "Matched DAPI and SND masks: \n SND: {}\n DAPI: {}\n".format(
                    fileName, fileNameDAPImask
                ),
                "info",
            )
            numberFilesProcessed += 1

        else:

            log1.report(
                "Could not find expected DAPI mask: {}".format(fileNameDAPImask),
                "warning",
            )

    return resultsTable

    # barcodesCoordinates.write(outputFile,format='ascii.ecsv',overwrite=True)


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

    print(
        "\n--------------------------------------------------------------------------"
    )
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

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName = "processSNDchannel"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder, fileNameRoot=sessionName)
    log1.addSimpleText(
        "\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(
            sessionName
        )
    )
    log1.report("Process SND channel MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD,
        "# Process SND channel {}".format(now.strftime("%d/%m/%Y %H:%M:%S")),
        "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # processes Secondary masks
        if label == "RNA":
            numberMaskedFiles = processesUserMasks(
                param, log1, session1, processingList
            )
            if numberMaskedFiles > 0:
                log1.report("Files processed: {}".format(numberMaskedFiles), "info")
            else:
                log1.report("Nothing processed...", "warning")
        print("\n")
        del param

    # exits
    session1.save(log1)
    log1.addSimpleText(
        "\n===================={}====================\n".format("Normal termination")
    )

    del log1, session1
    print("Elapsed time: {}".format(datetime.now() - begin_time))
