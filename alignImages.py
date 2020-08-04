#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 09:22:18 2020

@author: marcnol

Sets of functions that do alignment of 2D fiducial images. It also contains
code to apply these alignments to other channels (DAPI/ barcodes)

For the time being alignment is purely based on optimized sub-pixed accuracy
image cross correlation

"""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import os, glob
from skimage.feature.register_translation import _upsampled_dft

# from scipy.ndimage import fourier_shift
from imageProcessing import (
    Image,
    save2imagesRGB,
    saveImage2Dcmd,
    align2ImagesCrossCorrelation,
)
from fileManagement import folders, writeString2File, saveJSON, loadJSON, RT2fileName
from astropy.table import Table
from scipy.ndimage import shift as shiftImage

# =============================================================================
# FUNCTIONS
# =============================================================================


def displaysEqualizationHistograms(I_histogram, lower_threshold, outputFileName, log1, verbose=False):
    # hist1_before, hist1_after,hist2_before, hist2_after, , vebose=False, fileName='test'):
    fig = plt.figure(figsize=(6, 3))
    ax1 = plt.subplot(2, 2, 1)
    ax2 = plt.subplot(2, 2, 2)
    ax3 = plt.subplot(2, 2, 3)
    ax4 = plt.subplot(2, 2, 4)

    ax1.plot(I_histogram["Im1"][0][1], I_histogram["Im1"][0][0])
    ax2.plot(I_histogram["Im2"][0][1], I_histogram["Im2"][0][0])
    ax3.plot(I_histogram["Im1"][1][1], I_histogram["Im1"][1][0])
    ax4.plot(I_histogram["Im2"][1][1], I_histogram["Im2"][1][0])
    ax3.set_yscale("log")
    ax4.set_yscale("log")
    ax1.vlines(lower_threshold["Im1"], 0, I_histogram["Im1"][0][0].max(), colors="r")
    ax2.vlines(lower_threshold["Im2"], 0, I_histogram["Im2"][0][0].max(), colors="r")
    plt.savefig(outputFileName + "_intensityHist.png")
    writeString2File(
        log1.fileNameMD, "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_intensityHist.png"), "a",
    )

    if not verbose:
        plt.close(fig)


def showCCimage(image1_uncorrected, image2_uncorrected, outputFileName, shift, verbose=False):
    image_product = np.fft.fft2(image1_uncorrected) * np.fft.fft2(image2_uncorrected).conj()
    cc_image = _upsampled_dft(image_product, 150, 100, (shift * 100) + 75).conj()
    if verbose:
        plt.figure(figsize=(60, 30))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
        RGB_falsecolor_image = np.dstack([image1_uncorrected, image2_uncorrected, np.zeros([2048, 2048])])
        ax1.imshow(RGB_falsecolor_image, origin="lower", interpolation="nearest")
        ax1.set_axis_off()
        ax1.set_title("Super-imposed images")

        ax2.imshow(cc_image.real)
        ax2.set_axis_off()
        ax2.set_title("Supersampled XC sub-area")
    else:
        plt.figure(figsize=(30, 30))
        plt.imsave(outputFileName + "_CC.png", cc_image)
        plt.close()


def saveImageAdjusted(fileName, fileNameMD, image1):
    plt.figure(figsize=(30, 30))
    plt.imsave(fileName + "_adjusted.png", image1, cmap="hot")
    writeString2File(
        fileNameMD, "{}\n ![]({})\n".format(os.path.basename(fileName), fileName + "_adjusted.png"), "a",
    )
    plt.close()


def align2Files(fileName, imReference, param, log1, session1, dataFolder, verbose):
    """
    Uses preloaded ImReference Object and aligns it against filename

    Parameters
    ----------
    fileName : npy 2D array
        file of image to be aligned
    imReference : Image Class
        Object type <Image> with image reference
    param : Parameters Class
        Running parameters
    log1 : Log Class
        Logging object
    session1 : session Class
        info for session update
    dataFolder : folders Class
        DESCRIPTION.
    verbose : boolean
        True for display images

    Returns
    -------
    shift : float list, 2 dimensions
        offset in Y and X
    tableEntry : Table Class
        results zipped in Table Class form

    """
    fileName1 = imReference.fileName
    fileName2 = fileName

    outputFileName = dataFolder.outputFolders["alignImages"] + os.sep + os.path.basename(fileName2).split(".")[0]

    # loads image
    Im2 = Image()
    Im2.loadImage2D(fileName2, log1, dataFolder.outputFolders["zProject"])
    image1_uncorrected = imReference.data_2D / imReference.data_2D.max()
    image2_uncorrected = Im2.data_2D / Im2.data_2D.max()

    (
        shift,
        error,
        diffphase,
        lower_threshold,
        I_histogram,
        image2_corrected,
        image1_adjusted,
        image2_adjusted,
    ) = align2ImagesCrossCorrelation(image1_uncorrected, image2_uncorrected)

    image2_corrected_raw = shiftImage(image2_uncorrected, shift)
    image2_corrected_raw[image2_corrected_raw < 0] = 0

    # saveImageAdjusted(outputFileName+'_ref',log1.fileNameMD, image1_adjusted)
    # saveImageAdjusted(outputFileName,log1.fileNameMD, image2_adjusted)

    log1.report("Lower threshold for {}: {:.2f}".format(os.path.basename(outputFileName + "_ref"), lower_threshold["Im1"]))
    log1.report("Lower threshold for {}: {:.2f}".format(os.path.basename(outputFileName), lower_threshold["Im2"]))

    log1.report(f"Detected subpixel offset (y, x): {shift} px\n")

    # displays intensity histograms
    displaysEqualizationHistograms(I_histogram, lower_threshold, outputFileName, log1, verbose)

    """
    # saves uncrrected images to file
    save2imagesRGB(image1_uncorrected,image2_uncorrected,
                   outputFileName+'_overlay_uncorrected.png')
    """

    # thresholds corrected images for better display and saves
    # image1_corrected=image1_adjusted>0.1
    image2_corrected = image2_corrected > 0.1
    image1_uncorrected[image1_uncorrected < 0] = 0
    save2imagesRGB(
        image1_uncorrected, image2_corrected_raw, outputFileName + "_overlay_corrected.png",
    )
    writeString2File(
        log1.fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_overlay_corrected.png"),
        "a",
    )

    alignmentOutput = dataFolder.outputFiles["alignImages"]
    list2output = "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
        os.path.basename(fileName2), os.path.basename(fileName1), shift[0], shift[1], error, diffphase,
    )
    writeString2File(alignmentOutput + ".bed", list2output, "a")

    tableEntry = [
        os.path.basename(fileName2),
        os.path.basename(fileName1),
        shift[0],
        shift[1],
        error,
        diffphase,
    ]

    # saves raw registered image of fiducial
    saveImage2Dcmd(image2_corrected_raw, outputFileName + "_2d_registered", log1)

    del Im2
    return shift, tableEntry


def alignImages(param, log1, session1):

    if param.param["alignImages"]["operation"] == "overwrite":

        verbose = False

        # session
        sessionName = "alignImages"

        # processes folders and files
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.setsFolders()
        log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        log1.report("folders read: {}".format(len(dataFolder.listFolders)))
        writeString2File(
            log1.fileNameMD, "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
        )  # initialises MD file

        alignmentResultsTable = Table(
            names=("aligned file", "reference file", "shift_x", "shift_y", "error", "diffphase",),
            dtype=("S2", "S2", "f4", "f4", "f4", "f4"),
        )

        for currentFolder in dataFolder.listFolders:
            # initializes variables
            filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
            dataFolder.createsFolders(currentFolder, param)
            dictShifts = {}  # defaultdict(dict) # contains dictionary of shifts for each folder

            # generates lists of files to process
            param.files2Process(filesFolder)
            log1.report("-------> Processing Folder: {}".format(currentFolder))
            log1.report("About to process {} files\n".format(len(param.fileList2Process)))
            writeString2File(
                dataFolder.outputFiles["alignImages"] + ".bed",
                "File1 \t File_reference \t shift_y \t shift_x \t error \t diffphase",
                "w",
            )

            # Finds and loads Reference fiducial
            positionROIinformation = param.param["acquisition"]["positionROIinformation"]
            referenceBarcode = param.param["alignImages"]["referenceFiducial"]
            fileNameReferenceList, ROIList = RT2fileName(param.fileList2Process, referenceBarcode, positionROIinformation)
            
            if len(fileNameReferenceList) > 0:
                # loops over reference fiducials one ROI at a time
                for fileNameReference in fileNameReferenceList:
                    ROI = ROIList[fileNameReference]
                    imReference = Image()
                    imReference.loadImage2D(fileNameReference, log1, dataFolder.outputFolders["zProject"])
                    # saves reference 2D image for later segmentation
                    imReference.saveImage2D(
                        log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
                    )
                    dictShiftROI = {}

                    # loops over files in file list
                    for fileName in param.fileList2Process:
                        # excludes the reference fiducial and processes files in the same ROI
                        label = os.path.basename(fileName).split("_")[2]
                        if (fileName not in fileNameReference) and os.path.basename(fileName).split("_")[
                            positionROIinformation
                        ] == ROI:
                            shift, tableEntry = align2Files(fileName, imReference, param, log1, session1, dataFolder, verbose,)
                            session1.add(fileName, sessionName)
                            dictShiftROI[label] = shift.tolist()
                            alignmentResultsTable.add_row(tableEntry)
                    dictShifts["ROI:" + ROI] = dictShiftROI
                    del imReference

                # saves dicShifts
                saveJSON(dataFolder.outputFiles["dictShifts"] + ".json", dictShifts)
            else:
                log1.report(
                    "Reference Barcode file does not exist: {}", format(referenceBarcode),
                )

        alignmentResultsTable.write(
            dataFolder.outputFiles["alignImages"].split(".")[0] + ".dat", format="ascii.ecsv", overwrite=True,
        )
        del dataFolder, dictShifts


def appliesRegistrations(param, log1, session1):
    """This function will 
    - load DAPI, RNA and barcode 2D projected images, 
    - apply registrations
    - save registered images as npy arrays 
    """

    # verbose=False
    if param.param["alignImages"]["operation"] == "overwrite":

        # session
        sessionName = "registersImages"

        # processes folders and files
        dataFolder = folders(param.param["rootFolder"])
        dataFolder.setsFolders()
        log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        log1.report("folders read: {}".format(len(dataFolder.listFolders)))
        # writeString2File(log1.fileNameMD,"## {}: {}\n".format(sessionName,param.param['acquisition']['label']),'a') # initialises MD file

        positionROIinformation = param.param["acquisition"]["positionROIinformation"]

        for currentFolder in dataFolder.listFolders:
            # currentFolder=dataFolder.listFolders[0] # only one folder processed so far...
            filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
            dataFolder.createsFolders(currentFolder, param)
            log1.report("-------> Processing Folder: {}".format(currentFolder))

            # loads dicShifts with shifts for all ROIs and all labels
            dictShifts = loadJSON(dataFolder.outputFiles["dictShifts"] + ".json")

            # generates lists of files to process
            param.files2Process(filesFolder)
            log1.report("About to process {} files\n".format(len(param.fileList2Process)))

            if len(param.fileList2Process) > 0:
                # loops over files in file list
                for fileName in param.fileList2Process:

                    # gets shift from dictionary
                    ROI = os.path.basename(fileName).split("_")[positionROIinformation]
                    label = os.path.basename(fileName).split("_")[2]
                    try:
                        shiftArray = dictShifts["ROI:" + ROI][label]
                    except KeyError:
                        shiftArray = None
                        log1.report(
                            "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(ROI, label),
                            "ERROR",
                        )

                    if shiftArray != None:

                        shift = np.asarray(shiftArray)
                        # loads 2D image and applies registration
                        Im = Image()
                        Im.loadImage2D(fileName, log1, dataFolder.outputFolders["zProject"])
                        Im.data_2D = shiftImage(Im.data_2D, shift)
                        log1.report(
                            "Image registered using ROI:{}, label:{}, shift={}".format(ROI, label, shift), "info",
                        )

                        # saves registered 2D image
                        Im.saveImage2D(
                            log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
                        )

                        # logs output
                        session1.add(fileName, sessionName)
                    elif shiftArray == None and label == param.param["alignImages"]["referenceFiducial"]:
                        Im = Image()
                        Im.loadImage2D(fileName, log1, dataFolder.outputFolders["zProject"])
                        Im.saveImage2D(
                            log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
                        )
                        log1.report(
                            "Saving image for referenceRT ROI:{}, label:{}".format(ROI, label), "Warning",
                        )

                    else:
                        log1.report(
                            "No shift found in dictionary for ROI:{}, label:{}".format(ROI, label), "Warning",
                        )

        del dataFolder
