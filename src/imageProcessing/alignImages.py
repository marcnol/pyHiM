#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 09:22:18 2020

@author: marcnol

Sets of functions that do alignment of 2D fiducial images. It also contains
code to apply these alignments to other channels (masks/ barcodes)

For the time being alignment is purely based on optimized sub-pixed accuracy
image cross correlation

"""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import os, glob
from dask.distributed import get_client

from skimage.registration._phase_cross_correlation import _upsampled_dft
from skimage.exposure import match_histograms

from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

from imageProcessing.imageProcessing import (
    Image,
    save2imagesRGB,
    saveImage2Dcmd,
    saveImageDifferences,
    align2ImagesCrossCorrelation,
    alignImagesByBlocks,
    plottingBlockALignmentResults,
)

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    saveJSON,
    loadJSON,
    RT2fileName,
    printLog,
    getDictionaryValue,
)

from astropy.table import Table
from scipy.ndimage import shift as shiftImage

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")
# =============================================================================
# FUNCTIONS
# =============================================================================


def displaysEqualizationHistograms(I_histogram, lower_threshold, outputFileName, fileNameMD, verbose=False):
    # hist1_before, hist1_after,hist2_before, hist2_after, , vebose=False, fileName='test'):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    # fig= plt.figure(figsize=(6, 3))

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
        fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_intensityHist.png"),
        "a",
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


def removesInhomogeneousBackground(im, param):

    sigma_clip = SigmaClip(sigma=param.param["alignImages"]["background_sigma"])
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)

    im1_bkg_substracted = im - bkg.background

    return im1_bkg_substracted


def align2Files(fileName, imReference, param, session1, dataFolder, verbose):
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
    session1 : session Class
        info for session update
    dataFolder : folders Class
        DESCRIPTION.
    verbose : boolean
        True for display images

    Returns are returned as arguments!
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
    Im2 = Image(param)
    Im2.loadImage2D(fileName2, dataFolder.outputFolders["zProject"])

    # Normalises images
    image1_uncorrected = imReference.data_2D / imReference.data_2D.max()
    image2_uncorrected = Im2.data_2D / Im2.data_2D.max()

    # removes inhomogeneous background
    image1_uncorrected = removesInhomogeneousBackground(image1_uncorrected, param)
    image2_uncorrected = removesInhomogeneousBackground(image2_uncorrected, param)

    lower_threshold = getDictionaryValue(param.param["alignImages"], "lower_threshold", default=0.999)
    higher_threshold = getDictionaryValue(param.param["alignImages"], "higher_threshold", default=0.9999999)
    alignByBlock = getDictionaryValue(param.param["alignImages"], "alignByBlock", default=False)
    tolerance = getDictionaryValue(param.param["alignImages"], "tolerance", default=0.1)
    blockSize = getDictionaryValue(param.param["alignImages"], "blockSize", default=256)

    # if "lower_threshold" in param.param["alignImages"].keys():
    #     lower_threshold = param.param["alignImages"]["lower_threshold"]
    # else:
    #     lower_threshold = 0.999

    # if "higher_threshold" in param.param["alignImages"].keys():
    #     higher_threshold = param.param["alignImages"]["higher_threshold"]
    # else:
    #     higher_threshold = 0.9999999

    # if "alignByBlock" in param.param["alignImages"].keys():
    #     alignByBlock = param.param["alignImages"]["alignByBlock"]
    # else:
    #     alignByBlock = False

    # if "tolerance" in param.param["alignImages"].keys():
    #     tolerance = param.param["alignImages"]["tolerance"]
    # else:
    #     tolerance = 0.1

    # if "blockSize" in param.param["alignImages"].keys():
    #     blockSize = param.param["alignImages"]["blockSize"]
    # else:
    #     blockSize = 256

    if not alignByBlock:
        # [calculates unique translation for the entire image using cross-correlation]
        (
            shift,
            error,
            diffphase,
            lower_threshold,
            I_histogram,
            image2_corrected,
            image1_adjusted,
            image2_adjusted,
        ) = align2ImagesCrossCorrelation(
            image1_uncorrected, image2_uncorrected, lower_threshold=lower_threshold, higher_threshold=higher_threshold
        )

        # displays intensity histograms
        displaysEqualizationHistograms(I_histogram, lower_threshold, outputFileName, param.param["fileNameMD"],verbose)

    else:
        # [calculates block translations by cross-correlation and gets overall shift by polling]

        # normalizes images
        image1_uncorrected, image2_uncorrected = np.float32(image1_uncorrected), np.float32(image2_uncorrected)

        # matches histograms
        image2_uncorrected = np.float32(match_histograms(image2_uncorrected, image1_uncorrected))

        # calculates block shifts and polls for most favourable shift
        upsample_factor = 100
        blockSize = (blockSize, blockSize)

        (shift, error, relativeShifts, rmsImage, contour,) = alignImagesByBlocks(
            image1_uncorrected,
            image2_uncorrected,
            blockSize,
            upsample_factor=upsample_factor,
            minNumberPollsters=4,
            tolerance=tolerance,
        )
        diffphase = 0

        plottingBlockALignmentResults(
            relativeShifts, rmsImage, contour, fileName=outputFileName + "_block_alignments.png"
        )

        writeString2File(
            param.param["fileNameMD"],
            "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_block_alignments.png"),
            "a",
        )

        # saves mask of valid regions with a correction within the tolerance
        saveImage2Dcmd(rmsImage, outputFileName + "_rmsBlockMap")
        saveImage2Dcmd(relativeShifts, outputFileName + "_errorAlignmentBlockMap")

    image2_corrected_raw = shiftImage(image2_uncorrected, shift)

    image2_corrected_raw[image2_corrected_raw < 0] = 0

    error = np.sum(np.sum(np.abs(image1_uncorrected - image2_corrected_raw), axis=1))

    printLog(f"$ Detected subpixel offset (y, x): {shift} px")

    # [displays and saves results]

    # thresholds corrected images for better display and saves
    image1_uncorrected[image1_uncorrected < 0] = 0
    image2_uncorrected[image2_uncorrected < 0] = 0

    save2imagesRGB(image1_uncorrected, image2_corrected_raw, outputFileName + "_overlay_corrected.png")

    saveImageDifferences(
        image1_uncorrected,
        image2_uncorrected,
        image1_uncorrected,
        image2_corrected_raw,
        outputFileName + "_referenceDifference.png",
    )

    # reports image in MD file
    writeString2File(
        param.param["fileNameMD"],
        "{}\n ![]({})\n ![]({})\n".format(
            os.path.basename(outputFileName),
            outputFileName + "_overlay_corrected.png",
            outputFileName + "_referenceDifference.png",
        ),
        "a",
    )

    # outputs results to logfile
    alignmentOutput = dataFolder.outputFiles["alignImages"]
    list2output = "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
        os.path.basename(fileName2), os.path.basename(fileName1), shift[0], shift[1], error, diffphase,
    )
    writeString2File(alignmentOutput, list2output, "a")

    # creates Table entry to return
    tableEntry = [
        os.path.basename(fileName2),
        os.path.basename(fileName1),
        shift[0],
        shift[1],
        error,
        diffphase,
    ]

    # saves registered fiducial image
    saveImage2Dcmd(image2_corrected_raw, outputFileName + "_2d_registered")

    del Im2
    return shift, tableEntry


def alignImagesInCurrentFolder(currentFolder, param, dataFolder, session1, fileName=None):
    # session
    sessionName = "alignImages"
    verbose = False

    alignmentResultsTable = Table(
        names=("aligned file", "reference file", "shift_x", "shift_y", "error", "diffphase",),
        dtype=("S2", "S2", "f4", "f4", "f4", "f4"),
    )

    # initializes variables
    filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
    dataFolder.createsFolders(currentFolder, param)
    dictShifts = {}  # defaultdict(dict) # contains dictionary of shifts for each folder

    # generates lists of files to process for currentFolder
    param.files2Process(filesFolder)
    printLog("> Processing Folder: {}".format(currentFolder))
    printLog("> About to process {} files\n".format(len(param.fileList2Process)))
    writeString2File(
        dataFolder.outputFiles["alignImages"],
        "File1 \t File_reference \t shift_y \t shift_x \t error \t diffphase",
        "w",
    )

    # Finds and loads Reference fiducial information
    # positionROIinformation = param.param["acquisition"]["positionROIinformation"]
    referenceBarcode = param.param["alignImages"]["referenceFiducial"]
    printLog("$ Reference fiducial {}".format(referenceBarcode))

    # retrieves the list of fiducial image files to be aligned
    fileNameReferenceList, ROIList = RT2fileName(param, referenceBarcode)

    if len(fileNameReferenceList) > 0:

        # loops over fiducials images one ROI at a time
        for fileNameReference in fileNameReferenceList:

            # loads reference fiducial image for this ROI
            ROI = ROIList[fileNameReference]
            imReference = Image(param)
            imReference.loadImage2D(fileNameReference, dataFolder.outputFolders["zProject"])
            printLog("> Loading reference Image {}".format(fileNameReference))

            # saves reference 2D image of fiducial
            if not os.path.exists(
                imReference.getImageFileName(dataFolder.outputFolders["alignImages"], tag="_2d_registered")
            ):
                imReference.saveImage2D(
                    dataFolder.outputFolders["alignImages"], tag="_2d_registered",
                )

            dictShiftROI = {}

            fileName2ProcessList = [
                x
                for x in param.fileList2Process
                if (x not in fileNameReference) and param.decodesFileParts(os.path.basename(x))["roi"] == ROI
            ]
            printLog("Found {} files in ROI: {}".format(len(fileName2ProcessList), ROI))
            printLog(
                "[roi:cycle] {}".format(
                    "|".join(
                        [
                            str(param.decodesFileParts(os.path.basename(x))["roi"])
                            + ":"
                            + str(param.decodesFileParts(os.path.basename(x))["cycle"])
                            for x in fileName2ProcessList
                        ]
                    )
                )
            )

            if param.param["parallel"]:
                # running in parallel mode
                client = get_client()
                futures = list()
                labels = []

                for fileName2Process in fileName2ProcessList:
                    # excludes the reference fiducial and processes files in the same ROI
                    labels.append(os.path.basename(fileName2Process).split("_")[2])
                    futures.append(
                        client.submit(
                            align2Files, fileName2Process, imReference, param, session1, dataFolder, verbose
                        )
                    )

                printLog("$ Waiting for {} results to arrive".format(len(futures)))

                results = client.gather(futures)

                printLog("$ Retrieving {} results from cluster".format(len(results)))

                for result, label in zip(results, labels):
                    shift, tableEntry = result
                    dictShiftROI[label] = shift.tolist()
                    alignmentResultsTable.add_row(tableEntry)
                    session1.add(fileName2Process, sessionName)
                    # printLog("Processed: {}".format(label))
            else:
                # running in sequential mode
                nFiles=len(param.fileList2Process)

                for iFile, fileName2Process in enumerate(param.fileList2Process):
                    # excludes the reference fiducial and processes files in the same ROI
                    label = os.path.basename(fileName2Process).split("_")[2]
                    roi = param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
                    printLog("\n$ About to process file {} \ {}".format(iFile,nFiles))

                    if (fileName2Process not in fileNameReference) and roi == ROI:
                        if fileName == None or (
                            fileName != None and os.path.basename(fileName) == os.path.basename(fileName2Process)
                        ):
                            # aligns files and saves results to database in dict format and to a Table
                            shift, tableEntry = align2Files(
                                fileName2Process, imReference, param, session1, dataFolder, verbose,
                            )
                            dictShiftROI[label] = shift.tolist()
                            alignmentResultsTable.add_row(tableEntry)
                            session1.add(fileName2Process, sessionName)
                    elif (fileName2Process in fileNameReference):
                        printLog("\n$ Skipping reference file: {} ".format(os.path.basename(fileName2Process)))
            # accumulates shifst for this ROI into global dictionary
            dictShifts["ROI:" + ROI] = dictShiftROI
            del imReference

        # saves dicShifts dictionary with shift results
        dictionaryFileName = os.path.splitext(dataFolder.outputFiles["dictShifts"])[0] + ".json"
        saveJSON(dictionaryFileName, dictShifts)
        printLog("$ Saved alignment dictionary to {}".format(dictionaryFileName))

    else:
        printLog(
            "# Reference Barcode file does not exist: {}".format(referenceBarcode))
        raise ValueError

    return alignmentResultsTable


def alignImages(param, session1, fileName=None):
    """
    From a given parameters class it aligns all the fiducial images

    Parameters
    ----------
    param : Parameters class
        running parameters
    session1 : Session Class
        logs session information.

    Returns
    -------
    None.

    """
    sessionName = "registersImages"


    # processes folders and adds information to log files
    dataFolder = folders(param.param["rootFolder"])
    dataFolder.setsFolders()
    printLog("\n===================={}====================\n".format(sessionName))
    printLog("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        param.param["fileNameMD"], "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
    )

    # loops over folders
    for currentFolder in dataFolder.listFolders:
        alignmentResultsTable = alignImagesInCurrentFolder(
            currentFolder, param, dataFolder, session1, fileName
        )

    # saves Table with all shifts
    alignmentResultsTable.write(
        dataFolder.outputFiles["alignImages"].split(".")[0] + ".table", format="ascii.ecsv", overwrite=True,
    )

    del dataFolder


def appliesRegistrations2fileName(fileName2Process, param, dataFolder, session1, dictShifts):
    """
    Applies registration of fileName2Process

    Parameters
    ----------
    fileName2Process : string
        file to apply registration to
    param : Parameters class
    dataFolder : dataFolder class
    session1 : Session class
    dictShifts : Dictionnary
        contains the shifts to be applied to all ROIs

    Returns
    -------
    None.

    """
    # session
    sessionName = "registersImages"

    # gets shift from dictionary
    # ROI = os.path.basename(fileName2Process).split("_")[positionROIinformation]
    ROI = param.decodesFileParts(os.path.basename(fileName2Process))["roi"]

    label = os.path.basename(fileName2Process).split("_")[2]  # to FIX

    try:
        shiftArray = dictShifts["ROI:" + ROI][label]
    except KeyError:
        shiftArray = None
        printLog("$ Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                "ROI:" + ROI, label))

    if shiftArray != None:

        shift = np.asarray(shiftArray)
        # loads 2D image and applies registration
        Im = Image(param)
        Im.loadImage2D(fileName2Process, dataFolder.outputFolders["zProject"])
        Im.data_2D = shiftImage(Im.data_2D, shift)
        printLog("$ Image registered using ROI:{}, label:{}, shift={}".format(ROI, label, shift))

        # saves registered 2D image
        Im.saveImage2D(
            dataFolder.outputFolders["alignImages"], tag="_2d_registered",
        )

        # logs output
        session1.add(fileName2Process, sessionName)
    elif shiftArray == None and label == param.param["alignImages"]["referenceFiducial"]:
        Im = Image(param)
        Im.loadImage2D(fileName2Process, dataFolder.outputFolders["zProject"])
        Im.saveImage2D(
            dataFolder.outputFolders["alignImages"], tag="_2d_registered",
        )
        printLog(
            "$ Saving image for referenceRT ROI:{}, label:{}".format(ROI, label))

    else:
        printLog(
            "# No shift found in dictionary for ROI:{}, label:{}".format(ROI, label),status='WARN')


def appliesRegistrations2currentFolder(currentFolder, param, dataFolder, session1, fileName=None):
    """
    applies registrations to all files in currentFolder

    Parameters
    ----------
    currentFolder : TYPE
        DESCRIPTION.
    param : Parameters class
    dataFolder : dataFolder class
    session1 : Session class
    fileName : string, optional
        File to process. The default is None.

    Returns
    -------
    None.

    """

    # currentFolder=dataFolder.listFolders[0] # only one folder processed so far...
    filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
    dataFolder.createsFolders(currentFolder, param)
    printLog("> Processing Folder: {}".format(currentFolder))

    # loads dicShifts with shifts for all ROIs and all labels
    dictFileName = os.path.splitext(dataFolder.outputFiles["dictShifts"])[0] + ".json"

    # dictFileName = dataFolder.outputFiles["dictShifts"] + ".json"
    dictShifts = loadJSON(dictFileName)
    if len(dictShifts) == 0:
        printLog("# File with dictionary not found!: {}".format(dictFileName))
    else:
        printLog("$ Dictionary File loaded: {}".format(dictFileName))

    # generates lists of files to process
    param.files2Process(filesFolder)
    nFiles=len(param.fileList2Process)
    printLog("\n$ About to process {} files\n".format(nFiles))

    if len(param.fileList2Process) > 0:
        # loops over files in file list
        for i, fileName2Process in enumerate(param.fileList2Process):
            if fileName == None or (
                fileName != None and os.path.basename(fileName) == os.path.basename(fileName2Process)
            ):
                printLog("\n$ About to process file {} \ {}".format(i,nFiles))
                appliesRegistrations2fileName(fileName2Process, param, dataFolder, session1, dictShifts)


def appliesRegistrations(param, session1, fileName=None):
    """This function will
    - load masks, RNA and barcode 2D projected images,
    - apply registrations
    - save registered images as npy arrays
    """

    sessionName = "registersImages"

    # verbose=False

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    dataFolder.setsFolders()
    printLog("\n===================={}====================\n".format(sessionName))
    printLog("$ folders read: {}".format(len(dataFolder.listFolders)))

    for currentFolder in dataFolder.listFolders:
        appliesRegistrations2currentFolder(currentFolder, param, dataFolder, session1, fileName)

    del dataFolder
