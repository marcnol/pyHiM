#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:32 2020

@author: marcnol


Purpose: Loads masks and performs a local drift correction in the region surrounding the mask.

steps for each ROI:
    - load reference fiducial 2D projection
    - iterate over cycles
    - load fiducial for cycle <i>
    - check if masks are available
    - load mask file
    - iterate over masks <j>
    - obtain subvolume from reference and cycle <i> fiducial images for mask <j>
    - cross-correlate and obtain a second-order correction value
    - determine if we keep or not
    - store in database.


"""

""
# =============================================================================
# IMPORTS
# =============================================================================


import glob, os
import matplotlib.pylab as plt
import numpy as np

from astropy.table import Table
from tqdm import trange
from skimage.util import montage
import cv2

from dask.distributed import Client, LocalCluster, get_client, as_completed, wait


from imageProcessing.imageProcessing import Image
from fileProcessing.fileManagement import folders, writeString2File, ROI2FiducialFileName
from fileProcessing.fileManagement import daskCluster

from imageProcessing.alignImages import align2ImagesCrossCorrelation

from stardist import random_label_cmap

np.random.seed(6)
lbl_cmap = random_label_cmap()


def loadsFiducial(param, fileName, log1, dataFolder):
    """
    Finds the filename of the reference fiducial to use given a DAPI image file

    Parameters
    ----------
    param : Parameters Class
        DESCRIPTION.
    fileName : string
        Full Name of file

    Returns
    -------
    TYPE string
        ROI
    TYPE Image Class
        class with the 2D projection of the reference fiducial

    """
    # finds name of reference fiducial file
    # positionROIinformation = param.param["acquisition"]["positionROIinformation"]
    # ROI = os.path.basename(fileName).split("_")[positionROIinformation]
    ROI = param.decodesFileParts(os.path.basename(fileName))["roi"]

    referenceBarcode = param.param["alignImages"]["referenceFiducial"]
    fiducialFilename = ROI2FiducialFileName(param, fileName, referenceBarcode)

    # loads reference fiducial file
    if len(fiducialFilename) < 1:
        print(
            "Error, no reference candidate found for ROI:{} | filename:{}\n\n**Check you are using the correct reference fiducial!\n".format(
                ROI, os.path.basename(fileName)
            )
        )
        return -1
    elif len(fiducialFilename) > 1:
        print(
            "Error, too many reference candidates found for ROI:{} | filename:{}\nCandidates detected:{}".format(
                ROI, os.path.basename(fileName), fiducialFilename
            )
        )
        return -1
    else:
        print("Using reference fiducial> {}\n".format(os.path.basename(fiducialFilename[0])))
        # fullFiducialFilename = currentFolder+os.sep+fiducialFilename[0]
        fullFiducialFilename = fiducialFilename[0]

    imReference = Image(param, log1)
    imReference.loadImage2D(fullFiducialFilename, log1, dataFolder.outputFolders["zProject"])
    # imReference.imageShow(show=True)

    return ROI, imReference


def retrieveBarcodeList(param, fileName):
    """
    retrieves list of barcodes for which a fiducial is available in this ROI

    Parameters
    ----------
    param : Parameters Class
        DESCRIPTION.
    fileName : string
        Full Name of file

    Returns
    -------
    barcodeList : list
        list of barcodes retrieved in rootFolder.
    fiducialFileNames : list
        list of fiducialFileNames retrieved in rootFolder.
    """
    rootFolder = os.path.dirname(fileName)
    # positionROIinformation = param.param["acquisition"]["positionROIinformation"]
    # ROI = os.path.basename(fileName).split("_")[positionROIinformation]
    ROI = param.decodesFileParts(os.path.basename(fileName))["roi"]

    channelFiducial = param.param["acquisition"]["fiducialBarcode_channel"]

    listFiles = glob.glob(rootFolder + os.sep + "*.tif")

    # if (ROI in os.path.basename(x).split("_")[positionROIinformation])
    fiducialFileNames = [
        x
        for x in listFiles
        if (ROI in param.decodesFileParts(os.path.basename(x))["roi"])
        and ("RT" in os.path.basename(x))
        and (channelFiducial in os.path.basename(x))
    ]

    barcodeList = [os.path.basename(x).split("_")[2] for x in fiducialFileNames]

    print("barcodes to process: {}".format(barcodeList))

    return barcodeList, fiducialFileNames


def alignsSubVolumes(imageReference, imageBarcode, Masks, bezel=20, iMask=1):
    maskSize = Masks.shape

    # defines sub-volume around Mask

    minx, miny = np.nonzero(Masks == iMask)[0].min(), np.nonzero(Masks == iMask)[1].min()
    maxx, maxy = np.nonzero(Masks == iMask)[0].max(), np.nonzero(Masks == iMask)[1].max()

    minx, miny = np.max([minx - bezel, 0]), np.max([miny - bezel, 0])
    maxx, maxy = np.min([maxx + bezel, maskSize[0]]), np.min([maxy + bezel, maskSize[1]])
    boundingBox = np.array([minx, maxx, miny, maxy])
    del minx, maxx, miny, maxy

    # obtain subvolume from reference and cycle <i> fiducial images for mask <iMask>
    subVolumeReference = imageReference[boundingBox[0] : boundingBox[1], boundingBox[2] : boundingBox[3]]
    subVolume = imageBarcode[boundingBox[0] : boundingBox[1], boundingBox[2] : boundingBox[3]]

    # adjusts levels
    subVolume = subVolume / subVolume.max()
    subVolumeReference = subVolumeReference / subVolumeReference.max()

    # calculates shift using cross-correlation
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
        subVolumeReference, subVolume, lower_threshold=0.1, higher_threshold=0.9999999, upsample_factor=10
    )

    # print("Shift = {}".format(shift))
    # subVolumeCorrected = shiftImage(image2_adjusted, shift)

    result = shift, error, diffphase, image1_adjusted, image2_adjusted, image2_corrected

    return result


def pad_images_to_same_size(images):
    """
    :param images: sequence of images
    :return: list of images padded so that all images have same width and height (max width and height are used)
    """
    width_max = 0
    height_max = 0
    for img in images:
        h, w = img.shape[:2]
        width_max = max(width_max, w)
        height_max = max(height_max, h)

    images_padded = []
    for img in images:
        h, w = img.shape[:2]
        diff_vert = height_max - h
        pad_top = diff_vert // 2
        pad_bottom = diff_vert - pad_top
        diff_hori = width_max - w
        pad_left = diff_hori // 2
        pad_right = diff_hori - pad_left
        img_padded = cv2.copyMakeBorder(img, pad_top, pad_bottom, pad_left, pad_right, cv2.BORDER_CONSTANT, value=0)
        assert img_padded.shape[:2] == (height_max, width_max)
        images_padded.append(img_padded)

    return images_padded


def plotMontageImage(
    montage2DReference, montage2DCorrected, outputFileName, fileNameMD, fileInformation, tag="_Corrected.png"
):
    fig, ax1 = plt.subplots()
    fig.set_size_inches((30, 30))

    montage2DReference = montage2DReference / montage2DReference.max()
    montage2DCorrected = montage2DCorrected / montage2DCorrected.max()

    sz = montage2DReference.shape

    nullImage = np.zeros(sz)
    RGB = np.dstack([montage2DReference, montage2DCorrected, nullImage])

    ax1.imshow(RGB)
    ax1.axis("off")
    fig.savefig(outputFileName + tag)
    plt.close(fig)

    writeString2File(fileNameMD, fileInformation + "\n![]({})\n".format(outputFileName + tag), "a")


def localDriftCorrection_plotsLocalAlignments(
    imageListCorrected, imageListunCorrected, imageListReference, log1, dataFolder, ROI, barcode
):
    """
    converts list of images into mosaic and saves results

    Parameters
    ----------
    imageListCorrected : TYPE
        DESCRIPTION.
    imageListunCorrected : TYPE
        DESCRIPTION.
    imageListReference : TYPE
        DESCRIPTION.
    outputFileName : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    outputFileName = (
        dataFolder.outputFolders["alignImages"] + os.sep + "localDriftCorrection_" + "ROI:" + ROI + "barcode:" + barcode
    )

    imageListCorrectedPadded = pad_images_to_same_size(imageListCorrected)
    imageListunCorrectedPadded = pad_images_to_same_size(imageListunCorrected)
    imageListReferencePadded = pad_images_to_same_size(imageListReference)

    montage2DReference = montage(imageListReferencePadded)
    montage2DCorrected = montage(imageListCorrectedPadded)
    montage2DunCorrected = montage(imageListunCorrectedPadded)

    fileInformation = "**uncorrected** drift for ROI: {} barcode:{}".format(ROI, barcode)
    plotMontageImage(
        montage2DReference,
        montage2DunCorrected,
        outputFileName,
        log1.fileNameMD,
        fileInformation,
        tag="_uncorrected.png",
    )

    fileInformation = "**corrected** drift for ROI: {} barcode:{}".format(ROI, barcode)
    plotMontageImage(
        montage2DReference, montage2DCorrected, outputFileName, log1.fileNameMD, fileInformation, tag="_corrected.png"
    )

    del montage2DReference, montage2DCorrected, imageListReferencePadded, imageListCorrectedPadded
    del imageListReference, imageListCorrected


def localDriftCorrection_savesResults(dictShift, alignmentResultsTable, dataFolder, log1):

    # plots shift violins
    for ROI in dictShift.keys():
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
        shifts, labels = [], []

        for barcode in dictShift[ROI].keys():
            shiftsBarcode = np.array([0, 0])
            labels.append(barcode)
            for iMask in dictShift[ROI][barcode].keys():
                shiftsBarcode = np.concatenate((shiftsBarcode, dictShift[ROI][barcode][iMask]))
            shifts.append(shiftsBarcode)

        plt.violinplot(shifts)

        ax1.get_xaxis().set_tick_params(direction="out")
        ax1.xaxis.set_ticks_position("bottom")
        ax1.set_xticks(np.arange(1, len(labels) + 1))
        ax1.set_xticklabels(labels)
        ax1.set_xlim(0.25, len(labels) + 0.75)
        ax1.set_xlabel("ROI #{}".format(ROI))
        ax1.set_ylabel("local drift shift, px")
        outputFileName = dataFolder.outputFolders["alignImages"] + os.sep + "LocalShiftsViolinPlot_" + "ROI:" + ROI
        plt.savefig(outputFileName + ".png")
        plt.close()
        writeString2File(
            log1.fileNameMD, "Local drift for ROI: {}\n ![]({})\n".format(ROI, outputFileName + ".png"), "a",
        )

    # saves Table with all shifts
    alignmentResultsTable.write(
        dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat",
        format="ascii.ecsv",
        overwrite=True,
    )


def localDriftforRT(
    barcode,
    fileNameFiducial,
    imReferenceFileName,
    imageReferenceBackgroundSubstracted,
    Masks,
    bezel,
    shiftTolerance,
    ROI,
    alignmentResultsTable,
    log1,
    dataFolder,
    parallel=False,
):

    # loads 2D image and applies registration
    param = {}
    Im = Image(param, log1)
    Im.loadImage2D(fileNameFiducial, log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered")
    imageBarcode = Im.removesBackground2D(normalize=True)

    imageListCorrected, imageListunCorrected, imageListReference, errormessage = [], [], [], []

    dictShiftBarcode = {}

    # - iterate over masks <iMask>
    log1.info("Looping over {} masks for barcode: {}\n".format(Masks.max(), barcode))
    if parallel:
        Maskrange = range(1, Masks.max())
        log1.report("See progress in http://localhost:8787 ")
    else:
        Maskrange = trange(1, Masks.max())

    parallel = False

    if parallel:

        # need to fix the memory leakage!!
        # thre problem is that each mask submits a worker (500 masks) and they all return
        # 3 subVolumes that need to be all accumulated by a single worker (that is handling the barcode).
        # To solve the issue I need to recode so that images are not returned by the workers.
        # The issue now is that CURRENTLY these images need to be collected in order to make the mosaic...

        log1.info("Launching {} threads using dask".format(len(Maskrange)))
        futures = []
        client = get_client()

        for iMask in Maskrange:
            # calculates shift
            futures.append(
                client.submit(
                    alignsSubVolumes, imageReferenceBackgroundSubstracted, imageBarcode, Masks, bezel=bezel, iMask=iMask
                )
            )

        log1.info("Waiting for {} results to arrive".format(len(futures)))

        for batch in as_completed(futures, with_results=True).batches():
            log1.info("Reading batch from {} workers".format(len(batch)))
            for future, iResult in batch:

                shift, error, diffphase, subVolumeReference, subVolume, subVolumeCorrected = iResult
                del iResult

                # stores images in list
                imageListReference.append(subVolumeReference)
                imageListunCorrected.append(subVolume)

                # evaluates shift to determine if we keep or not
                if np.nonzero(np.absolute(shift) > shiftTolerance)[0].shape[0] > 0:
                    errormessage.append(
                        "ROI:{} | barcode:{}| Mask:{}> local shift = {} not kept as it is over the tolerance of {} px".format(
                            ROI, barcode, iMask, shift, shiftTolerance
                        )
                    )

                    shift = np.array([0, 0])
                    imageListCorrected.append(subVolume)
                else:
                    imageListCorrected.append(subVolumeCorrected)

                del subVolume, subVolumeCorrected, subVolumeReference
                # stores result in database
                dictShiftBarcode[str(iMask)] = shift

                # creates Table entry to return
                tableEntry = [
                    os.path.basename(fileNameFiducial),
                    os.path.basename(imReferenceFileName),
                    int(ROI),
                    int(barcode.split("RT")[1]),
                    iMask,
                    shift[0],
                    shift[1],
                    error,
                    diffphase,
                ]

                alignmentResultsTable.add_row(tableEntry)

        del futures

        resultAll = [dictShiftBarcode, imageListCorrected, imageListunCorrected, imageListReference, errormessage]

    else:

        for iMask in Maskrange:
            # calculates shift

            result = alignsSubVolumes(
                imageReferenceBackgroundSubstracted, imageBarcode, Masks, bezel=bezel, iMask=iMask
            )
            shift, error, diffphase, subVolumeReference, subVolume, subVolumeCorrected = result

            # stores images in list
            imageListReference.append(subVolumeReference)
            imageListunCorrected.append(subVolume)

            # evaluates shift to determine if we keep or not
            if np.nonzero(np.absolute(shift) > shiftTolerance)[0].shape[0] > 0:
                errormessage.append(
                    "ROI:{} | barcode:{}| Mask:{}> local shift = {} not kept as it is over the tolerance of {} px".format(
                        ROI, barcode, iMask, shift, shiftTolerance
                    )
                )

                shift = np.array([0, 0])
                imageListCorrected.append(subVolume)
            else:
                imageListCorrected.append(subVolumeCorrected)

            # stores result in database
            dictShiftBarcode[str(iMask)] = shift

            # creates Table entry to return
            tableEntry = [
                os.path.basename(fileNameFiducial),
                os.path.basename(imReferenceFileName),
                int(ROI),
                int(barcode.split("RT")[1]),
                iMask,
                shift[0],
                shift[1],
                error,
                diffphase,
            ]

            alignmentResultsTable.add_row(tableEntry)

        resultAll = [dictShiftBarcode, imageListCorrected, imageListunCorrected, imageListReference, errormessage]

    return resultAll


def localDriftallBarcodes(
    param,
    log1,
    dataFolder,
    fileNameDAPI,
    localDriftforRT,
    imReferenceFileName,
    imageReferenceBackgroundSubstracted,
    Masks,
    bezel,
    shiftTolerance,
    ROI,
    alignmentResultsTable,
):

    # - retrieves list of barcodes for which a fiducial is available in this ROI
    barcodeList, fiducialFileNames = retrieveBarcodeList(param, fileNameDAPI)
    if imReferenceFileName in fiducialFileNames:
        fiducialFileNames.remove(imReferenceFileName)

    # print(">>>Image <i> {} | sumImage: {}".format(imReference.data_2D.max(), imageReference.max()))

    dictShift, errormessage = {}, []

    if param.param["parallel"]:

        futures = list()
        client = get_client()

        remote_imReference = client.scatter(imageReferenceBackgroundSubstracted, broadcast=True)
        remote_Masks = client.scatter(Masks, broadcast=True)

        for barcode, fileNameFiducial in zip(barcodeList, fiducialFileNames):

            futures.append(
                client.submit(
                    localDriftforRT,
                    barcode,
                    fileNameFiducial,
                    imReferenceFileName,
                    remote_imReference,
                    remote_Masks,
                    bezel,
                    shiftTolerance,
                    ROI,
                    alignmentResultsTable,
                    log1,
                    dataFolder,
                    parallel=True,
                )
            )

        wait(futures)

        results = client.gather(futures)

        del remote_imReference, remote_Masks, futures

        print("Retrieving {} results from cluster".format(len(results)))
        for result, barcode in zip(results, barcodeList):
            dictShift[barcode], imageListCorrected, imageListunCorrected, imageListReference, errormessage1 = result
            errormessage += errormessage1
            # output mosaics with global and local alignments
            localDriftCorrection_plotsLocalAlignments(
                imageListCorrected, imageListunCorrected, imageListReference, log1, dataFolder, ROI, barcode
            )
    else:
        # - load fiducial for cycle <i>
        for barcode, fileNameFiducial in zip(barcodeList, fiducialFileNames):

            # calculates local drift for barcode by looping over Masks
            result = localDriftforRT(
                barcode,
                fileNameFiducial,
                imReferenceFileName,
                imageReferenceBackgroundSubstracted,
                Masks,
                bezel,
                shiftTolerance,
                ROI,
                alignmentResultsTable,
                log1,
                dataFolder,
            )

            dictShift[barcode], imageListCorrected, imageListunCorrected, imageListReference, errormessage1 = result
            errormessage += errormessage1

            # output mosaics with global and local alignments
            localDriftCorrection_plotsLocalAlignments(
                imageListCorrected, imageListunCorrected, imageListReference, log1, dataFolder, ROI, barcode
            )

    return dictShift, alignmentResultsTable, errormessage


def localDriftCorrection(param, log1, session1):
    sessionName = "localDriftCorrection"
    # tqdm._instances.clear()

    # processes folders and files
    log1.addSimpleText(
        "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
    )
    dataFolder = folders(param.param["rootFolder"])
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(session1.name, param.param["acquisition"]["label"]), "a",
    )

    if "localShiftTolerance" in param.param["alignImages"].keys():
        shiftTolerance = param.param["alignImages"]["localShiftTolerance"]
    else:
        shiftTolerance = 1

    if "bezel" in param.param["alignImages"].keys():
        bezel = param.param["alignImages"]["bezel"]
    else:
        bezel = 20

    print("Parameters> shiftTolerance = {} | bezel = {}".format(shiftTolerance, bezel))
    errormessage = []

    alignmentResultsTable = Table(
        names=(
            "reference file",
            "aligned file",
            "ROI #",
            "Barcode #",
            "CellID #",
            "shift_x",
            "shift_y",
            "error",
            "diffphase",
        ),
        dtype=("S2", "S2", "int", "int", "int", "f4", "f4", "f4", "f4"),
    )
    dictShift = {}

    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        dictShift = {}

        # iterates over ROIs
        for fileNameDAPI in param.fileList2Process:

            # label = param.param["acquisition"]["label"]

            # appends filename to session and defines ROI
            session1.add(fileNameDAPI, sessionName)
            print("processing> {}\n".format(os.path.basename(fileNameDAPI)))

            # - loads masks for ROI
            fileNameROImasks = os.path.basename(fileNameDAPI).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = os.sep.join(
                [os.path.dirname(fileNameDAPI), param.param["segmentedObjects"]["folder"], fileNameROImasks]
            )

            if os.path.exists(fullFileNameROImasks):
                Masks = np.load(fullFileNameROImasks)
                # fig = plt.figure(), plt.imshow(Masks, origin="lower", cmap=lbl_cmap,alpha=1)
                print("Masks read> {}".format(Masks.max()))

            else:
                print("Error> No Mask found! File expected: {}".format(fullFileNameROImasks))
                raise FileNotFoundError("I cannot find DAPI mask: {}".format(fullFileNameROImasks))

            # - loads reference fiducial file
            ROI, imReference = loadsFiducial(param, fileNameDAPI, log1, dataFolder)
            imageReferenceBackgroundSubstracted = imReference.removesBackground2D(normalize=True)

            dictShift[ROI], alignmentResultsTable, errormessage1 = localDriftallBarcodes(
                param,
                log1,
                dataFolder,
                fileNameDAPI,
                localDriftforRT,
                imReference.fileName,
                imageReferenceBackgroundSubstracted,
                Masks,
                bezel,
                shiftTolerance,
                ROI,
                alignmentResultsTable,
            )

            errormessage = errormessage + errormessage1
            # produces shift violin plots and saves results Table
            localDriftCorrection_savesResults(dictShift, alignmentResultsTable, dataFolder, log1)
            log1.report("\n".join(errormessage))

    return 0, dictShift, alignmentResultsTable
