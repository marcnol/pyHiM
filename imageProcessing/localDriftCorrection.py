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
from time import sleep

from astropy.table import Table
from tqdm import trange
from skimage.util import montage
import cv2

from imageProcessing.imageProcessing import Image
from fileProcessing.fileManagement import (
    folders, 
    writeString2File,
    ROI2FiducialFileName)

from imageProcessing.alignImages import align2ImagesCrossCorrelation

from stardist import random_label_cmap

np.random.seed(6)
lbl_cmap = random_label_cmap()

    
def loadsFiducial(param, fileName, log1,dataFolder):
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
    ROI = param.decodesFileParts(os.path.basename(fileName))['roi']

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
            "Error, too many reference candidates found for ROI:{} | filename:{}\n".format(
                ROI, os.path.basename(fileName)
            )
        )
        return -1
    else:
        print("Using reference fiducial> {}\n".format(os.path.basename(fiducialFilename[0])))
        # fullFiducialFilename = currentFolder+os.sep+fiducialFilename[0]
        fullFiducialFilename = fiducialFilename[0]

    imReference = Image()
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
    ROI = param.decodesFileParts(os.path.basename(fileName))['roi']
    
    channelFiducial = param.param["acquisition"]["fiducialBarcode_channel"]

    listFiles = glob.glob(rootFolder + os.sep + "*.tif")

    # if (ROI in os.path.basename(x).split("_")[positionROIinformation])
    fiducialFileNames = [
        x
        for x in listFiles
        if (ROI in param.decodesFileParts(os.path.basename(x))['roi'])
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
    subVolumeReference = imageReference[
        boundingBox[0] : boundingBox[1], boundingBox[2] : boundingBox[3]
    ]
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
    ) = align2ImagesCrossCorrelation(subVolumeReference, 
                                     subVolume,
                                     lower_threshold=0.1, 
                                     higher_threshold=0.9999999,
                                     upsample_factor=10)
    
    # print("Shift = {}".format(shift))
    # subVolumeCorrected = shiftImage(image2_adjusted, shift)
    
    return shift, error, diffphase, image1_adjusted, image2_adjusted, image2_corrected

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
        pad_top = diff_vert//2
        pad_bottom = diff_vert - pad_top
        diff_hori = width_max - w
        pad_left = diff_hori//2
        pad_right = diff_hori - pad_left
        img_padded = cv2.copyMakeBorder(img, pad_top, pad_bottom, pad_left, pad_right, cv2.BORDER_CONSTANT, value=0)
        assert img_padded.shape[:2] == (height_max, width_max)
        images_padded.append(img_padded)

    return images_padded

def localDriftCorrection_plotsLocalAlignments(imageListCorrected,imageListunCorrected,imageListReference,log1,dataFolder,ROI,barcode):
    '''
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

    '''
    outputFileName = dataFolder.outputFolders["alignImages"] + os.sep + "localDriftCorrection_" + "ROI:"+ROI +"barcode:"+barcode

    imageListCorrectedPadded=pad_images_to_same_size(imageListCorrected)
    imageListunCorrectedPadded=pad_images_to_same_size(imageListunCorrected)
    imageListReferencePadded=pad_images_to_same_size(imageListReference)            

    montage2DReference=montage(imageListReferencePadded)
    montage2DCorrected=montage(imageListCorrectedPadded)
    montage2DunCorrected=montage(imageListunCorrectedPadded)
    
    fig=plt.figure()
    fig.set_size_inches((30, 30))
    plt.imshow(montage2DReference,cmap = 'Blues', alpha = 0.5)
    plt.imshow(montage2DunCorrected,cmap = 'Reds', alpha = 0.5)
    plt.axis("off")
    plt.savefig(outputFileName + "_unCorrected.png")
    plt.close()
    writeString2File(
                log1.fileNameMD,
                "Corrected local drift for ROI: {} barcode:{} \n ![]({})\n".format(ROI, barcode,outputFileName + "_unCorrected.png"),
                "a",
            )

    fig=plt.figure()
    fig.set_size_inches((30, 30))
    plt.imshow(montage2DReference,cmap = 'Blues', alpha = 0.5)
    plt.imshow(montage2DCorrected,cmap = 'Reds', alpha = 0.5)
    plt.axis("off")
    plt.savefig(outputFileName + "_Corrected.png")
    plt.close()
    writeString2File(
                log1.fileNameMD,
                "Corrected local drift for ROI: {} barcode:{} \n ![]({})\n".format(ROI, barcode,outputFileName + "_Corrected.png"),
                "a",
            )
   
    del montage2DReference,montage2DCorrected,imageListReferencePadded,imageListCorrectedPadded
    del imageListReference,imageListCorrected


def localDriftCorrection_savesResults(dictShift,alignmentResultsTable,dataFolder,log1):


    # plots shift violins 
    for ROI in dictShift.keys():
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
        shifts,labels=[],[]
        
        for barcode in dictShift[ROI].keys():
            shiftsBarcode=np.array([0,0])
            labels.append(barcode)
            for iMask in dictShift[ROI][barcode].keys():
                shiftsBarcode=np.concatenate((shiftsBarcode,dictShift[ROI][barcode][iMask]))
            shifts.append(shiftsBarcode)
        
        plt.violinplot(shifts)

        ax1.get_xaxis().set_tick_params(direction='out')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.set_xticks(np.arange(1, len(labels) + 1))
        ax1.set_xticklabels(labels)
        ax1.set_xlim(0.25, len(labels) + 0.75)
        ax1.set_xlabel("ROI #{}".format(ROI))
        ax1.set_ylabel('local drift shift, px')
        outputFileName = dataFolder.outputFolders["alignImages"] + os.sep + "LocalShiftsViolinPlot_" + "ROI:"+ROI 
        plt.savefig(outputFileName + ".png")
        plt.close()                
        writeString2File(
                        log1.fileNameMD,
                        "Local drift for ROI: {}\n ![]({})\n".format(ROI, outputFileName + ".png"),
                        "a",
                    )
                        
    # saves Table with all shifts
    alignmentResultsTable .write(
        dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat", format="ascii.ecsv", overwrite=True,
    )
            
def localDriftCorrection(param, log1, session1):
    sessionName = "localDriftCorrection"
    # tqdm._instances.clear()
    
    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText(
        "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
    )
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(session1.name, param.param["acquisition"]["label"]), "a",
    )


    if 'localShiftTolerance' in param.param['alignImages'].keys():
        shiftTolerance = param.param['alignImages']['localShiftTolerance']
    else:
        shiftTolerance=1

    if 'bezel' in param.param['alignImages'].keys():
        bezel = param.param['alignImages']['bezel']
    else:
        bezel = 20
        
    print("Parameters> shiftTolerance = {} | bezel = {}".format(shiftTolerance,bezel))
    errormessage=[]
    
    alignmentResultsTable = Table(
            names=("aligned file", "reference file", "ROI #", "Barcode #", "CellID #","shift_x", "shift_y", "error", "diffphase",),
            dtype=("S2", "S2", 'int', 'int', 'int',"f4", "f4", "f4", "f4"),
        )
    dictShift={}
    
    for currentFolder in dataFolder.listFolders:
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        dictShift={}
        # iterates over ROIs
        for fileNameDAPI in param.fileList2Process:

            label = param.param["acquisition"]["label"]

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
                raise FileNotFoundError ("I cannot find DAPI mask: {}".format(fullFileNameROImasks))

            # - loads reference fiducial file
            ROI, imReference = loadsFiducial(param, fileNameDAPI,log1,dataFolder)
            imageReference = imReference.removesBackground2D(normalize=True)
            sumImage = np.copy(imageReference)

            print(">>>Image <i> {} | sumImage: {}".format(imReference.data_2D.max(), imageReference.max()))

            # - retrieves list of barcodes for which a fiducial is available in this ROI
            barcodeList, fiducialFileNames = retrieveBarcodeList(param, fileNameDAPI)
            if imReference.fileName in fiducialFileNames:
                fiducialFileNames.remove(imReference.fileName)

            dictShift[ROI]={}
            
            # - load fiducial for cycle <i>
            for barcode, fileNameFiducial in zip(barcodeList, fiducialFileNames):

                # loads 2D image and applies registration
                Im = Image()
                Im.loadImage2D(fileNameFiducial, log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered")
                imageBarcode = Im.removesBackground2D(normalize=True)
                sumImage += imageBarcode

                # - iterate over masks <iMask>
                dictShift[ROI][barcode]={}
                imageListCorrected,imageListunCorrected,imageListReference=[],[],[]
                # for iMask in tqdm(range(1, Masks.max()),mininterval=2,desc="Looping over {} masks for barcode: {}\n".format(Masks.max(),barcode)):
                print("\nprocessing> Looping over {} masks for barcode: {}\n".format(Masks.max(),barcode))
                sleep(0.2)
                for iMask in trange(1, Masks.max()):
                    # calculates shift                    

                    shift, error, diffphase, subVolumeReference, subVolume, subVolumeCorrected =  alignsSubVolumes(imageReference, 
                                                                                                  imageBarcode, 
                                                                                                  Masks, 
                                                                                                  bezel=bezel, 
                                                                                                  iMask=iMask)
                    # stores images in list
                    imageListReference.append(subVolumeReference)
                    imageListunCorrected.append(subVolume)
                    
                    # evaluates shift to determine if we keep or not
                    if  np.nonzero(np.absolute(shift)>shiftTolerance)[0].shape[0]>0:
                        errormessage.append("ROI:{} | barcode:{}| Mask:{}> local shift = {} not kept as it is over the tolerance of {}".format(
                            ROI,
                            barcode,
                            iMask,
                            shift,
                            shiftTolerance))
                        
                        shift=np.array([0,0])
                        imageListCorrected.append(subVolume)
                    else:
                        imageListCorrected.append(subVolumeCorrected)
                    
                    # stores result in database
                    dictShift[ROI][barcode][str(iMask)]=shift

                    # creates Table entry to return
                    tableEntry = [
                        os.path.basename(fileNameFiducial),
                        os.path.basename(imReference.fileName),
                        int(ROI),
                        int(barcode.split('RT')[1]),
                        iMask,
                        shift[0],
                        shift[1],
                        error,
                        diffphase,
                    ]
                    alignmentResultsTable.add_row(tableEntry)

                # output mosaics with global and local alignments
                localDriftCorrection_plotsLocalAlignments(imageListCorrected,imageListunCorrected,imageListReference,log1,dataFolder,ROI,barcode)

            del sumImage

            # produces shift violin plots and saves results Table 
            localDriftCorrection_savesResults(dictShift,alignmentResultsTable,dataFolder,log1)
            log1.report("\n".join(errormessage))

            
    return 0, dictShift, alignmentResultsTable


# # =============================================================================
# # MAIN
# # =============================================================================

# if __name__ == "__main__":

#     parser = argparse.ArgumentParser()
#     parser.add_argument("-F", "--rootFolder", help="Folder with images")
#     args = parser.parse_args()

#     print("\n--------------------------------------------------------------------------")

#     if args.rootFolder:
#         rootFolder = args.rootFolder
#     else:
#         # rootFolder = "/home/marcnol/data/Experiment_20/Embryo_1"
#         # rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
#         rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug"
#         # rootFolder='/home/marcnol/data/Embryo_debug_dataset/rawImages'

#     print("parameters> rootFolder: {}".format(rootFolder))
#     sessionName = "localDriftCorrection"

#     labels2Process = [
#         {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
#         {"label": "barcode", "parameterFile": "infoList_barcode.json"},
#         {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
#     ]

#     # session
#     session1 = session(rootFolder, sessionName)

#     # setup logs
#     log1 = log(rootFolder)
#     # labels2Process indeces: 0 fiducial, 1:
#     labelParameterFile = labels2Process[2]["parameterFile"]
#     param = Parameters(rootFolder, labelParameterFile)

#     dataFolder = folders(param.param["rootFolder"])

#     for currentFolder in dataFolder.listFolders:
#         filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
#         dataFolder.createsFolders(currentFolder, param)

#         # generates lists of files to process
#         param.files2Process(filesFolder)

#         for fileName in param.fileList2Process:
#             session1.add(fileName, sessionName)

#     errorCode, dictShift, alignmentResultsTable = localDriftCorrection(param, log1, session1)

#     if errorCode != 0:
#         print("Error code reported: {}".format(errorCode))
#     else:
#         print("normal termination")

#     # for fileName in param.fileList2Process:
#     #     session1.add(fileName, sessionName)
