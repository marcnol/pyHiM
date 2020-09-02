#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:50:49 2020

@author: marcnol

Purpose: Performs 3D Gaussian fits for barcodes

steps: 
    - iterate over ROIs
    - load mask for ROI <i>
    - check if mask available and load mask file 
    - iterate over cycles <i>
    - load 3D file for barcode <i>
    - iterate over masks <j>
    - 3D gaussian fit of barcodes <k> in mask <j>
    - determine if we keep or not
    - store in database.
    

"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob, os
import matplotlib.pylab as plt
import numpy as np
import uuid
import argparse
from datetime import datetime
from scipy.ndimage import shift as shiftImage
from scipy.optimize import curve_fit
from tqdm import trange

from astropy.stats import sigma_clipped_stats, SigmaClip, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, vstack, Column
from photutils import DAOStarFinder, CircularAperture, detect_sources
from photutils import detect_threshold, deblend_sources
from photutils import Background2D, MedianBackground
from photutils.segmentation.core import SegmentationImage

from imageProcessing.imageProcessing import Image, saveImage2Dcmd
from fileProcessing.fileManagement import (
    folders, writeString2File,loadJSON)
from imageProcessing.projectsBarcodes import projectsBarcodes
from fileProcessing.fileManagement import (session, log, Parameters)

# =============================================================================
# FUNCTIONS
# =============================================================================

def loadsBarcodeMap(dataFolder):
    fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"

    if os.path.exists(fileNameBarcodeCoordinates):
        barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
    else:
        print("\n\n *** ERROR: could not found coordinates file: {}".format(fileNameBarcodeCoordinates))
        return Table(), -1

    # Adds columns to Table to hold the results
    BarcodeMapLength=len(barcodeMap)
    colzPositionGaussian = Column(np.zeros((BarcodeMapLength,1)), name="zcentroidGauss", dtype=float)
    colzPositionMoment = Column(np.zeros((BarcodeMapLength,1)), name="zcentroidMoment", dtype=float)
    colSigmaGaussianFit = Column(np.zeros((BarcodeMapLength,1)), name="sigmaGaussFit", dtype=float)
    colResidualsGaussianFit = Column(np.zeros((BarcodeMapLength,1)), name="residualGaussFit", dtype=float)
    colFitKeep = Column(np.zeros((BarcodeMapLength,1)), name="3DfitKeep", dtype=int)

    barcodeMap.add_column(colzPositionGaussian, index=7)
    barcodeMap.add_column(colzPositionMoment, index=8)
    barcodeMap.add_column(colSigmaGaussianFit, index=17)
    barcodeMap.add_column(colResidualsGaussianFit, index=18)
    barcodeMap.add_column(colFitKeep, index=19)

    return barcodeMap, 0

def loadsShifts3Dimage(param, dataFolder, barcodeMapSinglebarcode):
    # nBarcode = barcodeMapROI_barcodeID.groups.keys[iBarcode][0]  # need to iterate over the first index
    nBarcode = np.unique(barcodeMapSinglebarcode['Barcode #'].data)[0]
    nROI= np.unique(barcodeMapSinglebarcode['ROI #'].data)[0]

    print("Loading 3D image for ROI # {}, barcode # {}".format(nROI,nBarcode))
    imageFile = dataFolder.masterFolder+os.sep+'scan_001_RT41_001_ROI_converted_decon_ch01.tif'
    Im3D=Image()
    Im3D.loadImage(imageFile)
    
    # apply drift correction to 3D image
    dictShifts = loadJSON(dataFolder.outputFiles["dictShifts"] + ".json")

    ROI = param.decodesFileParts(os.path.basename(imageFile))['roi']
    label = os.path.basename(imageFile).split("_")[2]
    try:
        shiftArray = dictShifts["ROI:" + ROI][label]
    except KeyError:
        shiftArray = None
        log1.report(
            "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(ROI, label), "ERROR",
        )
        
    Im3DShifted=Image()
    shift = np.asarray(shiftArray)
    imageShape=Im3D.data.shape
    numberZplanes=imageShape[0]
    Im3DShifted.data=np.zeros(imageShape) # creates array that will hold new shifted 3D data
    print("Shifting 3D image for ROI # {}, barcode # {}".format(nROI,nBarcode))
    for z in range(numberZplanes):
        Im3DShifted.data[z,:,:] = shiftImage(Im3D.data[z,:,:], shift)
        
    # makes 2D projection
    Im3DShifted.maxIntensityProjection()
    Im3DShifted.imageShow(show=True)
           
    return Im3DShifted


def showsImageNsources(im, xcentroids2D, ycentroids2D):
    # show results
    fig = plt.figure()
    fig.set_size_inches((30, 30))

    # positions = np.transpose((xcentroids2D + 0.5, ycentroids2D+ 0.5))
    positions = np.transpose((xcentroids2D, ycentroids2D))
    
    apertures = CircularAperture(positions, r=4.0)
    norm = simple_norm(im, "sqrt", percent=99.9)
    # norm = ImageNormalize(stretch=SqrtStretch())
    # plt.imshow(im1_bkg_substracted, clim=(0, 1), cmap="Greys", origin="lower", norm=norm)
    plt.imshow(im, cmap="Greys", origin="lower", norm=norm)
    apertures.plot(color="blue", lw=1.5, alpha=0.35)
    plt.xlim(0, im.shape[1] - 1)
    plt.ylim(0, im.shape[0] - 1)
    # plt.savefig(outputFileName + "_segmentedSources.png")
    plt.axis("off")
    # plt.close()
    # writeString2File(
    #     log1.fileNameMD,
    #     "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_segmentedSources.png"),
    #     "a",
    # )
    
def getFOV(x,y,window,imageShape):            
    
    fov = {'xleft': int(np.max([1,x-window])),
           'xright': int(np.min([imageShape[1],x+window])),
           'yleft':  int(np.max([1,y-window])),
           'yright': int(np.min([imageShape[2],y+window]))}
    return fov


def getSubVolume(imageData,fov):
    subVolume=imageData[:,fov['yleft']:fov['yright'],fov['xleft']:fov['xright']]
    return subVolume

def subVolume2Trace(subVolumeZscan):
    numberZplanes = subVolumeZscan.shape[0] 
    zTrace=np.zeros(numberZplanes)
    for z in range(numberZplanes):
        zTrace[z] = np.sum(subVolumeZscan[z,:,:],axis=0)

def sum2Dmatrix(matrix):
    sum=np.zeros(1)
    for i in range(matrix.shape[0]):
        sum+=np.sum(matrix[i,:])
    return sum

def getzTrace(subVolume):
    numberZplanes = subVolume.shape[0]
    zTrace=np.zeros((numberZplanes))
    for z in range(numberZplanes):
        zTrace[z]=sum2Dmatrix(subVolume[z,:,:])
    return zTrace
            
def weightedSum1D(zTrace):
    background = np.min(zTrace)
    zTrace -= background
    dims=zTrace.shape[0]
    weightedSum=np.zeros(1)
    
    for z in range(dims):
        weightedSum+=z*zTrace[z]
    weightedSum=weightedSum/np.sum(zTrace)
    return weightedSum

def gaussian(x, a, b, c):
    return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))

def getzPosition(Im3DShifted,x,y,imageShape,window=10):
        
    # retrieve subVolume
    fovSubVolume = getFOV(x,y,window,imageShape)
    subVolume = getSubVolume(Im3DShifted.data,fovSubVolume)
    
    # construct 1D z-scan
    zTrace = getzTrace(subVolume)
            
    # find z-position using center of gravity
    zPositionMoment = weightedSum1D(zTrace)
  
    # Gaussian fitting
    xdata=np.arange(0,imageShape[0])
    amplitude=np.max(zTrace)-np.min(zTrace)
    sigma=4
    pars, cov = curve_fit(f=gaussian, xdata=xdata, ydata=zTrace, p0=[amplitude, zPositionMoment[0], sigma], bounds=(0, [2*amplitude,imageShape[0],imageShape[0]]))
    zPositionGaussian = pars[1]
    sigma = pars[2]
    residual = np.linalg.norm(zTrace-gaussian(xdata,pars[0],pars[1],pars[2])) / pars[0]

    fitKeep = True
    
    # # plotting    
    # plt.plot(zTrace)
    # plt.axvline(x=zMax,ymin=0,ymax=np.max(zTrace), color='g')
    # plt.plot(xdata, gaussian(xdata,pars[0],pars[1],pars[2]), '-r', label='gaussian',alpha=0.5)    
    # plt.axvline(x=zPositionGaussian,ymin=0,ymax=np.max(zTrace), color='r')
    
    return zPositionMoment, zPositionGaussian, sigma, residual, fitKeep

def fitsZpositions(Im3DShifted,barcodeMapSinglebarcode,window):            
    imageShape=Im3DShifted.data.shape
    xcentroids2D=barcodeMapSinglebarcode['xcentroid'].data
    ycentroids2D=barcodeMapSinglebarcode['ycentroid'].data
    numberSpots = len(xcentroids2D)
    
    listZpositionsMoment,listZpositionsGaussian, listSigma, listResiduals, listFitKeep = [], [] , [], [], []
    print("Looping over {} spots...".format(numberSpots))
    for iSpot in trange(numberSpots):
        zPositionMoment, zPositionGaussian,sigma, res, fitKeep = getzPosition(Im3DShifted,xcentroids2D[iSpot],ycentroids2D[iSpot],imageShape,window)
        listZpositionsMoment.append(zPositionMoment)
        listZpositionsGaussian.append(zPositionGaussian)
        listSigma.append(sigma)
        listResiduals.append(res)
        listFitKeep.append(fitKeep)
        
    barcodeMapSinglebarcode['zcentroidGauss']=listZpositionsGaussian
    barcodeMapSinglebarcode['zcentroidMoment']=listZpositionsMoment
    barcodeMapSinglebarcode['sigmaGaussFit']=listSigma
    barcodeMapSinglebarcode['residualGaussFit']=listResiduals

    return barcodeMapSinglebarcode

def shows3DfittingResults(barcodeMapSinglebarcode,numberZplanes=60):
    nBarcode = np.unique(barcodeMapSinglebarcode['Barcode #'].data)[0]
    ROI= np.unique(barcodeMapSinglebarcode['ROI #'].data)[0]

    plt.subplot(211)
    plt.title('ROI: {} | barcode: {}'.format(ROI,nBarcode))
    cs1=plt.scatter(barcodeMapSinglebarcode['residualGaussFit'],barcodeMapSinglebarcode['sigmaGaussFit'], c=barcodeMapSinglebarcode['zcentroidGauss'], alpha=0.3)
    plt.xlabel('sigma')
    plt.ylabel('residuals/amplitude')
    cbar1=plt.colorbar(cs1)
    cbar1.set_label('zPosition')

    plt.subplot(212)
    cs2=plt.scatter(barcodeMapSinglebarcode['zcentroidGauss'],barcodeMapSinglebarcode['zcentroidMoment'], s=(barcodeMapSinglebarcode['sigmaGaussFit']), c=barcodeMapSinglebarcode['residualGaussFit'], alpha=0.5)
    plt.xlabel('zPosition, Gaussian')
    plt.ylabel('zPosition, Moment')
    plt.xlim(0,numberZplanes)
    plt.ylim(0,numberZplanes)
    cbar2=plt.colorbar(cs2)
    cbar2.set_label('residuals/amplitude')
    plt.clim(0,10) 


#%%
def makes3Dfits(param, dataFolder, log1):
    
    # Loads coordinate Tables for all barcodes and ROIs
    barcodeMap, errorCode = loadsBarcodeMap(dataFolder)

    window=3

    # loops over ROIs
    barcodeMapROI = barcodeMap.group_by("ROI #")
    SCmatrixCollated, uniqueBarcodes = [], []
    numberROIs = len(barcodeMapROI.groups.keys)
    print("\nROIs detected: {}".format(barcodeMapROI.groups.keys))

    for iROI in range(numberROIs):

        nROI = barcodeMapROI.groups.keys[iROI][0]  # need to iterate over the first index
        print("Working on ROI# {}".format(nROI))

        # loops over barcodes in that ROI
        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[iROI]
        barcodeMapROI_barcodeID = barcodeMapSingleROI.group_by("Barcode #")
        numberBarcodes=len(barcodeMapROI_barcodeID.groups.keys)
        print("\nNumber of barcodes detected: {}".format(numberBarcodes))

        for iBarcode in range(numberBarcodes):

            # find coordinates for this ROI and barcode
            barcodeMapSinglebarcode = barcodeMapROI_barcodeID.group_by("Barcode #").groups[iBarcode]

            # load 3D image: NEED TO MAKE IT RETRIEVE IMAGE NAME FROM BARCODE AND ROI !
            Im3DShifted = loadsShifts3Dimage(param, dataFolder, barcodeMapSinglebarcode)
            numberZplanes = Im3DShifted.data.shape[0]
            
            # shows 2D images and detected sources
            # showsImageNsources(Im3DShifted.data_2D, xcentroids2D, ycentroids2D)
           
            # loop over spots
            barcodeMapSinglebarcode = fitsZpositions(Im3DShifted,barcodeMapSinglebarcode,window)

            # displays results
            shows3DfittingResults(barcodeMapSinglebarcode,numberZplanes=numberZplanes)

            # record results by appending the ASTROPY table *** use index first then match BUIDs in barcodeMapSinglebarcode to
            # those in barcodeMapROI and replace values of the row in barcodeMapROI by those in barcodeMapSinglebarcode
            
            


#%%

def refitBarcodes3D(param, log1, session1):
    sessionName="refitBarcodes3D"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD, "## {}\n".format(sessionName), "a")
    
    for currentFolder in dataFolder.listFolders:
        # filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder, param)
        log1.report("-------> Processing Folder: {}".format(currentFolder))

        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]

        makes3Dfits(param,dataFolder, log1)
        session1.add(currentFolder, sessionName)

        log1.report("HiM matrix in {} processed".format(currentFolder), "info")
    

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--rootFolder", help="Folder with images")
    args = parser.parse_args()
    now = datetime.now()

    print("\n--------------------------------------------------------------------------")

    if args.rootFolder:
        rootFolder = args.rootFolder
    else:
        rootFolder ='/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0'

    print("parameters> rootFolder: {}".format(rootFolder))

    labels2Process = [
        {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
        {"label": "barcode", "parameterFile": "infoList_barcode.json"},
        {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
        {"label": "RNA", "parameterFile": "infoList_RNA.json"},
    ]

    # session
    sessionName="refitBarcodes3D"
    session1 = session(rootFolder, sessionName)

    # setup logs
    log1 = log(rootFolder)
    log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format("processingPipeline"))
    log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
    writeString2File(
        log1.fileNameMD, "# Hi-M {}: {}".format(sessionName,now.strftime("%Y/%m/%d %H:%M:%S")), "w",
    )  # initialises MD file

    for ilabel in range(len(labels2Process)):
        label = labels2Process[ilabel]["label"]
        labelParameterFile = labels2Process[ilabel]["parameterFile"]
        log1.addSimpleText("**Analyzing label: {}**".format(label))

        # sets parameters
        param = Parameters(rootFolder, labelParameterFile)

        # [builds PWD matrix for all folders with images]
        if label == "barcode":
            refitBarcodes3D(param, log1, session1)

    

    
    
    
    
    
    
    












   
 