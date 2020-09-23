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
from datetime import datetime
from scipy.ndimage import shift as shiftImage
from scipy.optimize import curve_fit
from shutil import copyfile

from tqdm import trange, tqdm
from astropy.visualization import simple_norm
from astropy.table import Table, Column
from photutils import CircularAperture
from dask.distributed import Client, LocalCluster, get_client, as_completed

from numba import jit

from imageProcessing.imageProcessing import Image
from fileProcessing.fileManagement import folders, writeString2File, loadJSON
from fileProcessing.fileManagement import daskCluster

# =============================================================================
# FUNCTIONS
# =============================================================================


class refitBarcodesClass:
    def __init__(self, param, log1, session1, parallel=False):
        self.param = param
        self.session1=session1
        # self.dataFolder = []
        self.log1 = log1
        self.window = 3
        self.parallel=parallel

    def loadsBarcodeMap(self):
        fileNameBarcodeCoordinates = self.dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"

        if os.path.exists(fileNameBarcodeCoordinates):
            barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
        else:
            print("\n\n *** ERROR: could not found coordinates file: {}".format(fileNameBarcodeCoordinates))
            return Table(), -1

        # Adds columns to Table to hold the results
        BarcodeMapLength = len(barcodeMap)
        if "zcentroidGauss" not in barcodeMap.keys():
            colzPositionGaussian = Column(np.zeros((BarcodeMapLength)), name="zcentroidGauss", dtype=float)
            barcodeMap.add_column(colzPositionGaussian, index=7)
        if "zcentroidMoment" not in barcodeMap.keys():
            colzPositionMoment = Column(np.zeros((BarcodeMapLength)), name="zcentroidMoment", dtype=float)
            barcodeMap.add_column(colzPositionMoment, index=8)
        if "sigmaGaussFit" not in barcodeMap.keys():
            colSigmaGaussianFit = Column(np.zeros((BarcodeMapLength)), name="sigmaGaussFit", dtype=float)
            barcodeMap.add_column(colSigmaGaussianFit, index=17)
        if "residualGaussFit" not in barcodeMap.keys():
            colResidualsGaussianFit = Column(np.zeros((BarcodeMapLength)), name="residualGaussFit", dtype=float)
            barcodeMap.add_column(colResidualsGaussianFit, index=18)
        if "3DfitKeep" not in barcodeMap.keys():
            colFitKeep = Column(np.zeros((BarcodeMapLength)), name="3DfitKeep", dtype=int)
            barcodeMap.add_column(colFitKeep, index=19)

        return barcodeMap, 0

    def findsFile2Process(self, nBarcode, nROI):
        Barcode = "RT" + str(nBarcode)
        ROI = str(nROI) + "_ROI"
        channelbarcode = self.param.setsChannel("barcode_channel", "ch01")

        filesFolder = glob.glob(self.dataFolder.masterFolder + os.sep + "*.tif")
        imageFile = [x for x in filesFolder if ROI in x and Barcode in x and channelbarcode in x]

        return imageFile

    def loadsShifts3Dimage(self, barcodeMapSinglebarcode):
        nBarcode = np.unique(barcodeMapSinglebarcode["Barcode #"].data)[0]
        nROI = np.unique(barcodeMapSinglebarcode["ROI #"].data)[0]

        imageFile = self.findsFile2Process(nBarcode, nROI)

        if len(imageFile) > 0:
            self.log1.report("Loading 3D image for ROI # {}, barcode # {}".format(nROI, nBarcode))

            Im3D = Image(self.param,self.log1)
            Im3D.loadImage(imageFile[0])

            # corrects drift for all barcodes, except the fiducial
            if "RT" + str(nBarcode) not in self.param.param["alignImages"]["referenceFiducial"]:

                # apply drift correction to 3D image
                dictShifts = loadJSON(self.dataFolder.outputFiles["dictShifts"] + ".json")

                ROI = self.param.decodesFileParts(os.path.basename(imageFile[0]))["roi"]
                label = os.path.basename(imageFile[0]).split("_")[2]
                try:
                    shiftArray = dictShifts["ROI:" + ROI][label]
                except KeyError:
                    shiftArray = None
                    self.log1.report(
                        "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(ROI, label), "ERROR",
                    )

                Im3DShifted = Image(self.param,self.log1)
                shift = np.asarray(shiftArray)
                imageShape = Im3D.data.shape
                numberZplanes = imageShape[0]
                Im3DShifted.data = np.zeros(imageShape)  # creates array that will hold new shifted 3D data
                self.log1.report("Shifting 3D image for ROI # {}, barcode # {}".format(nROI, nBarcode))
                
                if self.parallel:
                    R = range(numberZplanes)
                else:
                    R = trange(numberZplanes)
                    
                for z in R:
                    Im3DShifted.data[z, :, :] = shiftImage(Im3D.data[z, :, :], shift)

            else:
                # this is run for the fiducial
                Im3DShifted = Im3D

            # makes 2D projection
            Im3DShifted.maxIntensityProjection()
            # Im3DShifted.imageShow(show=True)

            return Im3DShifted
        else:
            return Image(self.param,self.log1)

    def showsImageNsources(self, im, xcentroids2D, ycentroids2D):
        # show results
        fig = plt.figure()
        fig.set_size_inches((30, 30))

        positions = np.transpose((xcentroids2D, ycentroids2D))

        apertures = CircularAperture(positions, r=4.0)
        norm = simple_norm(im, "sqrt", percent=99.9)
        plt.imshow(im, cmap="Greys", origin="lower", norm=norm)
        apertures.plot(color="blue", lw=1.5, alpha=0.35)
        plt.xlim(0, im.shape[1] - 1)
        plt.ylim(0, im.shape[0] - 1)
        plt.axis("off")

    def getFOV(self, x, y, imageShape):

        fov = {
            "xleft": int(np.max([1, x - self.window])),
            "xright": int(np.min([imageShape[1], x + self.window])),
            "yleft": int(np.max([1, y - self.window])),
            "yright": int(np.min([imageShape[2], y + self.window])),
        }
        return fov

    def getSubVolume(self, imageData, fov):
        subVolume = imageData[:, fov["yleft"] : fov["yright"], fov["xleft"] : fov["xright"]]
        return subVolume

    def subVolume2Trace(self, subVolumeZscan):
        numberZplanes = subVolumeZscan.shape[0]
        zTrace = np.zeros(numberZplanes)
        for z in range(numberZplanes):
            zTrace[z] = np.sum(subVolumeZscan[z, :, :], axis=0)

    def sum2Dmatrix(self, matrix):
        sum = np.zeros(1)
        for i in range(matrix.shape[0]):
            sum += np.sum(matrix[i, :])
        return sum

    def getzTrace(self, subVolume):
        numberZplanes = subVolume.shape[0]
        zTrace = np.zeros((numberZplanes))
        for z in range(numberZplanes):
            zTrace[z] = self.sum2Dmatrix(subVolume[z, :, :])
        return zTrace

    def weightedSum1D(self, zTrace):
        background = np.min(zTrace)
        zTrace -= background
        dims = zTrace.shape[0]
        weightedSum = np.zeros(1)

        for z in range(dims):
            weightedSum += z * zTrace[z]
        weightedSum = weightedSum / np.sum(zTrace)
        return weightedSum

    def gaussian(self, x, a, b, c):
        return a * np.exp(-np.power(x - b, 2) / (2 * np.power(c, 2)))

    # @jit(nopython=True)
    def getzPosition(self, Im3DShifted, x, y, imageShape):
        """
        Finds z position of PSF using moment and gaussian fitting    
        Given an image and (xy) coordinates, it extracts a subvolume of size 2*window+1
        calculates the trace of intensity in z,
        calculates moment
        uses as a seed to perform gaussian fitting
    
        Parameters
        ----------
        Im3DShifted : Image Class
            3D image
        x : float
            x coordinate.
        y : float
            y coordinate.
        imageShape : np.array
            size of the image in Im3DShifted.data
        window : int
            Size of the window for taking a subVolume. The default is 10.
    
        Returns
        -------
        zPositionMoment: float, position of z centroid based on moment
        zPositionGaussian: float, position of z centroid based on gaussian fit
        sigma: width of the gaussian
        residual: residual of the fit
        fitKeep: flag indicating if fit is to be used or not
        fitSuccess: indicates whether there were problems calling curve_fit
        """

        # retrieve subVolume
        fovSubVolume = self.getFOV(x, y, imageShape)
        subVolume = self.getSubVolume(Im3DShifted.data, fovSubVolume)

        # construct 1D z-scan
        zTrace = self.getzTrace(subVolume)

        # find z-position using center of gravity
        zPositionMoment = self.weightedSum1D(zTrace)

        # Gaussian fitting
        xdata = np.arange(0, imageShape[0])
        amplitude = np.max(zTrace) - np.min(zTrace)

        width = 4

        try:
            pars, cov = curve_fit(
                f=self.gaussian,
                xdata=xdata,
                ydata=zTrace,
                p0=[amplitude, zPositionMoment[0], width],
                bounds=(0, [2 * amplitude, imageShape[0], imageShape[0]]),
            )
            fitSuccess = True
        except ValueError:
            self.log1.report("Gaussian fitting: out of bounds")
            fitSuccess = False
        except RuntimeError:
            self.log1.report("Gaussian fitting: did not converge")
            fitSuccess = False

        if "residual_max" in self.param.param.keys():
            residualThreshold = self.param.param["residual_max"]
        else:
            residualThreshold = 2.5
        if "sigma_max" in self.param.param.keys():        
            sigmaThreshold = self.param.param["sigma_max"]
        else:
            sigmaThreshold = 5 
        if "centroidDifference_max" in self.param.param.keys():        
            centroidMaxDifference = self.param.param["centroidDifference_max"] 
        else:
            centroidMaxDifference = 5
        
        if fitSuccess:
            zPositionGaussian = pars[1]
            sigma = pars[2]
            residual = np.linalg.norm(zTrace - self.gaussian(xdata, pars[0], pars[1], pars[2])) / pars[0]
            
            if residual < residualThreshold\
                and sigma < sigmaThreshold and\
                    (zPositionGaussian - zPositionMoment) < centroidMaxDifference : 
                fitKeep = True
            else:
                fitKeep = False                
        else:
            zPositionGaussian, sigma, residual = zPositionMoment, np.nan, np.nan
            fitKeep = False

        # # plotting
        # plt.plot(zTrace)
        # plt.axvline(x=zMax,ymin=0,ymax=np.max(zTrace), color='g')
        # plt.plot(xdata, gaussian(xdata,pars[0],pars[1],pars[2]), '-r', label='gaussian',alpha=0.5)
        # plt.axvline(x=zPositionGaussian,ymin=0,ymax=np.max(zTrace), color='r')

        return zPositionMoment, zPositionGaussian, sigma, residual, fitKeep, fitSuccess

    def fitsZpositions(self, Im3DShifted, barcodeMapSinglebarcode):
        imageShape = Im3DShifted.data.shape
        xcentroids2D = barcodeMapSinglebarcode["xcentroid"].data
        ycentroids2D = barcodeMapSinglebarcode["ycentroid"].data
        numberSpots = len(xcentroids2D)

        listZpositionsMoment, listZpositionsGaussian, listSigma, listResiduals, listFitKeep, listfitSuccess = (
            [],
            [],
            [],
            [],
            [],
            [],
        )
        
        self.log1.report("Looping over {} spots...".format(numberSpots))
 
        if self.parallel:
            R = range(numberSpots)      
        else:
            R = trange(numberSpots)

        for iSpot in R:
            results = self.getzPosition(Im3DShifted, xcentroids2D[iSpot], ycentroids2D[iSpot], imageShape)
            zPositionMoment, zPositionGaussian, sigma, res, fitKeep, fitSuccess = results
            listZpositionsMoment.append(zPositionMoment)
            listZpositionsGaussian.append(zPositionGaussian)
            listSigma.append(sigma)
            listResiduals.append(res)
            listFitKeep.append(fitKeep)
            listfitSuccess.append(fitSuccess)

        self.log1.report("Unfailed fittings: {} out of {} ".format(sum(listfitSuccess), numberSpots))
        barcodeMapSinglebarcode["zcentroidGauss"] = listZpositionsGaussian
        barcodeMapSinglebarcode["zcentroidMoment"] = listZpositionsMoment
        barcodeMapSinglebarcode["sigmaGaussFit"] = listSigma
        barcodeMapSinglebarcode["residualGaussFit"] = listResiduals
        barcodeMapSinglebarcode["3DfitKeep"] = listFitKeep

        return barcodeMapSinglebarcode

    def shows3DfittingResults(self, barcodeMapSinglebarcode, numberZplanes=60, show=False):
        nBarcode = np.unique(barcodeMapSinglebarcode["Barcode #"].data)[0]
        ROI = np.unique(barcodeMapSinglebarcode["ROI #"].data)[0]
        outputFileName = self.outputFileName + "_3Drefit_ROI:" + str(ROI) + "_barcode:" + str(nBarcode) + ".png"        

        # filters list and keeps only those with 3DfitKeep == True
        filteredCentroidMoment = []
        filteredCentroidGauss = []
        filteredResidual = []
        filteredSigma = []
        if "centroidDifference_max" in self.param.param.keys():        
            fluxThreshold = self.param.param["fluxThreshold"]
        else:
            fluxThreshold = 200
            
        for item in barcodeMapSinglebarcode:
            if item["3DfitKeep"] and item["flux"] > fluxThreshold :
                filteredCentroidMoment.append(item["zcentroidMoment"])    
                filteredCentroidGauss.append(item["zcentroidGauss"])    
                filteredResidual.append(item["residualGaussFit"])    
                filteredSigma.append(item["sigmaGaussFit"])    
            
        # plots first subpanel with centroid vs residuals
        fig, (ax2, ax1) = plt.subplots(2, 1)
        cs1a = ax2.scatter(
            barcodeMapSinglebarcode["zcentroidGauss"],
            barcodeMapSinglebarcode["residualGaussFit"],
            s=barcodeMapSinglebarcode["sigmaGaussFit"],
            c=barcodeMapSinglebarcode["sigmaGaussFit"],
            alpha=0.3,
        )
        cs1b = ax2.scatter(
            filteredCentroidGauss,
            filteredResidual,
            marker=".",c='k',
            alpha=0.5,
        )

        ax1.set_xlabel("zPosition, Gaussian")
        ax2.set_ylabel("residuals/amplitude")
        # cbar1 = ax1.colorbar(cs1)
        cbar1 = fig.colorbar(cs1a, ax=ax2)
        cbar1.set_label("sigmaGaussFit")
        cs1a.set_clim(0, 10)

        # plots second subpanel with zCentroid from moment and from Gaussian fits
        cs2 = ax1.scatter(
            barcodeMapSinglebarcode["zcentroidGauss"],
            barcodeMapSinglebarcode["zcentroidMoment"],
            s=(barcodeMapSinglebarcode["sigmaGaussFit"]),
            c=barcodeMapSinglebarcode["residualGaussFit"],
            alpha=0.5,
        )
        cs1b = ax1.scatter(
            filteredCentroidGauss,
            filteredCentroidMoment,
            marker=".",c='k',
            alpha=0.5,
        )
        ax1.set_xlabel("zPosition, Gaussian")
        ax1.set_ylabel("zPosition, Moment")
        ax1.set_xlim(0, numberZplanes)
        ax1.set_ylim(0, numberZplanes)
        cbar2 = fig.colorbar(cs2,ax=ax1)        
        cbar2.set_label("residuals/amplitude")
        cs2.set_clim(0, 10)

        fig.suptitle("ROI: {} | barcode: {}".format(ROI, nBarcode), fontsize=12)

        # saves output
        plt.savefig(outputFileName)
        
        # closes plot
        if not show:
            plt.close()

        # write line in MD file pointing to plot
        writeString2File(
            self.log1.fileNameMD, "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName), "a",
        )

    def refitsBarcode(self, barcodeMapSinglebarcode):
        '''
        Refits a barcode encoded in barcodeMapSinglebarcode
        
        Parameters
        ----------
        barcodeMapSinglebarcode : ASTROPY table
            List of coordinates detected in an ROI for a specific barcode.

        Returns
        -------
        barcodeMapSinglebarcode : ASTROPY table
            ASTROPY table with the fitted z centroids.

        '''
        # load 3D image
        Im3DShifted = self.loadsShifts3Dimage(barcodeMapSinglebarcode)
        numberZplanes = Im3DShifted.data.shape[0]

        # shows 2D images and detected sources
        # self.showsImageNsources(Im3DShifted.data_2D, xcentroids2D, ycentroids2D)
          
        # loop over spots
        barcodeMapSinglebarcode = self.fitsZpositions(Im3DShifted, barcodeMapSinglebarcode)

        # displays results
        self.shows3DfittingResults(barcodeMapSinglebarcode, numberZplanes=numberZplanes)

        return barcodeMapSinglebarcode
    
    def applyResults(self,barcodeMap, result):
        for iSpot in range(len(result)):
            barcodeMap.loc[result['Buid'][iSpot]]=result[iSpot]

    def rewritesBarcodeMap(self,barcodeMap):            
        fileNameBarcodeCoordinates = self.dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"
        fileNameBarcodeCoordinatesOld = self.dataFolder.outputFiles["segmentedObjects"] + "_barcode2D.dat"

        copyfile(fileNameBarcodeCoordinates, fileNameBarcodeCoordinatesOld)
                
        barcodeMap.write(fileNameBarcodeCoordinates, format="ascii.ecsv", overwrite=True)
    
    def refitFilesinFolder(self):
        '''
        Refits all the barcode files found in rootFolder

        Returns
        -------
        None.

        '''
        now = datetime.now()

        # Loads coordinate Tables for all barcodes and ROIs
        barcodeMap, errorCode = self.loadsBarcodeMap()

        # retrieves parameters
        if "3DGaussianfitWindow" in self.param.param["segmentedObjects"].keys():
            self.window = self.param.param["segmentedObjects"]["3DGaussianfitWindow"]
        else:
            self.window = 3

        # loops over ROIs
        barcodeMapROI = barcodeMap.group_by("ROI #")
        numberROIs = len(barcodeMapROI.groups.keys)
        self.log1.info("\nROIs detected: {}".format(barcodeMapROI.groups.keys))

        # self.param.param['parallel']=False
        availableBarcodes = np.unique(barcodeMap["Barcode #"].data)
        maxnumberBarcodes = availableBarcodes.shape[0]
        self.log1.info("Max number of barcodes detected: {}".format(maxnumberBarcodes))
        
        if self.parallel:
            futures = list()
            
            client=get_client()
            
            self.log1.info("Go to http://localhost:8787/status for information on progress...")
            
            for iROI in range(numberROIs):
                nROI = barcodeMapROI.groups.keys[iROI][0]  # need to iterate over the first index
                self.log1.info("Working on ROI# {}".format(nROI))
    
                # loops over barcodes in that ROI
                barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[iROI]
                barcodeMapROI_barcodeID = barcodeMapSingleROI.group_by("Barcode #")
                numberBarcodes = len(barcodeMapROI_barcodeID.groups.keys)
                self.log1.info("\nNumber of barcodes detected: {}".format(numberBarcodes))
    
                for iBarcode in range(numberBarcodes):
                    # find coordinates for this ROI and barcode
                    barcodeMapSinglebarcode = barcodeMapROI_barcodeID.group_by("Barcode #").groups[iBarcode]
                    result = client.submit(self.refitsBarcode, barcodeMapSinglebarcode)
                    futures.append(result)    
            
            self.log1.info("Waiting for {} results to arrive".format(len(futures)))

            results = client.gather(futures)       

            self.log1.info("{} results retrieved from cluster".format(len(results)))

            del futures
                 
        else:
            results = []
            for iROI in range(numberROIs):
                nROI = barcodeMapROI.groups.keys[iROI][0]  # need to iterate over the first index
                self.log1.report("Working on ROI# {}".format(nROI))
    
                # loops over barcodes in that ROI
                barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[iROI]
                barcodeMapROI_barcodeID = barcodeMapSingleROI.group_by("Barcode #")
                numberBarcodes = len(barcodeMapROI_barcodeID.groups.keys)
                self.log1.report("\nNumber of barcodes detected: {}".format(numberBarcodes))
    
                for iBarcode in range(numberBarcodes):
                    # find coordinates for this ROI and barcode
                    barcodeMapSinglebarcode = barcodeMapROI_barcodeID.group_by("Barcode #").groups[iBarcode]
                    result = self.refitsBarcode(barcodeMapSinglebarcode)
                    results.append(result)

        # record results by appending the ASTROPY table *** use index first then match BUIDs in barcodeMapSinglebarcode to
        # those in barcodeMapROI and replace values of the row in barcodeMapROI by those in barcodeMapSinglebarcode
        self.log1.report("Recording results...")
        barcodeMap.add_index('Buid')

        for result in tqdm(results):
            # iterates over bbuid's
            for iSpot in range(len(result)):
                # matches bbuid's
                barcodeMap.loc[result[iSpot]['Buid']]=result[iSpot]
                
        print("RefitFilesinFolder time: {}".format(datetime.now() - now))

        self.rewritesBarcodeMap(barcodeMap)
        
    def refitFolders(self):
        '''
        runs refitting routine in rootFolder

        Returns
        -------
        None.

        '''
        sessionName = "refitBarcodes3D"
    
        # processes folders and files
        self.log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        self.log1.report("folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.log1.fileNameMD, "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        currentFolder = self.dataFolder.listFolders[0]
        
        self.dataFolder.createsFolders(currentFolder, self.param)
        self.outputFileName = self.dataFolder.outputFiles["segmentedObjects"]

        self.log1.report("-------> Processing Folder: {}".format(currentFolder))
        self.log1.parallel = self.parallel
        
        self.refitFilesinFolder()

        self.session1.add(currentFolder, sessionName)

        self.log1.report("HiM matrix in {} processed".format(currentFolder), "info")

        return 0