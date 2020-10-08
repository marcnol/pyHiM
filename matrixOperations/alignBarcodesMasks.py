#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

These scripts assign barcodes to DAPI masks, calculates the pair-wise distances 
for each barcode set in a cell, and computes the single-cell PWD matrix and 
its ensemble (which is represented).

This file contains as well some tools for representation of matrices.
    

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os, sys
import uuid
import re
import numpy as np
from tqdm.contrib import tzip
from tqdm import trange

from sklearn.metrics import pairwise_distances

from astropy.table import Table

from photutils.segmentation import SegmentationImage

from fileProcessing.fileManagement import (
    folders,
    writeString2File,
    )

from matrixOperations.HIMmatrixOperations import plotMatrix, plotDistanceHistograms

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# CLASSES
# =============================================================================


class cellID:
    def __init__(self, param, barcodeMapROI, Masks, ROI,ndims=2):
        self.param=param
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0
        self.ndims=ndims
        self.dictErrorBlockMasks={} # contains the results from blockAlignment, if existing
        
        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI
        self.alignmentResultsTable=Table()
        self.barcodesinMask = dict()
        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def initializeLists(self):
        self.ROIs, self.cellID, self.nBarcodes, self.barcodeIDs, self.p, self.cuid, self.buid = [], [], [], [], [], [], []

    # def visualize(self):
    #     pass
    #     # imageBarcodes = np.zeros([2048, 2048])
    #     # MasksBarcodes = Masks
    #     # R = []

    #     # for i in range(len(self.barcodeMapROI.groups[0])):
    #     #     y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
    #     #     x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
    #     #     barcodeID = self.barcodeMapROI.groups[0]["Barcode #"][i]
    #     #     imageBarcodes[x_int][y_int] = barcodeID
    #     #     MasksBarcodes[x_int][y_int] += 20 * barcodeID
    #     #     R.append([y_int, x_int, barcodeID])

    #     # # Shows results
    #     # Ra = np.array(R)
    #     # # plt.imshow(Masks, origin="lower", cmap="jet")
    #     # plt.scatter(Ra[:, 0], Ra[:, 1], s=5, c=Ra[:, 2], alpha=0.5)

    def filterLocalizations_Quality(self,i,flux_min):
        """
        [filters barcode localizations either by brigthness or 3D localization accuracy]

        Parameters
        ----------
        i : int
            index in barcodeMap Table
        flux_min : float
            Minimum flux to keep barcode localization

        Returns
        -------
        keep : Boolean
            True if the test is passed.

        """
        if "3DfitKeep" in self.barcodeMapROI.groups[0].keys() and self.ndims==3:
            # [reading the flag in barcodeMapROI assigned by the 3D localization routine]
            keep = self.barcodeMapROI.groups[0]["3DfitKeep"][i]
        else:
            # [or by reading the flux from 2D localization]
            keep = self.barcodeMapROI.groups[0]["flux"][i] > flux_min 

        return keep


    def filterLocalizations_BlockAlignment(self,i,toleranceDrift,blockSize):
        """
        [filters barcode per blockAlignmentMask, if existing]            
        runs only if localAligment was not run!     

        Parameters
        ----------
        i : int
            index in barcodeMap Table
        toleranceDrift : float
            tolerance to keep barcode localization, in pixel units
        blockSize : int
            size of blocks used for blockAlignment.

        Returns
        -------
        keepAlignment : Boolean
            True if the test is passed.

        """
        y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
        x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
        keepAlignment = True
        if not self.alignmentResultsTableRead:
            barcodeID = 'barcode:' + str(self.barcodeMapROI.groups[0]["Barcode #"][i])
            barcodeROI = 'ROI:' + str(self.barcodeMapROI.groups[0]["ROI #"][i])
            
            if len(self.dictErrorBlockMasks)>0:
                if barcodeROI in self.dictErrorBlockMasks.keys():
                    if barcodeID in self.dictErrorBlockMasks[barcodeROI].keys():
                        errorMask = self.dictErrorBlockMasks[barcodeROI][barcodeID]
                        keepAlignment = errorMask[int(np.floor(x_int/blockSize)),int(np.floor(y_int/blockSize))] < toleranceDrift 
            
            # keeps it always if barcode is fiducial
            if  "RT" + str(self.barcodeMapROI.groups[0]["Barcode #"][i]) in self.param.param["alignImages"]["referenceFiducial"]:
                keepAlignment=True
        return keepAlignment 


    def alignByMasking(self):
        '''
        Assigns barcodes to masks and creates <NbarcodesinMask>


        Returns
        -------
        self.barcodesinMask # dictionnary with the identities of barcodes contained in each mask.
            Keys: 'maskID_1', 'maskID_2', and so on
            
        self.NbarcodesinMask # vector containing the number of barcodes for each mask
        self.NcellsAssigned # number of cells assigned
        self.NcellsUnAssigned # number of cells unassigned
        '''

        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        if "flux_min" in self.param.param["segmentedObjects"]:
            flux_min = self.param.param["segmentedObjects"]["flux_min"]
        else:
            flux_min = 0
        if "toleranceDrift" in self.param.param["segmentedObjects"]:
            toleranceDrift = self.param.param["segmentedObjects"]["toleranceDrift"]
        else:
            toleranceDrift= 100

        print("Flux min = {} | ToleranceDrift = {} px".format(flux_min,toleranceDrift))
        
        # Produces images of distribution of fluxes.
        
        
        blockSize=256
        keepAll, keepAlignmentAll, NbarcodesROI=[],[], 0
        # loops over barcode Table rows in a given ROI
        for i in trange(len(self.barcodeMapROI.groups[0])):

            # [filters barcode localizations either by]
            keep = self.filterLocalizations_Quality(i,flux_min)
            
            # [filters barcode per blockAlignmentMask, if existing]            
            keepAlignment = self.filterLocalizations_BlockAlignment(i,toleranceDrift,blockSize)
            
            # applies all filters
            if keep and keepAlignment:
                # keeps the particle if the test passed
                y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
                x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
                
                #finds what mask label this barcode is sitting on
                maskID = self.Masks[x_int][y_int]
    
                # attributes CellID to a barcode            
                self.barcodeMapROI["CellID #"][i] = maskID
                
                # if it is not background,
                if maskID > 0:
                    # increments counter of number of barcodes in the cell mask attributed
                    NbarcodesinMask[maskID] += 1
                    
                    # stores the identify of the barcode to the mask
                    self.barcodesinMask["maskID_" + str(maskID)].append(i)

            # keeps statistics
            if int(self.barcodeMapROI.groups[0]["ROI #"][i]) == int(self.nROI):
                keepAll.append(keep)
                keepAlignmentAll.append(keepAlignment)
                NbarcodesROI+=1

        # Total number of masks assigned and not assigned
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned

        # this list contains which barcodes are allocated to which masks
        self.NbarcodesinMask = NbarcodesinMask
        
        print("KeepQuality: {}| keepAlignment: {}| Total: {}".format(
            sum(keepAll),
            sum(keepAlignmentAll),
            NbarcodesROI))        

        print("Number cells assigned: {} | discarded: {}".format(
            self.NcellsAssigned,
            self.NcellsUnAssigned))        
    
    def searchLocalShift(self,ROI,CellID,x_uncorrected,y_uncorrected):
        '''
        Searches for local drift for current mask. If it exists then id adds it to the uncorrected coordinates

        Parameters
        ----------
        ROI : string
            ROI used
        CellID: string
            ID of the cell
        x_uncorrected : float
            x coordinate.
        y_uncorrected : float
            y coordinate.

        Returns
        -------
        x_corrected : float
            corrected x coordinate.
        y_corrected : float
            corrected y coordinate.

        '''
        _foundMatch=False
        x_corrected, y_corrected = [], []
        for row in self.alignmentResultsTable:
            if row["ROI #"]==ROI and row["CellID #"]==CellID:
                _foundMatch=True
                x_corrected, y_corrected = x_uncorrected + row["shift_x"], y_uncorrected + row["shift_y"]
                
        # keeps uncorrected values if no match is found                            
        if not _foundMatch:
            print("Did not find match for CellID #{} in ROI #{}".format(CellID,ROI))
            x_corrected,y_corrected = x_uncorrected,y_uncorrected
            self.foundMatch.append(False)                    
        else:
            self.foundMatch.append(True)
            
        return x_corrected, y_corrected

    def buildsVector(self,groupKeys,x,y,z):
        '''
        Builds vector from coordinates
        
        Parameters
        ----------
        groupKeys : list
            list of headers in the barcodes table
        x : float
            x coordinates
        y : float
            y coordinates
        z : float
            z coordinates
            
        Returns
        -------
        R : np array
            vector with coordinates.

        '''
        if self.ndims==3 and "zcentroidGauss" in groupKeys:
            R = np.column_stack((x,y,z,)) 
        else:
            R = np.column_stack((x,y,)) 

        return R

    def calculatesPWDsingleMask(self,ROI,CellID,groupKeys,x_uncorrected, y_uncorrected,z_uncorrected):
        '''
        Calculates PWD between barcodes detected in a given mask

        Parameters
        ----------
        ROI : string
            ROI used
        CellID: string
            ID of the cell
        x_uncorrected: float
            x coordinates uncorrected
        y_uncorrected: float
            y coordinates uncorrected
        z_uncorrected: float
            z coordinates uncorrected
            
        Returns
        -------
        Returns pairwise distance matrix between corrected barcodes

        '''

        if len(self.alignmentResultsTable)>0:                  
 
            # searches for local alignment shift for this mask in this ROI
            x_corrected, y_corrected = self.searchLocalShift(ROI,CellID,x_uncorrected,y_uncorrected)
              
            # applies local drift correction
            R = self.buildsVector(groupKeys,x_corrected, y_corrected,z_uncorrected )
                    
        else:
            # does not apply local drift correction
            R = self.buildsVector(groupKeys, x_uncorrected, y_uncorrected,z_uncorrected )

        return pairwise_distances(R)


    def buildsSCdistanceTable(self):
        '''
        iterates over all masks, calculates PWD for each mask, assigns them to SCdistanceTable

        Returns
        -------
        SCdistanceTable

        '''
        # sorts Table by cellID
        barcodeMapROI = self.barcodeMapROI
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")  # ROI data sorted by cellID
        
        self.initializeLists()
        
        # iterates over all cell masks in an ROI
        self.foundMatch=[]
        for key, group in tzip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):
            if key["CellID #"] > 1:  # excludes cellID 0 as this is background

                x_uncorrected, y_uncorrected= np.array(group["xcentroid"].data), np.array(group["ycentroid"].data)
                groupKeys, CellID, ROI = group.keys(), key["CellID #"], group["ROI #"].data[0]
                
                if self.ndims==3 and "zcentroidGauss" in group.keys():
                    z_uncorrected = np.array(group["zcentroidGauss"].data)
                else:
                    z_uncorrected = []

                PWD=self.calculatesPWDsingleMask(ROI,CellID,groupKeys,x_uncorrected, y_uncorrected,z_uncorrected)

                self.ROIs.append(group["ROI #"].data[0])
                self.cellID.append(key["CellID #"])
                self.nBarcodes.append(len(group))
                self.barcodeIDs.append(group["Barcode #"].data)
                self.buid.append(group["Buid"].data)
                self.p.append(PWD)
                self.cuid.append(str(uuid.uuid4()))  # creates cell unique identifier
                # print("CellID #={}, nBarcodes={}".format(key['CellID #'],len(group)))
                
        print("Local correction applied to {}/{} barcodes in ROI {}".format(np.nonzero(self.foundMatch)[0].shape[0],len(self.foundMatch),group["ROI #"].data[0]))
        if self.ndims==3 and "zcentroidGauss" in group.keys():
            print("Coordinates dimensions: 3")
        else:
            print("Coordinates dimensions: 2")
            
        SCdistanceTable = Table()
        SCdistanceTable["Cuid"] = self.cuid
        SCdistanceTable["ROI #"] = self.ROIs
        SCdistanceTable["CellID #"] = self.cellID
        SCdistanceTable["nBarcodes"] = self.nBarcodes
        SCdistanceTable["Barcode #"] = self.barcodeIDs
        SCdistanceTable["Buid"] = self.buid
        SCdistanceTable["PWDmatrix"] = self.p

        self.SCdistanceTable=SCdistanceTable
    
    def buildsdistanceMatrix(self, mode="mean"):
        """
        Builds pairwise distance matrix from a coordinates table

        Parameters
        ----------
        mode : string, optional
            The default is "mean": calculates the mean distance if there are several combinations possible.
            "min": calculates the minimum distance if there are several combinations possible.
            "last": keeps the last distance calculated
            
        Returns
        -------
        self.SCmatrix the single-cell PWD matrix
        self.meanSCmatrix the ensamble PWD matrix (mean of SCmatrix without nans)
        self.uniqueBarcodes list of unique barcodes
        
        """
        # [ builds SCdistanceTable ]
        self.buildsSCdistanceTable()
        print("Cells with barcodes found: {}".format(len(self.SCdistanceTable)))

        # [ builds SCmatrix ]

        numberMatrices = len(self.SCdistanceTable)  # z dimensions of SCmatrix
        uniqueBarcodes = np.unique(self.barcodeMapROI["Barcode #"].data)
        # number of unique Barcodes for xy dimensions of SCmatrix
        numberUniqueBarcodes = uniqueBarcodes.shape[0]
        SCmatrix = np.zeros((numberUniqueBarcodes, numberUniqueBarcodes, numberMatrices))
        SCmatrix[:] = np.NaN

        # loops over cell masks
        for iCell, scPWDitem in zip(range(numberMatrices), self.SCdistanceTable):
            barcodes2Process = scPWDitem["Barcode #"]
            
            # loops over barcodes detected in cell mask: barcode1
            for barcode1, ibarcode1 in zip(barcodes2Process, range(len(barcodes2Process))):
                indexBarcode1 = np.nonzero(uniqueBarcodes == barcode1)[0][0]
                
                # loops over barcodes detected in cell mask: barcode2
                for barcode2, ibarcode2 in zip(barcodes2Process, range(len(barcodes2Process))):
                    indexBarcode2 = np.nonzero(uniqueBarcodes == barcode2)[0][0]
                    
                    if barcode1 != barcode2:
                    
                        # attributes distance from the PWDmatrix field in the scPWDitem table
                        newdistance = scPWDitem["PWDmatrix"][ibarcode1][ibarcode2]
                        
                        # inserts value into SCmatrix
                        if mode == "last":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = newdistance
                        elif mode == "mean":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = np.nanmean(
                                [newdistance, SCmatrix[indexBarcode1][indexBarcode2][iCell],]
                            )
                        elif mode == "min":
                            SCmatrix[indexBarcode1][indexBarcode2][iCell] = np.nanmin(
                                [newdistance, SCmatrix[indexBarcode1][indexBarcode2][iCell],]
                            )

        self.SCmatrix = SCmatrix
        self.meanSCmatrix = np.nanmean(SCmatrix, axis=2)
        self.uniqueBarcodes = uniqueBarcodes

# =============================================================================
# FUNCTIONS
# =============================================================================


def loadsLocalAlignment(dataFolder):
    '''
    reads and returns localAlignmentTable, if it exists

    Parameters
    ----------
    dataFolder : folder()
        DESCRIPTION.

    Returns
    -------
    alignmentResultsTable : Table()
        DESCRIPTION.
    alignmentResultsTableRead : Boolean
        DESCRIPTION.

    '''
    localAlignmentFileName=dataFolder.outputFiles["alignImages"].split(".")[0] + "_localAlignment.dat"
    if os.path.exists(localAlignmentFileName):
        alignmentResultsTable= Table.read(localAlignmentFileName, format="ascii.ecsv")
        alignmentResultsTableRead=True
        print("LocalAlignment file loaded")
    else:
        print("\n\n*** Warning: could not find localAlignment: {}\n Proceeding with only global alignments...".format(localAlignmentFileName))
        alignmentResultsTableRead=False
        alignmentResultsTable=Table()

    return alignmentResultsTable, alignmentResultsTableRead

def loadsBarcodeMap(fileNameBarcodeCoordinates, ndims):
    '''
    Loads barcodeMap

    Parameters
    ----------
    fileNameBarcodeCoordinates : string
        filename with barcodeMap
    ndims : int
        either 2 or 3.

    Returns
    -------
    barcodeMap : Table()
    localizationDimension : int
        either 2 or 3.

    '''
    if os.path.exists(fileNameBarcodeCoordinates):
        barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
        if ndims==3 and "zcentroidGauss" in barcodeMap.keys():
            localizationDimension = 3
        else:
            localizationDimension = 2
    else:
        print("\n\n*** ERROR: could not find coordinates file: {}".format(fileNameBarcodeCoordinates))
        sys.exit()
        
    return  barcodeMap, localizationDimension

def buildsDictionaryErrorAlignmentMasks(param,dataFolder):
    """
    Builds and returns dictionary with error alignment block masks produced during the alignment process if 
    the 'blockAlignment' option was used

    Parameters
    ----------
    param : Parameters()
    dataFolder : folder()

    Returns
    -------
    dictErrorBlockMasks : dict

    """
    folder = dataFolder.outputFolders["alignImages"] 
    fileList = glob.glob(folder + os.sep + "*_errorAlignmentBlockMap.npy")

    # decodes files and builds dictionnary
    fileNameRegExp = param.param["acquisition"]["fileNameRegExp"]
    fileNameRegExp = fileNameRegExp.split('.')[0]
    listRE = [re.search(fileNameRegExp,os.path.basename(x).split("_errorAlignmentBlockMap.npy")[0]) for x in fileList]
    
    dictErrorBlockMasks = dict()

    for file, regExp in zip(fileList,listRE):
        if 'ROI:'+str(int(regExp['roi'])) not in dictErrorBlockMasks.keys():
            dictErrorBlockMasks['ROI:'+str(int(regExp['roi']))]= {}
        if 'barcode:'+regExp['cycle'].split('RT')[-1] not in dictErrorBlockMasks.keys():
            newMask = np.load(file)
            dictErrorBlockMasks['ROI:'+str(int(regExp['roi']))]['barcode:'+regExp['cycle'].split('RT')[-1]]= newMask

    return dictErrorBlockMasks

def buildsPWDmatrix(param,
    currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize=0.1, logNameMD="log.md", ndims=2
):

    # Loads localAlignment if it exists
    alignmentResultsTable, alignmentResultsTableRead = loadsLocalAlignment(dataFolder)

    # Loads coordinate Tables      
    barcodeMap, localizationDimension = loadsBarcodeMap(fileNameBarcodeCoordinates, ndims)


    # Builds dictionnary with filenames of errorAlignmentBlockMasks for each ROI and each barcode
    dictErrorBlockMasks = buildsDictionaryErrorAlignmentMasks(param,dataFolder)
   
    # processes tables 
    barcodeMapROI = barcodeMap.group_by("ROI #")
    numberROIs = len(barcodeMapROI.groups.keys)
    print("\nROIs detected: {}".format(numberROIs))
    
    # loops over ROIs
    filesinFolder = glob.glob(currentFolder + os.sep + "*.tif")
    SCmatrixCollated, uniqueBarcodes, processingOrder = [], [], 0

    for ROI in range(numberROIs):
        nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

        print("Loading masks and pre-processing barcodes for ROI# {}".format(nROI))

        barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

        # finds file with cell masks
        fileList2Process = [
            file
            for file in filesinFolder
            if file.split("_")[-1].split(".")[0] == "ch00"
            and "DAPI" in os.path.basename(file).split("_")
            and int(os.path.basename(file).split("_")[3]) == nROI
        ]

        if len(fileList2Process) > 0:
            
            # loads file with cell masks
            fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
            fullFileNameROImasks = os.path.dirname(fileNameBarcodeCoordinates) + os.sep + fileNameROImasks
            if os.path.exists(fullFileNameROImasks):
                Masks = np.load(fullFileNameROImasks)

                # Assigns barcodes to Masks for a given ROI
                cellROI = cellID(param,barcodeMapSingleROI, Masks, ROI,ndims=localizationDimension)
                cellROI.ndims=ndims
                cellROI.nROI = nROI
                
                if alignmentResultsTableRead:
                    cellROI.alignmentResultsTable = alignmentResultsTable
                    
                cellROI.dictErrorBlockMasks = dictErrorBlockMasks
                cellROI.alignmentResultsTableRead = alignmentResultsTableRead
                
                # finds what barcodes are in each cell mask    
                cellROI.alignByMasking()

                # builds the single cell distance Matrix
                cellROI.buildsdistanceMatrix("min")  # mean min last

                print(
                    "ROI: {}, N cells assigned: {} out of {}\n".format(ROI, cellROI.NcellsAssigned-1, cellROI.numberMasks)
                )

                uniqueBarcodes = cellROI.uniqueBarcodes

                # saves Table with results per ROI
                cellROI.SCdistanceTable.write(
                    outputFileName + "_order:" + str(processingOrder) + "_ROI:" + str(nROI) + ".ecsv",
                    format="ascii.ecsv",
                    overwrite=True,
                )

                if len(SCmatrixCollated) > 0:
                    SCmatrixCollated = np.concatenate((SCmatrixCollated, cellROI.SCmatrix), axis=2)
                else:
                    SCmatrixCollated = cellROI.SCmatrix
                del cellROI

                processingOrder += 1

            else:
                print(
                    "Error, no DAPI mask file found for ROI: {}, segmentedMasks: {}\n".format(
                        nROI, fileNameBarcodeCoordinates
                    )
                )
                print("File I was searching for: {}".format(fullFileNameROImasks))
                print("Debug: ")
                for file in filesinFolder:
                    if (
                        file.split("_")[-1].split(".")[0] == "ch00"
                        and "DAPI" in file.split("_")
                        and int(os.path.basename(file).split("_")[3]) == nROI
                    ):
                        print("Hit found!")
                    print(
                        "fileSplit:{}, DAPI in filename: {}, ROI: {}".format(
                            file.split("_")[-1].split(".")[0],
                            "DAPI" in os.path.basename(file).split("_"),
                            int(os.path.basename(file).split("_")[3]),
                        )
                    )

    # saves output
    np.save(outputFileName + "_HiMscMatrix.npy", SCmatrixCollated)
    np.savetxt(outputFileName + "_uniqueBarcodes.ecsv", uniqueBarcodes, delimiter=" ", fmt="%d")

    # plots outputs
    
    # adapts clim depending on whether 2 or 3 dimensions are use for the barcode localizations
    if localizationDimension==2:
        clim=1.6
    else:
        clim=2.2
        
    plotMatrix(
        SCmatrixCollated, uniqueBarcodes, pixelSize, numberROIs, outputFileName, logNameMD, mode="median", clim=clim, cm='terrain'
    )  # need to validate use of KDE. For the moment it does not handle well null distributions

    plotDistanceHistograms(SCmatrixCollated, pixelSize, outputFileName, logNameMD)


def processesPWDmatrices(param, log1, session1):
    sessionName = "buildsPWDmatrix"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText("\n===================={}====================\n".format(sessionName))
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(log1.fileNameMD, "## {}\n".format(sessionName), "a")

    for currentFolder in dataFolder.listFolders:
        # filesFolder=glob.glob(currentFolder+os.sep+'*.tif')
        dataFolder.createsFolders(currentFolder, param)
        log1.report("-------> Processing Folder: {}".format(currentFolder))


        fileNameBarcodeCoordinates = dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"
        
        # 2D
        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]
        log1.report("-------> 2D processing: {}".format(outputFileName ))        
        
        if "pixelSizeXY" in param.param['acquisition'].keys():
            pixelSize = param.param['acquisition']['pixelSizeXY']
        else:
            pixelSize = 0.1
            
        buildsPWDmatrix(
            param,currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize, log1.fileNameMD,
        )

        # 3D
        outputFileName = dataFolder.outputFiles["buildsPWDmatrix"]+"_3D"
        log1.report("-------> 3D processing: {}".format(outputFileName ))        

        if ("pixelSizeZ" in param.param['acquisition'].keys()) and ("pixelSizeXY" in param.param['acquisition'].keys()):
            pixelSizeXY = param.param['acquisition']['pixelSizeXY']
            pixelSizeZ = param.param['acquisition']['pixelSizeZ']
            
            # need to solve the issue of voxelsize...
            # pixelSize = [pixelSizeXY,pixelSizeXY,pixelSizeZ]
            pixelSize = pixelSizeXY
        else:
            pixelSize = 0.1            
            
        buildsPWDmatrix(
            param,currentFolder, fileNameBarcodeCoordinates, outputFileName, dataFolder, pixelSize, log1.fileNameMD, ndims=3
        )

        # loose ends
        session1.add(currentFolder, sessionName)

        log1.report("HiM matrix in {} processed".format(currentFolder), "info")


