#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

test fitting barcode spots to masks

"""

# =============================================================================
# IMPORTS
# =============================================================================


import glob, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from imageProcessing import Image
from fileManagement import folders
from fileManagement import session, writeString2File

from astropy.table import Table, vstack, Column
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clip, sigma_clipped_stats

from photutils.segmentation import SegmentationImage
from sklearn.metrics import pairwise_distances

# =============================================================================
# CLASSES
# =============================================================================


class cellID:
    def __init__(self, barcodeMapROI, Masks, ROI):
        self.barcodeMapROI = barcodeMapROI
        self.Masks = Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask = 0

        self.SegmentationMask = SegmentationImage(self.Masks)
        self.numberMasks = self.SegmentationMask.nlabels
        self.ROI = ROI

        self.barcodesinMask = dict()
        for mask in range(self.numberMasks + 1):
            self.barcodesinMask["maskID_" + str(mask)] = []

    def visualize(self):

        imageBarcodes = np.zeros([2048, 2048])
        MasksBarcodes = Masks
        R = []

        for i in range(len(self.barcodeMapROI.groups[0])):
            y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
            x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
            barcodeID = self.barcodeMapROI.groups[0]["Barcode #"][i]
            imageBarcodes[x_int][y_int] = barcodeID
            MasksBarcodes[x_int][y_int] += 20 * barcodeID
            R.append([y_int, x_int, barcodeID])

        # Shows results
        Ra = np.array(R)
        plt.imshow(Masks, origin="lower", cmap="jet")
        plt.scatter(Ra[:, 0], Ra[:, 1], s=5, c=Ra[:, 2], alpha=0.5)

    def alignByMasking(self):
        # [ Assigns barcodes to masks and creates <NbarcodesinMask> ]
        NbarcodesinMask = np.zeros(self.numberMasks + 2)
        print("ROI:{}".format(self.ROI))
        for i in range(len(self.barcodeMapROI.groups[0])):
            y_int = int(self.barcodeMapROI.groups[0]["xcentroid"][i])
            x_int = int(self.barcodeMapROI.groups[0]["ycentroid"][i])
            # barcodeID = self.barcodeMapROI.groups[ROI]['Barcode #'][i]
            maskID = self.Masks[x_int][y_int]
            self.barcodeMapROI["CellID #"][i] = maskID
            if maskID > 0:
                NbarcodesinMask[maskID] += 1
                self.barcodesinMask["maskID_" + str(maskID)].append(i)

        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask > 0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned
        self.NbarcodesinMask = NbarcodesinMask

    def buildsdistanceMatrix(self, mode="mean"):
        """
        

        
        """
        print("building distance matrix")
        barcodeMapROI = self.barcodeMapROI

        # [ builds SCdistanceTable ]
        barcodeMapROI_cellID = barcodeMapROI.group_by("CellID #")  # ROI data sorted by cellID
        ROIs, cellID, nBarcodes, barcodeIDs, p = [], [], [], [], []

        for key, group in zip(barcodeMapROI_cellID.groups.keys, barcodeMapROI_cellID.groups):

            if key["CellID #"] > 1:  # excludes cellID 0 as this is background
                R = np.column_stack((np.array(group["xcentroid"].data), np.array(group["ycentroid"].data),))
                ROIs.append(group["ROI #"].data[0])
                cellID.append(key["CellID #"])
                nBarcodes.append(len(group))
                barcodeIDs.append(group["Barcode #"].data)
                p.append(pairwise_distances(R))
                # print("CellID #={}, nBarcodes={}".format(key['CellID #'],len(group)))

        SCdistanceTable = Table()  # [],names=('CellID', 'barcode1', 'barcode2', 'distances'))
        SCdistanceTable["ROI #"] = ROIs
        SCdistanceTable["CellID #"] = cellID
        SCdistanceTable["nBarcodes"] = nBarcodes
        SCdistanceTable["Barcode #"] = barcodeIDs
        SCdistanceTable["PWDmatrix"] = p

        self.SCdistanceTable = SCdistanceTable

        print("Cells with barcodes found: {}".format(len(SCdistanceTable)))

        # [ builds SCmatrix ]
        numberMatrices = len(SCdistanceTable)  # z dimensions of SCmatrix
        uniqueBarcodes = np.unique(barcodeMapROI["Barcode #"].data)
        numberUniqueBarcodes = uniqueBarcodes.shape[0]  # number of unique Barcodes for xy dimensions of SCmatrix
        SCmatrix = np.zeros((numberUniqueBarcodes, numberUniqueBarcodes, numberMatrices))
        SCmatrix[:] = np.NaN

        for iCell, scPWDitem in zip(range(numberMatrices), SCdistanceTable):
            barcodes2Process = scPWDitem["Barcode #"]
            for barcode1, ibarcode1 in zip(barcodes2Process, range(len(barcodes2Process))):
                indexBarcode1 = np.nonzero(uniqueBarcodes == barcode1)[0][0]
                for barcode2, ibarcode2 in zip(barcodes2Process, range(len(barcodes2Process))):
                    indexBarcode2 = np.nonzero(uniqueBarcodes == barcode2)[0][0]
                    if barcode1 != barcode2:
                        newdistance = scPWDitem["PWDmatrix"][ibarcode1][ibarcode2]
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

# loads coordinate file

rootFolder = "/home/marcnol/data/Experiment_20/Embryo_1"
fullFolder = rootFolder + "/rawData/segmentedObjects/"
fileNameBarcodeCoordinates = fullFolder + "segmentedObjects_barcode.dat"

"""
rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
fullFolder=rootFolder +'/rawData/segmentedObjects/'
fileNameBarcodeCoordinates =  fullFolder+'segmentedObjects_barcode.dat'
fileNameROImasks = fullFolder +'scan_001_DAPI_018_ROI_converted_decon_ch00_Masks.npy'

rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
fullFolder=rootFolder +'/rawImages/segmentedObjects/'
fileNameBarcodeCoordinates =  fullFolder+'segmentedObjects_barcode.dat'
fileNameROImasks = fullFolder +'scan_001_DAPI_017_ROI_converted_decon_ch00_Masks.npy'
"""

# Processes Tables
barcodeMap = Table.read(fileNameBarcodeCoordinates, format="ascii.ecsv")
barcodeMapROI = barcodeMap.group_by("ROI #")

SCmatrixCollated = []
for ROI in range(len(barcodeMapROI.groups.keys)):
    nROI = barcodeMapROI.groups.keys[ROI][0]  # need to iterate over the first index

    print("ROIs detected: {}".format(barcodeMapROI.groups.keys))

    barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[ROI]

    # finds file for masks
    filesFolder = glob.glob(rootFolder + "/rawData/" + "*.tif")

    fileList2Process = [
        file
        for file in filesFolder
        if file.split("_")[-1].split(".")[0] == "ch00"
        and "DAPI" in file.split("_")
        and int(os.path.basename(file).split("_")[3]) == nROI
    ]

    if len(fileList2Process) > 0:

        # loads Masks
        fileNameROImasks = os.path.basename(fileList2Process[0]).split(".")[0] + "_Masks.npy"
        Masks = np.load(fullFolder + fileNameROImasks)

        # Assigns barcodes to Masks for a given ROI
        cellROI = cellID(barcodeMapSingleROI, Masks, ROI)

        cellROI.alignByMasking()

        cellROI.buildsdistanceMatrix("min")  # mean min last

        print("ROI: {}, N cells assigned: {} out of {}".format(ROI, cellROI.NcellsAssigned, cellROI.numberMasks))

        uniqueBarcodes = cellROI.uniqueBarcodes

        if len(SCmatrixCollated) > 0:
            SCmatrixCollated = np.concatenate((SCmatrixCollated, cellROI.SCmatrix), axis=2)
        else:
            SCmatrixCollated = cellROI.SCmatrix
        del cellROI

#%%
pixelSize = 0.1
meanSCmatrix = pixelSize * np.nanmedian(SCmatrixCollated, axis=2)
print("ROIs detected: {}".format(barcodeMapROI.groups.keys))

fig = plt.figure()
pos = plt.imshow(meanSCmatrix, cmap="seismic")
plt.xlabel("barcode #")
plt.ylabel("barcode #")
plt.title("PWD matrix" + " | n=" + str(SCmatrixCollated.shape[2]) + " | ROIs=" + str(len(barcodeMapROI.groups.keys)))
plt.xticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
plt.yticks(np.arange(SCmatrixCollated.shape[0]), uniqueBarcodes)
cbar = plt.colorbar(pos)
cbar.minorticks_on()
cbar.set_label("distance, um")
plt.clim(0, 1.4)

#%%

NplotsX = 3
sizeX, sizeY = NplotsX * 4, 4
fig, (ax1, ax2, ax3) = plt.subplots(figsize=(sizeX, sizeY), ncols=NplotsX)

pos1 = ax1.hist(pixelSize * SCmatrixCollated[0, 1, :])
pos2 = ax2.hist(pixelSize * SCmatrixCollated[0, 2, :])
pos3 = ax3.hist(pixelSize * SCmatrixCollated[1, 2, :])
plt.xlabel("distance, um")
plt.ylabel("counts")


#%%
"""
rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset'
fullFolder=rootFolder +'/rawImages/segmentedObjects/'
fileNameBarcodeCoordinates =  fullFolder+'segmentedObjects_barcode.dat'
fileNameROImasks = fullFolder +'scan_001_DAPI_017_ROI_converted_decon_ch00_Masks.npy'

barcodeMap = Table.read(fileNameBarcodeCoordinates,format='ascii.ecsv')
Masks=np.load(fileNameROImasks)
# loads Masks

# Processes Tables 
ROI=0
barcodeMapROI=barcodeMap.group_by('ROI #').groups[ROI]

print('ROIs detected: {}'.format(barcodeMapROI.groups.keys))

# Assigns barcodes to Masks for a given ROI
cellROI2 = cellID(barcodeMapROI,Masks)

cellROI2.alignByMasking()
print('N cells assigned: {} out of {}'.format(cellROI.NcellsAssigned,cellROI.numberMasks))

cellROI2.buildsdistanceMatrix('min') # mean min last

# collates SC matrices
SCmatrixCollated=np.concatenate((cellROI.SCmatrix,cellROI2.SCmatrix),axis=2)
"""
