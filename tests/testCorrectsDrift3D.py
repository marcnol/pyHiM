#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 22:03:28 2021

@author: marcnol

- provide zcentroid also for 2D fitting routines
- make routine that segments fiducials so they can be used to compare accuracy of drift correction.
- finish coding in section in alignBarcodesMask

"""

import os
from astropy.table import Table
from skimage import io
import numpy as np
import matplotlib.pylab as plt
from sklearn.metrics import pairwise_distances

folder = "/home/marcnol/data/Embryo_debug_dataset/test_dataset/segmentedObjects"
file =folder+os.sep+"segmentedObjects_barcode.dat"
barcodeMap = Table.read(file, format="ascii.ecsv")

imageFile='/home/marcnol/data/Embryo_debug_dataset/test_dataset/alignImages/'+'scan_001_RT27_001_ROI_converted_decon_ch01_2d_registered.npy'
img = np.load(imageFile)

#%% selects localizations for a given ROI and barcode
barcodeMapROI = barcodeMap.group_by("ROI #")

numberROIs = len(barcodeMapROI.groups.keys)
iROI = barcodeMapROI.groups.keys[0][0]  # need to iterate over the first index
i=0
barcode =barcodeMapROI.groups[0]["Barcode #"][i]
ROI = barcodeMapROI.groups[0]["ROI #"][i]
Barcode =  31

# load dict of 3D shifts
localAlignmentFileName ='/home/marcnol/data/Embryo_debug_dataset/test_dataset/alignImages/'+"alignImages_block3Dalignment.dat"
alignmentResultsTable = Table.read(localAlignmentFileName, format="ascii.ecsv")
blockSizeXY = alignmentResultsTable[0]["blockXY"]
zxy_uncorrectedList,zxy_correctedList=[],[]
foundMatch=[]
for i in range(len(barcodeMapROI.groups[0])): # i is the index of the barcode in barcodeMapROI
    barcode = barcodeMapROI.groups[0]["Barcode #"][i]
    roi = barcodeMapROI.groups[0]["ROI #"][i]

    if roi==ROI and barcode==Barcode:
        x_corrected = barcodeMapROI.groups[0]["ycentroid"][i] # control inversion between x-y
        y_corrected = barcodeMapROI.groups[0]["xcentroid"][i]
        z_corrected = 0.0

        zxy_uncorrected = [z_corrected, x_corrected , y_corrected]

        zxy_uncorrectedList.append(zxy_uncorrected)

        # correct 3D shift
        _foundMatch = False
        zxyBlock = [np.floor(a/blockSizeXY).astype(int) for a in zxy_uncorrected]
        for row in alignmentResultsTable:
            if row["ROI #"] == ROI and row["label"] == "RT" + str(barcode) and row["block_i"] == zxyBlock[1] and row["block_j"] == zxyBlock[2]:
                _foundMatch = True
                shifts = [row["shift_z"],row["shift_x"],row["shift_y"]]
                zxy_corrected = [a+shift for a,shift in zip(zxy_uncorrected,shifts)]
                zxy_correctedList.append(zxy_corrected)
        if not _foundMatch:
            print("Did not find match for ROI #{} barcode #{}, ZXY={}".format(ROI, barcode,zxyBlock))
            zxy_correctedList.append(zxy_uncorrected)
            foundMatch.append(False)
        else:
            foundMatch.append(True)
        # verify correction


#%% plots image and localizatons

y = [x[1] for x in zxy_uncorrectedList]
x = [x[2] for x in zxy_uncorrectedList]

yCorrected = [x[1] for x in zxy_correctedList]
xCorrected = [x[2] for x in zxy_correctedList]

plt.imshow(img,cmap="Greys",vmax=1000)
plt.plot(xCorrected,yCorrected,'o',color='r',alpha=.8)
plt.plot(x,y,'+',color='y',alpha=.9)

#%% plots histograms of differences of localizations between barcode 31 and barcode 27

# get localizations for reference barcode: RT27
RefenceBarcode =  27
zxy_uncorrectedList_reference=[]
for i in range(len(barcodeMapROI.groups[0])): # i is the index of the barcode in barcodeMapROI
    barcode = barcodeMapROI.groups[0]["Barcode #"][i]
    roi = barcodeMapROI.groups[0]["ROI #"][i]

    if roi==ROI and barcode==RefenceBarcode :
        x_corrected = barcodeMapROI.groups[0]["ycentroid"][i] # control inversion between x-y
        y_corrected = barcodeMapROI.groups[0]["xcentroid"][i]
        z_corrected = 0.0
        zxy_uncorrected = [z_corrected, x_corrected , y_corrected]
        zxy_uncorrectedList_reference.append(zxy_uncorrected)

y_reference = [x[1] for x in zxy_uncorrectedList_reference]
x_reference = [x[2] for x in zxy_uncorrectedList_reference]

# plots differences between localizations
R = np.column_stack((x, y,))
R_corrected = np.column_stack((xCorrected, yCorrected,))
R_reference = np.column_stack((x_reference, y_reference,))
PWD = pairwise_distances(R,R_reference )
PWD_corrected = pairwise_distances(R_corrected,R_reference )

deltaX_uncorrected=[]
deltaX_corrected=[]
for i in range(len(x_reference)):
    minDistance_uncorrected = np.min(PWD[i,:])
    deltaX_uncorrected.append(minDistance_uncorrected)
    minDistance_corrected = np.min(PWD_corrected[i,:])
    deltaX_corrected.append(minDistance_corrected )

fig = plt.figure()
fig.subplots_adjust(top=0.8)
ax1 = fig.add_subplot(211)
ax1.set_xlabel('min distance, px')
ax1.set_ylabel('counts')
ax1.hist(deltaX_uncorrected,bins=20,range=(0, 3),alpha=0.5)

# ax2 = fig.add_axes([0.15, 0.1, 0.7, 0.3])
ax2 = fig.add_subplot(212)
ax2.set_xlabel('min distance, px')
ax2.set_ylabel('counts')
ax2.hist(deltaX_corrected,bins=20,range=(0, 3),alpha=0.5)
