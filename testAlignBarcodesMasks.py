barcodesinMask=dict()#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:23:36 2020

@author: marcnol

test fitting barcode spots to masks

"""

# =============================================================================
# IMPORTS
# =============================================================================


import glob,os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from imageProcessing import Image
from fileManagement import folders
from fileManagement import session,writeString2File

from astropy.table import Table, vstack, Column
from astropy.visualization import SqrtStretch,simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.segmentation import SegmentationImage

# =============================================================================
# CLASSES
# =============================================================================

class cellID():
    def __init__(self,barcodeMapROI,Masks):
        self.barcodeMapROI=barcodeMapROI
        self.Masks=Masks
        self.NcellsAssigned = 0
        self.NcellsUnAssigned = 0
        self.NbarcodesinMask=0

        self.SegmentationMask=SegmentationImage(self.Masks)
        self.numberMasks=self.SegmentationMask.nlabels
        
        self.barcodesinMask=dict()
        for mask in range(self.numberMasks):
            self.barcodesinMask['maskID_'+str(mask)]=[]        

    def visualize(self):

        imageBarcodes=np.zeros([2048,2048])
        MasksBarcodes=Masks
        R=[]

        for i in range(len(self.barcodeMapROI.groups[ROI])):
            y_int=int(self.barcodeMapROI.groups[ROI]['xcentroid'][i])
            x_int=int(self.barcodeMapROI.groups[ROI]['ycentroid'][i])
            barcodeID = self.barcodeMapROI.groups[ROI]['Barcode #'][i]
            imageBarcodes[x_int][y_int] = barcodeID
            MasksBarcodes[x_int][y_int] += 20*barcodeID
            R.append([y_int,x_int,barcodeID])

        # Shows results
        Ra=np.array(R)
        plt.imshow(Masks, origin='lower', cmap='jet')
        plt.scatter(Ra[:,0],Ra[:,1],s=5,c=Ra[:,2],alpha=0.5)
      
    def alignByMasking(self):
        NbarcodesinMask=np.zeros(self.numberMasks)

        for i in range(len(self.barcodeMapROI.groups[ROI])):
            y_int=int(self.barcodeMapROI.groups[ROI]['xcentroid'][i])
            x_int=int(self.barcodeMapROI.groups[ROI]['ycentroid'][i])
            #barcodeID = self.barcodeMapROI.groups[ROI]['Barcode #'][i]
            maskID = self.Masks[x_int][y_int]
            self.barcodeMapROI['CellID #'][i]=maskID
            if maskID>0:
                NbarcodesinMask[maskID]+=1
                self.barcodesinMask['maskID_'+str(maskID)].append(i)
            
        self.NcellsAssigned = np.count_nonzero(NbarcodesinMask>0)
        self.NcellsUnAssigned = self.numberMasks - self.NcellsAssigned
        self.NbarcodesinMask=NbarcodesinMask
        
    def buildsdistanceMatrix(self):
        print('building distance matrix')
        
        for icell in range(self.numberMasks):
            if self.NbarcodesinMask[icell]>1:
                X=[] # matrix with barcodes in cell
                pdist(X, metric='euclidian')
        
        
        
        
# =============================================================================
# FUNCTIONS
# =============================================================================

# loads coordinate file
rootFolder='/home/marcnol/data/Experiment_15/Embryo_006_ROI18'
fullFolder=rootFolder +'/rawData/segmentedObjects/'
fileNameBarcodeCoordinates =  fullFolder+'segmentedObjects_barcode.dat'
fileNameROImasks = fullFolder +'scan_001_DAPI_018_ROI_converted_decon_ch00_Masks.npy'

barcodeMap = Table.read(fileNameBarcodeCoordinates,format='ascii.ecsv')
Masks=np.load(fileNameROImasks)
# loads Masks


# Processes Tables 
ROI=0
barcodeMapROI=barcodeMap.group_by('ROI #')
print('ROIs detected: {}'.format(barcodeMapROI.groups.keys))

# Assigns barcodes to Masks for a given ROI
cellROI = cellID(barcodeMapROI,Masks)

cellROI.alignByMasking()
print('N cells assigned: {} out of {}'.format(cellROI.NcellsAssigned,cellROI.numberMasks))

cellROI.visualize()



