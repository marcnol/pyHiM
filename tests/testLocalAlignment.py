#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:10:39 2020

@author: marcnol
"""
import numpy as np
import matplotlib.pylab as plt
from stardist import random_label_cmap
from alignImages import align2ImagesCrossCorrelation
from scipy.ndimage import shift as shiftImage
from localDriftCorrection import alignsSubVolumes

np.random.seed(6)
lbl_cmap = random_label_cmap()

alignedFilesFolder='/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/alignImages/'
referenceFileName='scan_001_RT27_001_ROI_converted_decon_ch00_2d_registered.npy'
barcode29FileName='scan_001_RT29_001_ROI_converted_decon_ch00_2d_registered.npy'
barcode37FileName='scan_001_RT37_001_ROI_converted_decon_ch00_2d_registered.npy'

segmentedFilesFolder='/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/segmentedObjects/'
fullFileNameROImasks = 'scan_006_DAPI_001_ROI_converted_decon_ch00_Masks.npy'

Masks = np.load(segmentedFilesFolder+fullFileNameROImasks)
maskSize = Masks.shape


imageReference  = np.load(alignedFilesFolder+referenceFileName)
# imageReference = imReference.removesBackground2D(normalize=True)
imageBarcode  = np.load(alignedFilesFolder+barcode29FileName)

#%%

fig = plt.figure()
fig.set_size_inches((5, 5))
pos = plt.imshow(imageReference, cmap="magma", origin="lower")
plt.imshow(Masks, origin="lower", cmap=lbl_cmap, alpha=0.1)
fig.colorbar(pos)
plt.axis("off")

#%%
iMask=6
bezel=20

shift, subVolumeReference, subVolume, subVolumeCorrected =  alignsSubVolumes(imageReference, imageBarcode, Masks, bezel=bezel, iMask=iMask)
print("Shift = {}".format(shift))
# fig=plt.figure()
# fig.set_size_inches((5, 5))
# plt.imshow(image1_adjusted,cmap = 'Blues', alpha = 1)

                                 fig=plt.figure()
fig.set_size_inches((5, 5))
plt.imshow(subVolumeReference,cmap = 'Blues', alpha = 0.5)
plt.imshow(subVolume,cmap = 'Reds', alpha = 0.5)

fig=plt.figure()
fig.set_size_inches((5, 5))
plt.imshow(subVolumeReference,cmap = 'Blues', alpha = 0.5)
plt.imshow(subVolumeCorrected,cmap = 'Reds', alpha = 0.5)

