#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:49:13 2021

@author: marcnol
"""




import numpy as np
import matplotlib.pylab as plt


from skimage import io
import os
import numpy as np
import matplotlib.pylab as plt
from imageProcessing.imageProcessing import _reinterpolatesFocalPlane, imageShowWithValues
# from astropy.stats import SigmaClip
from scipy.stats import sigmaclip

from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift as shiftImage


rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18"
filename1 = rootFolder+os.sep+"scan_001_RT27_001_ROI_converted_decon_ch00.tif"
filename2 = rootFolder+os.sep+"scan_001_RT29_001_ROI_converted_decon_ch00.tif"

print("Reading files: \n{}\n{}".format(filename1,filename2))
img1=io.imread(filename1).squeeze()
img2=io.imread(filename2).squeeze()

#%%
print("Calculating shifts...")
upsample_factor=100

img1=img1/img1.max()
img2=img2/img2.max()

shift, error, diffphase = phase_cross_correlation(img1, img2, upsample_factor=upsample_factor)

img2_corrected = shiftImage(img2, shift)
img2_corrected=img2_corrected/img2_corrected.max()

#%% plots 2d projection
fig, axes = plt.subplots(1,2)
fig.set_size_inches((10, 5))
ax=axes.ravel()
# fig.suptitle(title)
subVolume=((0,60),(0,2048),(0,2048))
subVolume=((0,60),(1000,1250),(800,1000))

cmap='Greys'
ax[0].imshow(np.sum(img1[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]],axis=0),cmap=cmap,vmax=1)
# ax[1].imshow(img2_corrected[subVolume[0][0],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]].transpose(),cmap=cmap,vmax=1000)
ax[1].imshow(np.sum(img2_corrected[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]],axis=0),cmap=cmap,vmax=1)


#%% plots 2d single Z plane
fig, axes = plt.subplots(1,2)
fig.set_size_inches((10, 5))
ax=axes.ravel()
# fig.suptitle(title)
subVolume=((40,60),(0,2048),(0,2048))
cmap='jet'
ax[0].imshow((img1[subVolume[0][0],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]]),cmap=cmap,vmax=.11)
ax[1].imshow(img2_corrected[subVolume[0][0],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]].transpose(),cmap=cmap,vmax=1000)
# ax[1].imshow((img2[subVolume[0][0],subVolume[1][0]:subVolume[1][1],subVolume[2][0]:subVolume[2][1]]),cmap=cmap,vmax=.041)

#%% plots 3D

fig, axes = plt.subplots(1,3)
fig.set_size_inches((10, 5))
ax=axes.ravel()
# fig.suptitle(title)
subVolume=((0,60),(1000,1250),(1200,1000))
cmap='Greys'
ax[0].imshow(img1[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)
diffMatrix=img2-img1
ax[1].imshow(img2[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)
diffMatrixCorrected=img2_corrected-img1
ax[2].imshow(img2_corrected[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)



