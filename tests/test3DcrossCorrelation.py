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
filename2 = rootFolder+os.sep+"scan_001_RT41_001_ROI_converted_decon_ch00.tif"

print("Reading files: \n{}\n{}".format(filename1,filename2))
img1=io.imread(filename1).squeeze()
img2=io.imread(filename2).squeeze()

#%%
print("Calculating shifts...")
upsample_factor=100

shift, error, diffphase = phase_cross_correlation(img1, img2, upsample_factor=upsample_factor)

img2_corrected = shiftImage(img2, shift)

#%%

fig, axes = plt.subplots(1,3)
fig.set_size_inches((10, 5))
ax=axes.ravel()
# fig.suptitle(title)
subVolume=((0,60),(900,1100),(900,900))
cmap='jet'
ax[0].imshow(img1[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)
ax[1].imshow(img2[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)
ax[2].imshow(img2_corrected[subVolume[0][0]:subVolume[0][1],subVolume[1][0]:subVolume[1][1],subVolume[2][0]].transpose(),cmap=cmap)