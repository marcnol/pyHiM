#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 15:28:59 2020

@author: marcnol

# integrate spot and dapi segmentation into pipeline

# sum up all 2D registered RT images to see if I can see the 'region' labeled as another parameter
to better define barcode groups rather than just using DAPI.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters, morphology, measure, color, exposure
from skimage.segmentation import watershed, clear_border
from scipy import ndimage as ndi
import matplotlib.patches as mpatches
from imageProcessing import Image
from fileManagement import loadJSON


#%% segments using normal morphological operations from skimage
rootDir = "/home/marcnol/Documents/Images/Embryo_debug_dataset"

file = "scan_004_RT19_017_ROI_converted_decon_ch01.tif"
fileName = rootDir + "/rawImages/alignImages/" + file.split(".")[0] + "_2d_registered.npy"
im = np.load(fileName)
im = exposure.rescale_intensity(im, out_range=(0, 1))

threshold = filters.threshold_otsu(im)
# binary= im > threshold*1.0
binary = morphology.closing(im > threshold * 2, morphology.square(3))
cleared = clear_border(binary)

label_image = measure.label(cleared)

# image_label_overlay=color.label2rgb(label_image,image=im,bg_label=0)
image_label_overlay = color.label2rgb(label_image, bg_label=0)

# plt.imshow(image_label_overlay)
# segmentation = watershed(threshold, markers)

fig, ax = plt.subplots(figsize=(10, 6))
ax.imshow(image_label_overlay, alpha=1)

for region in measure.regionprops(label_image):
    # take regions with large enough areas
    if region.area >= 25 and region.area < 500:
        # draw rectangle around segmented coins
        minr, minc, maxr, maxc = region.bbox
        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr, fill=False, edgecolor="red", linewidth=1,)
        ax.add_patch(rect)


ax.set_axis_off()
plt.tight_layout()
plt.show()

#%% Will use astropy and photutils to get point sources

from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import DAOStarFinder, CircularAperture, find_peaks, detect_sources
from photutils import detect_threshold, deblend_sources
from photutils import background
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

rootDir = "/home/marcnol/Documents/Images/Embryo_debug_dataset"

file = "scan_004_RT18_017_ROI_converted_decon_ch01.tif"
fileName = rootDir + "/rawImages/alignImages/" + file.split(".")[0] + "_2d_registered.npy"

im = np.load(fileName)

brightest = 1000
threshold_over_std = 1.8
fwhm = 3.0

# removes background
# mean, median, std = sigma_clipped_stats(im, sigma=3.0)
# print((mean, median, std))

sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(im, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
threshold_image = bkg.background + (1.0 * bkg.background_rms)  # background-only error image

im1_bkg_substracted = im - bkg.background
mean, median, std = sigma_clipped_stats(im1_bkg_substracted, sigma=3.0)

# estimates sources
daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_over_std * std, brightest=brightest)
sources = daofind(im1_bkg_substracted)
"""for col in sources.colnames:  
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)  
"""

# show results
fig = plt.figure()
fig.set_size_inches((50, 50))
positions = np.transpose((sources["xcentroid"] + 0.5, sources["ycentroid"] + 0.5))
apertures = CircularAperture(positions, r=4.0)
# norm = ImageNormalize(stretch=SqrtStretch())
norm = simple_norm(im, "sqrt", percent=99.9)
plt.imshow(im1_bkg_substracted, cmap="Greys", origin="lower", norm=norm)
apertures.plot(color="blue", lw=1.5, alpha=0.5)
# apertures.plot(color='#0547f9', lw=1.5)
plt.xlim(0, im.shape[1] - 1)
plt.ylim(0, im.shape[0] - 1)

#%% will use astropy and photutils to segment DAPI images!


file = "scan_001_DAPI_017_ROI_converted_decon_ch00.tif"
fileName = rootDir + "/rawImages/alignImages/" + file.split(".")[0] + "_2d_registered.npy"
data = np.load(fileName)

threshold = detect_threshold(data, nsigma=2.0)

sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(data, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
threshold = bkg.background + (1.0 * bkg.background_rms)  # background-only error image

sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)

segm_deblend = deblend_sources(data, segm, npixels=50, filter_kernel=kernel, nlevels=32, contrast=0.001)

norm = ImageNormalize(stretch=SqrtStretch())

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(data, origin="lower", cmap="Greys_r", norm=norm)
ax1.set_title("Data")
cmap = segm.make_cmap(random_state=12345)
# ax2.imshow(segm, origin='lower', cmap=cmap)
ax2.imshow(segm_deblend, origin="lower", cmap=cmap)
ax2.set_title("Segmentation Image")


#%% will 3D fit sources from original image

dictShifts = loadJSON(rootDir + os.sep + "alignImages.bed.json")

fileNameTif = rootDir + "/rawImages/" + file
I_3D = Image()
I_3D.loadImage(fileNameTif)
ROI = file.split("_")[3]
barcode = file.split("_")[2]
shiftArray = dictShifts["ROI:" + ROI][barcode]
RTmatrix = []

for region in measure.regionprops(label_image, im):
    # take regions with large enough areas
    if region.area >= 25 and region.area < 500:
        # draw rectangle around segmented coins
        minr, minc, maxr, maxc = region.bbox
        x, y = region.weighted_centroid
        subImage3D = I_3D.data[:, minr:maxr, minc:maxc]

        RTmatrix.append(
            [
                ROI,
                barcode,
                x,
                y,
                region.equivalent_diameter,
                region.max_intensity,
                region.mean_intensity,
                region.convex_area,
            ]
        )


del I_3D
