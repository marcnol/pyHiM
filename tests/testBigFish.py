#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:35:52 2020

@author: marcnol

Download latest development branch as .zip from https://github.com/fish-quant/big-fish/tree/develop
copy into /home/marcnol/Repositories/bigfish
ln -s /home/marcnol/Repositories/bigfish/bigfish /home/marcnol/Repositories/pyHiM/bigfish

"""

import os
import bigfish
from bigfish import stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
import numpy as np

print("Big-FISH version: {0}".format(bigfish.__version__))
import matplotlib.pyplot as plt


#%%

sample = "droso"
sample = "example"

if os.path.exists("/home/marcnol/data"):
    # running in Atlantis
    rootFolder = "/home/marcnol/data/Embryo_debug_dataset/bigfish/"
    print("Running in Atlantis")
else:
    # running in the lab
    rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/bigfish/"
    print("Running in the lab")

path_input = rootFolder + sample + "/input"
path_output = rootFolder + sample + "/output"

if "droso" not in sample:
    stack.check_input_data(path_input)

#%%

recipe = {"fov": "fov_1", "c": ["dapi", "smfish"], "opt": "experiment_1", "ext": "tif", "pattern": "opt_c_fov.ext"}

imageStack = stack.build_stack(recipe, input_folder=path_input)
print("\r shape: {0}".format(imageStack.shape))
print("\r dtype: {0}".format(imageStack.dtype))

image3D = imageStack[0, 1, ...]
print("smfish channel")
print("\r shape: {0}".format(image3D.shape))
print("\r dtype: {0}".format(image3D.dtype))

image2D = stack.maximum_projection(image3D)
print("smfish channel (2D maximum projection)")
print("\r shape: {0}".format(image2D.shape))
print("\r dtype: {0}".format(image2D.dtype))


#%%
# parameters
voxel_size_z = 350
voxel_size_yx = 105
psf_z = 350
psf_yx = 100

# sigma
sigma_z, sigma_yx, sigma_yx = stack.get_sigma(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
print("standard deviation of the PSF (z axis): {:0.3f} pixels".format(sigma_z))
print("standard deviation of the PSF (yx axis): {:0.3f} pixels".format(sigma_yx))


spots, threshold = detection.detect_spots(
    image3D, return_threshold=True, voxel_size_z=voxel_size_z, voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx
)
print("detected spots")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("\r threshold: {0}".format(threshold))

# sigma
sigma_z, sigma_yx, sigma_yx = stack.get_sigma(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
sigma = (sigma_z, sigma_yx, sigma_yx)

# LoG filter
rna_log = stack.log_filter(image3D, sigma)

# local maximum detection
mask = detection.local_maximum_detection(rna_log, min_distance=sigma)

# thresholding
threshold = detection.automated_threshold_setting(rna_log, mask)
spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
print("detected spots")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("\r threshold: {0}".format(threshold))

(radius_z, radius_yx, radius_yx) = stack.get_radius(voxel_size_z, voxel_size_yx, psf_z, psf_yx)
print("radius z axis: {0:0.3f}".format(radius_z))
print("radius yx axis: {0:0.3f}".format(radius_yx))

plot.plot_detection(image2D, spots, radius=radius_yx, framesize=(10, 8), contrast=True)

#%% Dense regions decomposition

spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
    image3D,
    spots,
    voxel_size_z,
    voxel_size_yx,
    psf_z,
    psf_yx,
    alpha=0.7,  # alpha impacts the number of spots per candidate region
    beta=1,  # beta impacts the number of candidate regions to decompose
    gamma=5,
)  # gamma the filtering step to denoise the image

print("detected spots before decomposition")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("detected spots after decomposition")
print("\r shape: {0}".format(spots_post_decomposition.shape))
print("\r dtype: {0}".format(spots_post_decomposition.dtype))


plot.plot_detection(image2D, spots_post_decomposition, radius=radius_yx, framesize=(10, 8), contrast=True)

plot.plot_reference_spot(reference_spot, rescale=True)

#%%  Clusters detection¶

radius = 350
nb_min_spots = 4
spots_post_clustering, clusters = detection.detect_clusters(
    spots_post_decomposition, voxel_size_z, voxel_size_yx, radius, nb_min_spots
)
print("detected spots after clustering")
print("\r shape: {0}".format(spots_post_clustering.shape))
print("\r dtype: {0}".format(spots_post_clustering.dtype))
print("detected clusters")
print("\r shape: {0}".format(clusters.shape))
print("\r dtype: {0}".format(clusters.dtype))

# plot
plot.plot_detection(
    image2D,
    spots=[spots_post_decomposition, clusters[:, :3]],
    shape=["circle", "polygon"],
    radius=[radius_yx, radius_yx * 2],
    color=["red", "blue"],
    linewidth=[1, 2],
    fill=[False, True],
    framesize=(10, 8),
    contrast=True,
)

#%% 2D and 3D detection


spots, threshold = detection.detect_spots(
    image3D, return_threshold=True, voxel_size_z=voxel_size_z, voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx
)
print("detected spots")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("\r threshold: {0}".format(threshold))

(radius_yx, radius_yx) = stack.get_radius(None, voxel_size_yx, None, psf_yx)
print("radius yx axis: {0:0.3f}".format(radius_yx))

spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
    image3D,
    spots,
    voxel_size_z=voxel_size_z,
    voxel_size_yx=voxel_size_yx,
    psf_z=psf_z,
    psf_yx=psf_yx,
    alpha=0.7,
    beta=1,
    gamma=5,
)

print("detected spots before decomposition")
print("\r shape: {0}".format(spots.shape))
print("\r dtype: {0}".format(spots.dtype))
print("detected spots after decomposition")
print("\r shape: {0}".format(spots_post_decomposition.shape))
print("\r dtype: {0}".format(spots_post_decomposition.dtype))

radius = 350
nb_min_spots = 4
spots_post_clustering, clusters = detection.detect_clusters(
    spots_post_decomposition, voxel_size_z, voxel_size_yx, radius, nb_min_spots
)
print("detected spots after clustering")
print("\r shape: {0}".format(spots_post_clustering.shape))
print("\r dtype: {0}".format(spots_post_clustering.dtype))
print("detected clusters")
print("\r shape: {0}".format(clusters.shape))
print("\r dtype: {0}".format(clusters.dtype))


# plot
plot.plot_detection(
    image2D,
    spots=[spots_post_decomposition, clusters[:, :2]],
    shape=["circle", "polygon"],
    radius=[radius_yx, radius_yx * 2],
    color=["red", "blue"],
    linewidth=[1, 2],
    fill=[False, True],
    framesize=(10, 8),
    contrast=True,
)

from imageProcessing.segmentMasks import _showsImageSources

x = spots[:, 2]
y = spots[:, 1]
flux = spots[:, 0]

fig = _showsImageSources(image2D, image2D, x, y, flux, percent=99.5)
fig.savefig(path_output + "/_segmentedSources.png")


#%% subpixel fitting

spots_subpixel = detection.fit_subpixel(
    image3D, spots, voxel_size_z=voxel_size_z, voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx
)

x = spots_subpixel[:, 2]
y = spots_subpixel[:, 1]
z = spots_subpixel[:, 0]

fig = _showsImageSources(image2D, image2D, x, y, z, percent=99.5)
fig.savefig(path_output + "/_segmentedSourcesSubpixel.png")

# YZ plane plotting
image3D_log = stack.log_filter(image3D, sigma)

xPlane = int(image3D.shape[2] / 2)
window = 60
windowD = int(window / 1)
image3D_ZY = np.sum(image3D_log[:, :, xPlane - window : xPlane + window], axis=2)

selection = np.nonzero((x > xPlane - windowD) & (x < xPlane + windowD))
y0 = spots[selection, 1]
z0 = spots[selection, 0]

y = spots_subpixel[selection, 1]
z = spots_subpixel[selection, 0]

fig, ax = plt.subplots()
fig.set_size_inches((50, 50))

plt.imshow(image3D_ZY, cmap="nipy_spectral")
plt.scatter(y0, z0, color="r", alpha=0.4, marker="x")
plt.scatter(y, z, color="w", alpha=0.7, marker="+")

fig.savefig(path_output + "/_segmentedSourcesSubpixel_YZplane.png")


#%% localizes in 2d using astropy then gets 3D using bigfish

# 2D fitting
from imageProcessing.segmentMasks import _segmentSourceInhomogBackground
from astropy.stats import SigmaClip

threshold_over_std = 1  # param.param["segmentedObjects"]["threshold_over_std"]
fwhm = 3  # param.param["segmentedObjects"]["fwhm"]
brightest = 1100  # param.param["segmentedObjects"]["brightest"]  # keeps brightest sources

sigma_clip = SigmaClip(sigma=3)  # param.param["segmentedObjects"]["background_sigma"])

sources, im1_bkg_substracted = _segmentSourceInhomogBackground(image2D, threshold_over_std, fwhm, brightest, sigma_clip)

flux = sources["flux"]
x = sources["xcentroid"]
y = sources["ycentroid"]

fig = _showsImageSources(image2D, image2D, x, y, flux, percent=99.5, vmin=0, vmax=flux.max())
fig.savefig(path_output + "/_segmentedSourcesAstropyXY.png")

#%% compares ASTROPY and BIGFISH

x_ASTROPY = sources["xcentroid"]
y_ASTROPY = sources["ycentroid"]
x_BIGFISH = spots_subpixel[:, 2]
y_BIGFISH = spots_subpixel[:, 1]

fig, ax = plt.subplots()
fig.set_size_inches((50, 50))

plt.imshow(np.sum(image3D_log, axis=0), cmap="nipy_spectral")
plt.scatter(x_BIGFISH, y_BIGFISH, color="r", alpha=0.9, marker="+")
plt.scatter(x_ASTROPY, y_ASTROPY, color="w", alpha=0.9, marker="x")

fig.savefig(path_output + "/_segmentedSources_BIGFISH_vs_ASTROPY_XY.png")

#%% 3D fitting using bigfish

spotsAstropy = np.zeros((len(sources), 3))
spotsAstropy[:, 0] = float(image3D.shape[0] / 2)
spotsAstropy[:, 1] = sources["ycentroid"]
spotsAstropy[:, 2] = sources["xcentroid"]
spotsAstropy = spotsAstropy.astype(int)
spotsAstropy_subpixel = detection.fit_subpixel(
    image3D, spotsAstropy, voxel_size_z=voxel_size_z, voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx
)

x = spotsAstropy_subpixel[:, 2]
y = spotsAstropy_subpixel[:, 1]
x0 = sources["xcentroid"]
y0 = sources["ycentroid"]
z = spotsAstropy_subpixel[:, 0]

fig = _showsImageSources(image2D, image2D, x, y, z, percent=99.5, vmin=z.max(), vmax=z.max())
fig.savefig(path_output + "/_segmentedSourcesAstropyXYZ.png")

#%% plots YZ
xPlane = int(image3D.shape[2] / 2)
window = 60
image3D_ZY = np.sum(image3D_log[:, :, xPlane - window : xPlane + window], axis=2)
# image3D_ZY = np.sum(image3D[:,:,xPlane-window:xPlane+window], axis=2)

selection = np.nonzero((spotsAstropy[:, 2] > xPlane - window) & (spotsAstropy[:, 2] < xPlane + window))
y = spotsAstropy_subpixel[selection, 1]
z = spotsAstropy_subpixel[selection, 0]

fig, ax = plt.subplots()
fig.set_size_inches((50, 50))

image3D_ZY = np.sum(image3D_log[:, :, xPlane - window : xPlane + window], axis=2)

plt.imshow(image3D_ZY, cmap="nipy_spectral")
plt.scatter(y, z, color="w", alpha=0.8, marker="x")

fig.savefig(path_output + "/_segmentedSourcesAstropyXYZ_YZplane.png")

#%% refits z positions from XY astropy coordinates
from imageProcessing.imageProcessing import Image
from imageProcessing.refitBarcodes3D import refitBarcodesClass
from fileProcessing.fileManagement import Parameters, log, session
from astropy.table import Table, vstack, Column

parameterFile = "infoList_barcode.json"
param = Parameters(rootFolder, parameterFile)
log1 = log(rootFolder=rootFolder)
session1 = session(rootFolder, "refit3D")

fittingSession = refitBarcodesClass(param, log1, session1)

Im3D = Image(param, log1)

Im3D.loadImage(path_input + os.sep + "experiment_1_smfish_fov_1.tif")
Im3D.maxIntensityProjection()
numberZplanes = Im3D.data.shape[0]

barcodeMapSinglebarcode = sources.copy()
barcodeID = "1"
ROI = 1
colBarcode = Column(int(barcodeID) * np.ones(len(barcodeMapSinglebarcode)), name="Barcode #", dtype=int)
colROI = Column(int(ROI) * np.ones(len(barcodeMapSinglebarcode)), name="ROI #", dtype=int)

barcodeMapSinglebarcode.add_column(colBarcode, index=0)
barcodeMapSinglebarcode.add_column(colROI, index=0)

fittingSession.outputFileName = path_output + os.sep + "outputBarcodes"

barcodeMapSinglebarcode = fittingSession.fitsZpositions(Im3D, barcodeMapSinglebarcode)


fittingSession.shows3DfittingResults(barcodeMapSinglebarcode, numberZplanes=numberZplanes)

#%% refits Z using BIGFISH

# creates matrix from ASTROPY sources with format zyx
spotsAstropy2 = np.zeros((len(sources), 3))
spotsAstropy2[:, 0] = barcodeMapSinglebarcode["zcentroidGauss"]
spotsAstropy2[:, 1] = barcodeMapSinglebarcode["ycentroid"]
spotsAstropy2[:, 2] = barcodeMapSinglebarcode["xcentroid"]

# removes nans from the z-axis
theseRnotNans = ~np.isnan(spotsAstropy2[:, 0])
spotsAstropy2 = spotsAstropy2[theseRnotNans, :]

# converts to int
spotsAstropy2 = spotsAstropy2.astype(int)

# calls BIGFISH to refit in 3D
spotsAstropy2_subpixel = detection.fit_subpixel(
    image3D, spotsAstropy2, voxel_size_z=voxel_size_z, voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx
)

# retrives results in xyz
x2 = spotsAstropy2_subpixel[:, 2]
y2 = spotsAstropy2_subpixel[:, 1]
z2 = spotsAstropy2_subpixel[:, 0]

#%% Makes YZ slice and refits using ASTROPY
from photutils import Background2D, MedianBackground
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats


# sets parameters
flux_min = 5
window = 2

# makes 2D YZ image by projection
xPlane = int(image3D.shape[2] / 2)
image3D_ZY = np.sum(image3D[:, :, xPlane - window : xPlane + window], axis=2)

# finds the yz coordinates from ASTROPY XY fitting and z-gaussian reinterpolation
selection_ASTROPY_3D = np.nonzero(
    (barcodeMapSinglebarcode["xcentroid"] > xPlane - window) & (barcodeMapSinglebarcode["xcentroid"] < xPlane + window)
)
y_ASTROPY_3D = barcodeMapSinglebarcode["ycentroid"][selection_ASTROPY_3D]
z_ASTROPY_3D = barcodeMapSinglebarcode["zcentroidGauss"][selection_ASTROPY_3D]
flux_ASTROPY_3D = barcodeMapSinglebarcode["flux"][selection_ASTROPY_3D]

# finds the yz coordinates from ASTROPY XY fitting and BIGFISH 3D gaussian fitting
selection_3D_bigfish = np.nonzero(
    (spotsAstropy2_subpixel[:, 2] > xPlane - window) & (spotsAstropy2_subpixel[:, 2] < xPlane + window)
)
y_3D_bigfish = spotsAstropy2_subpixel[selection_3D_bigfish, 1]
z_3D_bigfish = spotsAstropy2_subpixel[selection_3D_bigfish, 0]

# re fits in YZ using ASTROPY
bkg_estimator = MedianBackground()
bkg = Background2D(image3D_ZY, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
image3D_ZY_bkg_substracted = image3D_ZY - bkg.background

mean, median, std = sigma_clipped_stats(image3D_ZY_bkg_substracted, sigma=3.0)
daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_over_std * std, brightest=brightest, exclude_border=False,)
sourcesYZ_fitting = daofind(image3D_ZY_bkg_substracted)
flux2 = sourcesYZ_fitting["flux"]
keepHighFlux = np.nonzero(flux2 > flux_min)
(Y_YZ_fitting, Z_YZ_fitting, flux_YZ_fitting,) = (
    sourcesYZ_fitting["xcentroid"][keepHighFlux],
    sourcesYZ_fitting["ycentroid"][keepHighFlux],
    sourcesYZ_fitting["flux"][keepHighFlux],
)

# plots results
fig, ax = plt.subplots()
fig.set_size_inches((50, 50))

# image3D_ZY = np.sum(image3D_log[:,:,xPlane-window:xPlane+window], axis=2)
# image3D_ZY = np.sum(image3D[:,:,xPlane-window:xPlane+window], axis=2)


plt.imshow(image3D_ZY, cmap="nipy_spectral")
plt.scatter(y_ASTROPY_3D, z_ASTROPY_3D, color="w", alpha=0.8, marker="x")
plt.scatter(y_3D_bigfish, z_3D_bigfish, color="y", alpha=0.8, marker="+")
plt.scatter(Y_YZ_fitting, Z_YZ_fitting, color="r", alpha=0.5, marker="o")
fig.savefig(path_output + "/_segmentedSourcesAstropyXYZ_YZplane_nipy.png")

# fig=_showsImageSources(image3D_ZY , image3D_ZY , y, z, flux, percent=99.5,vmin=flux.max(),vmax=flux.max())

fig = _showsImageSources(
    image3D_ZY,
    image3D_ZY,
    Y_YZ_fitting,
    Z_YZ_fitting,
    flux_YZ_fitting,
    percent=99.5,
    vmin=0,
    vmax=flux_YZ_fitting.max(),
)

fig.savefig(path_output + "/_segmentedSourcesAstropyXYZ_YZplane.png")
