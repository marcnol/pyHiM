#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 09:17:30 2020

@author: marcnol
"""

import os

import numpy as np

import matplotlib.pyplot as plt

from photutils import Background2D, MedianBackground
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats

from imageProcessing.imageProcessing import Image
from imageProcessing.refitBarcodes3D import refitBarcodesClass, refits_3D_ASTROPY
from fileProcessing.fileManagement import Parameters, log, session, folders
from astropy.table import Table, vstack, Column
from imageProcessing.segmentMasks import _segmentSourceInhomogBackground
from astropy.stats import SigmaClip
from imageProcessing.segmentMasks import _showsImageSources
from tqdm import trange, tqdm

#################
# sets paths
#################
sample = "droso"
# sample='example'
reloadsBarcodeFile = False

if os.path.exists("/home/marcnol/data"):
    # running in Atlantis
    rootFolder = "/home/marcnol/data/Embryo_debug_dataset/bigfish"
    print("Running in Atlantis")
else:
    # running in the lab
    rootFolder = "/mnt/grey/DATA/users/marcnol/test_HiM/bigfish"
    print("Running in the lab")

path_input = rootFolder + sample + "/input"
path_output = rootFolder + sample + "/output"


#################
# loads data
#################
parameterFile = "infoList_barcode.json"
param = Parameters(rootFolder, parameterFile)
log1 = log(rootFolder=rootFolder)
session1 = session(rootFolder, "refit3D")

fittingSession = refitBarcodesClass(param, log1, session1)

if not reloadsBarcodeFile:
    # this will load the image and recalculated the XY source positions
    #################
    # loads raw image
    #################
    Im3D = Image(param, log1)
    Im3D.loadImage(path_input + os.sep + "experiment_1_smfish_fov_1.tif")
    Im3D.maxIntensityProjection()
    image2D = Im3D.data_2D
    image3D = Im3D.data
    numberZplanes = Im3D.data.shape[0]

    # generates sources file
    threshold_over_std = 1  # param.param["segmentedObjects"]["threshold_over_std"]
    fwhm = 3  # param.param["segmentedObjects"]["fwhm"]
    brightest = 2500  # param.param["segmentedObjects"]["brightest"]  # keeps brightest sources

    sigma_clip = SigmaClip(sigma=3)  # param.param["segmentedObjects"]["background_sigma"])

    sources, im1_bkg_substracted = _segmentSourceInhomogBackground(
        image2D, threshold_over_std, fwhm, brightest, sigma_clip
    )
    barcodeMapSinglebarcode = sources.copy()
    barcodeID = "1"
    ROI = 1
    colBarcode = Column(int(barcodeID) * np.ones(len(barcodeMapSinglebarcode)), name="Barcode #", dtype=int)
    colROI = Column(int(ROI) * np.ones(len(barcodeMapSinglebarcode)), name="ROI #", dtype=int)

    barcodeMapSinglebarcode.add_column(colBarcode, index=0)
    barcodeMapSinglebarcode.add_column(colROI, index=0)

else:
    # this will load the image and load the XY source positions already calculated
    fittingSession.dataFolder = folders("/home/marcnol/data/Embryo_debug_dataset/Experiment_18")
    # fittingSession.dataFolder.setsFolders()
    fittingSession.dataFolder.outputFiles = dict()

    fittingSession.dataFolder.outputFiles[
        "segmentedObjects"
    ] = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18/segmentedObjects/segmentedObjects"
    fittingSession.dataFolder.outputFiles[
        "dictShifts"
    ] = "/home/marcnol/data/Embryo_debug_dataset/bigfish/droso/input/alignImages/alignImages.json"

    fileNameBarcodeCoordinates = fittingSession.dataFolder.outputFiles["segmentedObjects"] + "_barcode.dat"

    barcodeMap, _ = fittingSession.loadsBarcodeMap()

    barcodeMapROI = barcodeMap.group_by("ROI #")
    numberROIs = len(barcodeMapROI.groups.keys)
    iROI = 0
    barcodeMapSingleROI = barcodeMap.group_by("ROI #").groups[iROI]

    barcodeMapROI_barcodeID = barcodeMapSingleROI.group_by("Barcode #")
    numberBarcodes = len(barcodeMapROI_barcodeID.groups.keys)

    iBarcode = 1

    barcodeMapSinglebarcode = barcodeMapROI_barcodeID.group_by("Barcode #").groups[iBarcode]

    #################
    # loads raw image
    #################
    Im3D = fittingSession.loadsShifts3Dimage(barcodeMapSinglebarcode)
    image2D = Im3D.data_2D
    image3D = Im3D.data

# # refits 3D positions using gaussian z-profile fitting
# fittingSession.outputFileName = path_output+os.sep+"outputBarcodes"
# barcodeMapSinglebarcode = fittingSession.fitsZpositions(Im3D, barcodeMapSinglebarcode)
# fittingSession.shows3DfittingResults(barcodeMapSinglebarcode, numberZplanes=numberZplanes)

#%%
##################################
# refits z positions using ASTROPY
##################################
# sets parameters
distTolerance = 1  # in px
flux_min = 3
window = 5
bkg_estimator = MedianBackground()
addSources = False

threshold_over_std = 1  # param.param["segmentedObjects"]["threshold_over_std"]
fwhm = 3  # param.param["segmentedObjects"]["fwhm"]
brightest = 100  # param.param["segmentedObjects"]["brightest"]  # keeps brightest sources
sigma_clip = SigmaClip(sigma=3)  # param.param["segmentedObjects"]["background_sigma"])

barcodeMapNew, stats = refits_3D_ASTROPY(
    Im3D,
    barcodeMapSinglebarcode,
    distTolerance=distTolerance,
    flux_min=flux_min,
    window=window,
    bkg_estimator=bkg_estimator,
    addSources=addSources,
    threshold_over_std=threshold_over_std,
    fwhm=fwhm,
    brightest=brightest,
    sigma_clip=sigma_clip,
    path_output=path_output,
)

#%%


# z-gaussian fit
# selection_3Dgaussian = np.nonzero((barcodeMapSinglebarcode["xcentroid"]>xPlane-window) & (barcodeMapSinglebarcode["xcentroid"]<xPlane+window))
# y_3Dgaussian=barcodeMapSinglebarcode["ycentroid"][selection_3Dgaussian ]
# z_3Dgaussian=barcodeMapSinglebarcode["zcentroidGauss"][selection_3Dgaussian ]
# flux_3Dgaussian = barcodeMapSinglebarcode["flux"][selection_3Dgaussian ]

# fig=_showsImageSources(image3D_ZY , image3D_ZY , y_3Dgaussian, z_3Dgaussian, flux_3Dgaussian, percent=99.5,vmin=0,vmax=flux_3Dgaussian.max())
# fig.savefig(path_output + "/_segmentedSourcesAstropyXY_zgaussian.png")


# # initializes results image
# nZplanes=image3D.shape[0]
# resultsImage = np.zeros(((nZplanes)*len(range(window,image3D.shape[2],2*window)),image3D.shape[1]))
# iPlane=0

# for xPlane in trange(window,image3D.shape[2],2*window):

#     # makes 2D YZ image by projection
#     resultsImage[iPlane:iPlane+nZplanes,:] = np.sum(image3D[:,:,xPlane-window:xPlane+window], axis=2)
#     iPlane+=nZplanes

# fig, ax = plt.subplots()
# fig.set_size_inches((40, 200))

# plt.imshow(resultsImage,cmap='Greys')
# fig.savefig(path_output + "/_segmentedSourcesAstropyXYZ_combined.png")
