#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:41:58 2021

@author: marcnol

Purpose: Corrects drift in 3D

Often there is drift in the z-position from cycle to cycle.

The drift correction routines take care of the corrections in XY but not in Z.

steps:
    - iterate over ROIs
    - load 3D file for cycle <i>
    - substract background
    - load global alignment for this cycle
    - re-align 3D image using XY alignment
    - segment objects to get labeled masks
    - Get weighted moments and gaussian fits
    
    - display results:
        - overlap XY, XZ, YZ image projections with moment and gaussian fits
    - output results in a Table() using same formatting as that used in segmentMasks.py

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os
import matplotlib.pylab as plt
import numpy as np
from datetime import datetime
from scipy.ndimage import shift as shiftImage
from scipy.optimize import curve_fit
from shutil import copyfile
from skimage import io
import uuid

from tqdm import trange, tqdm
from astropy.visualization import simple_norm
from astropy.table import Table, Column
from photutils import CircularAperture
from dask.distributed import Client, LocalCluster, get_client, as_completed

from numba import jit

from imageProcessing.imageProcessing import Image
from imageProcessing.imageProcessing import (
    _reinterpolatesFocalPlane,
    imageShowWithValues,
    imageShowWithValuesSingle,
    imageAdjust,
    _removesInhomogeneousBackground,
    appliesXYshift3Dimages,
    imageBlockAlignment3D,
    plots3DshiftMatrices,
    combinesBlocksImageByReprojection,
    plots4images,
    makesShiftMatrixHiRes,
    preProcess3DImage,
    _segments3DvolumesByThresholding,
    display3D,
)
from fileProcessing.fileManagement import folders, writeString2File, loadJSON
from fileProcessing.fileManagement import daskCluster, RT2fileName, loadsAlignmentDictionary

from photutils import Background2D, MedianBackground
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats

from astropy.stats import SigmaClip
from imageProcessing.segmentMasks import _showsImageSources
from tqdm import trange, tqdm
from skimage import exposure
from skimage.registration import phase_cross_correlation

from bigfish.detection.spot_modeling import fit_subpixel
from skimage.measure import regionprops


# =============================================================================
# CLASSES
# =============================================================================

class segmentSources3D:
    def __init__(self, param, log1, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.log1 = log1
        self.window = 3
        self.parallel = parallel

    def createsOutputTable(self):
        output = Table(
                names=(
                    "Buid",
                    "ROI #",
                    "CellID #",
                    "Barcode #",
                    "id",
                    "zcentroid",
                    "xcentroid",
                    "ycentroid",                    
                    "sharpness",
                    "roundness1",
                    "roundness2",
                    "npix",
                    "sky",
                    "peak",
                    "flux",
                    "mag",
                ),
                dtype=("S2", 
                       "int", 
                       "int", 
                       "int", 
                       "int",
                       "f4", 
                       "f4", 
                       "f4", 
                       "f4", 
                       "f4", 
                       "f4", 
                       "int", 
                       "f4",
                       "f4",
                       "f4",
                       "f4",
                       ),
            )
        return output 
    
    def getMaskProperties(self, segmentedImage3D, image3D_aligned, threshold=10):
        """
        get object properties from labeled image and formats
        centroids in NPY array

        Parameters
        ----------
        segmentedImage3D : NPY 3D array
            labeled 3D image.
        image3D_aligned : NPY 3D array
            pre-processed 3D image.

        Returns
        -------
        spots : NPY int64 array
            list of spots with the format: zyx

        """
        properties = regionprops(segmentedImage3D, intensity_image=image3D_aligned)
        centroids=[x.weighted_centroid for x in properties]
        sharpness=[x.weighted_centroid for x in properties]
        roundness1=[x.eccentricity for x in properties]
        roundness2=[x.solidity for x in properties]
        npix=[x.area for x in properties]
        sky=[0.0 for x in properties]
        peak=[x.max_intensity for x in properties]
        flux=[x.max_intensity/threshold for x in properties] # peak intensity over the detection threshold
        mag=[-2.5*np.log10(x) for x in flux] # -2.5 log10(flux)

        z=[x[0] for x in centroids]
        y=[x[1] for x in centroids]
        x=[x[2] for x in centroids]
        spots = np.zeros((len(z),3))
        spots[:,0]=z
        spots[:,1]=y
        spots[:,2]=x
        spots=spots.astype('int64')

        return (
                spots,
                sharpness,
                roundness1,
                roundness2,
                npix,
                sky,
                peak,
                flux,
                mag,
                ) 
    
    def segmentSources3DinFolder(self):
        """
        Fits sources in all files in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()

        referenceBarcode = self.param.param["alignImages"]["referenceFiducial"]
        self.log1.info("\nReference Barcode: {}".format(referenceBarcode))
        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)

        threshold_over_std,sigma ,boxSize, filter_size=1, 3, (32, 32),(3, 3)
        area_min, area_max, nlevels, contrast = 3, 1000, 64, 0.001
        voxel_size_z, voxel_size_yx,psf_z, psf_yx = 250, 100, 500, 200 

        fileNameReferenceList, ROIList = RT2fileName(self.param, referenceBarcode)

        numberROIs = len(ROIList)
        self.log1.info("\nDetected {} ROIs".format(numberROIs))

        # loads dicShifts with shifts for all ROIs and all labels
        
        dictShifts, dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder, self.log1)

        # creates Table that will hold results
        output = self.createsOutputTable()

        if numberROIs > 0:

            # loops over ROIs
            for fileNameReference in fileNameReferenceList:
                # loads reference fiducial image for this ROI
                ROI = ROIList[fileNameReference]
                fileName2ProcessList = [x for x in self.param.fileList2Process\
                                        if self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI]

                print("Found {} files in ROI: {}".format(len(fileName2ProcessList), ROI))
                print("[roi:cycle] {}".format("|".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                                + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"]) for x in fileName2ProcessList])))

                for fileName2Process in self.param.fileList2Process:
                    # excludes the reference fiducial and processes files in the same ROI
                    roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
                    label = os.path.basename(fileName2Process).split("_")[2]  # to FIX

                    if roi == ROI:

                        # - load  and preprocesses 3D fiducial file
                        print("\n\nProcessing cycle {}".format(os.path.basename(fileName2Process)))
                        image3D0 = io.imread(fileName2Process).squeeze()
                        image3D = preProcess3DImage(image3D0, self.lower_threshold, self.higher_threshold)

                        # drifts 3D stack in XY
                        if dictShiftsAvailable:
                            # uses existing shift calculated by alignImages
                            try:
                                shift = dictShifts["ROI:" + roi][label]
                                print("Applying existing XY shift...")
                            except KeyError:
                                shift = None
                                raise SystemExit(
                                    "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                                        "ROI:" + ROI, label))
                                
                        # applies XY shift to 3D stack
                        print("shifts XY = {}".format(shift))
                        image3D_aligned = appliesXYshift3Dimages(image3D, shift)
                        
                        # segments 3D volumes
                        binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D_aligned,
                                                                                    threshold_over_std=threshold_over_std,
                                                                                    sigma = sigma,
                                                                                    boxSize=boxSize,
                                                                                    filter_size=filter_size,
                                                                                    area_min = area_min,
                                                                                    area_max=area_max,
                                                                                    nlevels=nlevels,
                                                                                    contrast=contrast,
                                                                                    deblend3D = True)
                        # gets centroids and converts to spot int64 NPY array
                        (
                            spots,
                            sharpness,
                            roundness1,
                            roundness2,
                            npix,
                            sky,
                            peak,
                            flux,
                            mag,
                            ) = self.getMaskProperties(segmentedImage3D, image3D_aligned,threshold = threshold_over_std)

                        # calls bigfish to get 3D sub-pixel coordinates based on 3D gaussian fitting
                        spots_subpixel = fit_subpixel(image3D_aligned, 
                                                      spots, 
                                                      voxel_size_z=voxel_size_z, 
                                                      voxel_size_yx=voxel_size_yx,
                                                      psf_z=psf_z, 
                                                      psf_yx=psf_yx)

                        # updates table
                        for i in range(spots_subpixel.shape[0]):   
                            z,x,y = spots_subpixel[i,:]
                            Table_entry = [str(uuid.uuid4()),
                                           roi,
                                           0,
                                           int(label.split('RT')[1]),
                                           i,
                                           z,
                                           x,
                                           y,
                                           sharpness[i],
                                           roundness1[i],
                                           roundness2[i],
                                           npix[i],
                                           sky[i],
                                           peak[i],
                                           flux[i],
                                           mag[i],
                                           ]
                            output.add_row(Table_entry)

                        # fig3 = plt.figure(constrained_layout=False)
                        # fig3.set_size_inches((20 * 2, 20))
                        # gs = fig3.add_gridspec(2, 2)
                        # ax = [fig3.add_subplot(gs[:, 0]), fig3.add_subplot(gs[0, 1]), fig3.add_subplot(gs[1, 1])]

                        # titles = ["Z-projection", "X-projection", "Y-projection"]

                        # for axis, output, i in zip(ax, outputs, self.axes2Plot):
                        #     axis.imshow(output[0])
                        #     axis.set_title(titles[i])

                        # fig3.tight_layout()

                        # fig4 = plots3DshiftMatrices(SSIM_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        # fig4.suptitle("SSIM block matrices")

                        # fig5 = plots3DshiftMatrices(MSE_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        # fig5.suptitle("mean square root block matrices")

                        # fig6 = plots3DshiftMatrices(NRMSE_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        # fig6.suptitle("normalized root mean square block matrices")

                        # # saves figures
                        # figTitles = ['_bkgSubstracted.png','_shiftMatrices.png',
                        #              '_3Dalignments.png','_SSIMblocks.png',
                        #              '_MSEblocks.png','_NRMSEblocks.png']
                        # outputFileNames = ['/home/marcnol/Documents'+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]
                        # outputFileNames = [self.dataFolder.outputFolders["alignImages"]+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]

                        # figs=[fig1,fig2,fig3,fig4,fig5,fig6]
                        # for fig, file in zip(figs,outputFileNames):
                        #     fig.savefig(file)

                        # dict with shiftMatrix and NRMSEmatrix: https://en.wikipedia.org/wiki/Root-mean-square_deviation
                        # These matrices can be used to apply and assess zxy corrections for any pixel in the 3D image
                        # reference file,aligned file,ROI,label,block_i,block_j,shift_z,shift_x,shift_y,quality_xy,quality_zy,quality_zx


        # saves Table with all shifts
        alignmentResultsTable.write(
            self.dataFolder.outputFiles["alignImages"].split(".")[0] + "_block3D.dat",
            format="ascii.ecsv",
            overwrite=True,
        )

        print("alignFiducials3D procesing time: {}".format(datetime.now() - now))


    def alignFiducials3D(self):
        """
        runs refitting routine in rootFolder

        Returns
        -------
        None.

        """
        sessionName = "alignFiducials3D"

        # processes folders and files
        self.log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        self.log1.report("folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.log1.fileNameMD, "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.outputFileName = self.dataFolder.outputFiles["alignImages"]

        self.log1.report("-------> Processing Folder: {}".format(self.currentFolder))
        self.log1.parallel = self.parallel

        self.blockSizeXY = 128
        self.upsample_factor=100
        self.lower_threshold = 0.9
        self.higher_threshold=0.9999
        self.axes2Plot = range(3)

        self.alignFiducials3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        self.log1.report("HiM matrix in {} processed".format(self.currentFolder), "info")

        return 0

    