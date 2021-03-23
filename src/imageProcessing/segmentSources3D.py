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
import numpy as np
from datetime import datetime
from skimage import io
import uuid

from tqdm import trange, tqdm
from astropy.table import Table
from dask.distributed import Client, get_client

from imageProcessing.imageProcessing import (
    appliesXYshift3Dimages,
    preProcess3DImage,
    _segments3DvolumesByThresholding,
    display3D_assembled,
    _reinterpolatesFocalPlane,
)
from fileProcessing.fileManagement import folders, writeString2File, try_get_client, restart_client
from fileProcessing.fileManagement import RT2fileName, loadsAlignmentDictionary

from skimage import exposure

from bigfish.detection.spot_modeling import fit_subpixel
from skimage.measure import regionprops
from skimage.util.apply_parallel import apply_parallel


# =============================================================================
# CLASSES
# =============================================================================

class segmentSources3D:
    def __init__(self, param, log1, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.log1 = log1
        # self.window = 3
        self.parallel = parallel

        self.p=dict()

        # parameters used for 3D segmentation and deblending
        self.p["threshold_over_std"]=1
        self.p["sigma"]=3
        self.p["boxSize"]=(32, 32)
        self.p["filter_size"]=(3, 3)
        self.p["area_min"]=3
        self.p["area_max"]=1000
        self.p["nlevels"]=64
        self.p["contrast"]=0.001

        # parameters used for 3D gaussian fitting
        self.p["voxel_size_z"] = 250
        self.p["voxel_size_yx"] = 100
        self.p["psf_z"] = 500
        self.p["psf_yx"] = 200

        # range used for adjusting image levels during pre-precessing
        self.p["lower_threshold"] = 0.9
        self.p["higher_threshold"] = 0.9999

        # parameters used for plotting 3D image

        # sets the number of planes around the center of the image used to represent localizations in XZ and ZY
        self.p["windowDisplay"]= 10

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

    def getMaskProperties(self, segmentedImage3D, image3D_aligned, threshold=10,nTolerance=1000):
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
        # gets object properties
        properties = regionprops(segmentedImage3D, intensity_image=image3D_aligned)

        # selects nTolerance brightest spots and keeps only these for further processing
        peak0=[x.max_intensity for x in properties]
        peakList = peak0.copy()
        peakList.sort()
        last2keep=np.min([nTolerance,len(peakList)])
        highestPeakValue  = peakList[-last2keep]
        selection = list(np.nonzero(peak0>highestPeakValue)[0])

        # attributes properties using the brightests spots selected
        peak=[properties[x].max_intensity for x in selection]
        centroids=[properties[x].weighted_centroid for x in selection]
        sharpness=[float(properties[x].filled_area/properties[x].bbox_area) for x in selection]
        roundness1=[properties[x].equivalent_diameter for x in selection]
        roundness2=[properties[x].extent for x in selection]
        npix=[properties[x].area for x in selection]
        sky=[0.0 for x in selection]
        peak=[properties[x].max_intensity for x in selection]
        flux=[properties[x].max_intensity/threshold for x in selection] # peak intensity over the detection threshold
        mag=[-2.5*np.log10(x) for x in flux] # -2.5 log10(flux)

        # centroids=[x.weighted_centroid for x in properties]
        # sharpness=[float(x.filled_area/x.bbox_area) for x in properties]
        # roundness1=[x.equivalent_diameter for x in properties]
        # roundness2=[x.extent for x in properties]
        # npix=[x.area for x in properties]
        # sky=[0.0 for x in properties]
        # peak=[x.max_intensity for x in properties]
        # flux=[100*x.max_intensity/threshold for x in properties] # peak intensity over the detection threshold
        # mag=[-2.5*np.log10(x) for x in flux] # -2.5 log10(flux)

        # converts centroids to spot coordinates for bigfish to run 3D gaussian fits
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
    def plotsImage3D(self,image3D,localizations=None):
        '''
        makes list with XY, XZ and ZY projections and sends for plotting

        Parameters
        ----------
        image3D : numpy array
            image in 3D.
        localizations : list
            list of localizations to overlay onto 3D image. The default is None.

        Returns
        -------
        figure handle

        '''
        img = image3D
        center = int(img.shape[1]/2)
        window = self.p["windowDisplay"]

        images = list()
        images.append(np.sum(img,axis=0))
        images.append(np.sum(img[:,:,center-window:center+window],axis=2))
        images.append(np.sum(img[:,center-window:center+window,:],axis=1))

        fig1 = display3D_assembled(images, localizations = localizations, plottingRange = [center,window])

        return fig1

    def segmentSources3DinFolder(self):
        """
        Fits sources in all files in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()

        referenceBarcode = self.param.param["alignImages"]["referenceFiducial"]
        brightest = self.param.param["segmentedObjects"]["brightest"]
        blockSizeXY = self.param.param["zProject"]["blockSize"]
        zWindow = self.param.param["zProject"]["zwindows"]
        # get blocksize & window from infoList !!!!!!!

        self.log1.info("\nReference Barcode: [{}]".format(referenceBarcode))

        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)

        p=self.p

        fileNameReferenceList, ROIList = RT2fileName(self.param, referenceBarcode)

        numberROIs = len(ROIList)
        self.log1.info("\nDetected {} ROIs".format(numberROIs))

        # loads dicShifts with shifts for all ROIs and all labels
        dictShifts, dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder, self.log1)

        # creates Table that will hold results
        outputTable = self.createsOutputTable()

        # restart_client()

        if numberROIs > 0:

            # loops over ROIs
            for fileNameReference in fileNameReferenceList:
                # loads reference fiducial image for this ROI
                ROI = ROIList[fileNameReference]
                fileName2ProcessList = [x for x in self.param.fileList2Process\
                                        if self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI]
                Nfiles2Process=len(fileName2ProcessList)
                print("Found {} files in ROI [{}]".format(Nfiles2Process, ROI))
                print("[roi:cycle] {}".format(" | ".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                                + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"]) for x in fileName2ProcessList])))

                for fileIndex, fileName2Process in enumerate(self.param.fileList2Process):

                    # excludes the reference fiducial and processes files in the same ROI
                    roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
                    # label = os.path.basename(fileName2Process).split("_")[2]  # to FIX
                    label = str(self.param.decodesFileParts(os.path.basename(fileName2Process))["cycle"])

                    if roi == ROI:

                        # - load  and preprocesses 3D fiducial file
                        print("\n\nProcessing roi:[{}] cycle:[{}] {}/{}".format(roi,label,fileIndex,Nfiles2Process))
                        print("File:{}".format(os.path.basename(fileName2Process)))
                        image3D0 = io.imread(fileName2Process).squeeze()

                        # restricts analysis to a sub volume containing sources
                        focalPlaneMatrix, zRange, _= _reinterpolatesFocalPlane(image3D0,blockSizeXY = blockSizeXY, window=zWindow)
                        # zRange = (40,range(30,50))
                        image3D = image3D0[zRange[1],:,:].copy()
                        zOffset = zRange[1][0]
                        print("Focal plane found: {}, zRange = {}, imageSize = {}".format(zRange[0],zRange[1],image3D.shape))

                        # preprocesses image by background substraction and level normalization
                        image3D = preProcess3DImage(image3D, self.p["lower_threshold"], self.p["higher_threshold"])

                        # drifts 3D stack in XY
                        if dictShiftsAvailable and  label != referenceBarcode:
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
                        if label != referenceBarcode:
                            print("Applies shift = {}".format(shift))
                            image3D_aligned = appliesXYshift3Dimages(image3D, shift)
                        else:
                            print("Running reference fiducial cycle: no shift applied!")
                            shift = np.array([0.,0.])
                            image3D_aligned = image3D

                        # segments 3D volumes
                        # restart_client()
                        binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D_aligned,
                                                                                    threshold_over_std=p["threshold_over_std"],
                                                                                    sigma = p["sigma"],
                                                                                    boxSize=p["boxSize"],
                                                                                    filter_size=p["filter_size"],
                                                                                    area_min = p["area_min"],
                                                                                    area_max=p["area_max"],
                                                                                    nlevels=p["nlevels"],
                                                                                    contrast=p["contrast"],
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
                            ) = self.getMaskProperties(segmentedImage3D, image3D_aligned,threshold = p["threshold_over_std"],nTolerance=brightest)

                        print(">Recovered {} localizations".format(len(peak)))

                        print(">Rescales image values after reinterpolation")
                        image3D_aligned = exposure.rescale_intensity(image3D_aligned, out_range=(0, 1)) # removes negative backgrounds

                        print(">Refits spots using gaussian 3D fittings...")
                        # calls bigfish to get 3D sub-pixel coordinates based on 3D gaussian fitting
                        spots_subpixel = fit_subpixel(image3D_aligned,
                                                      spots,
                                                      voxel_size_z=p["voxel_size_z"],
                                                      voxel_size_yx=p["voxel_size_yx"],
                                                      psf_z=p["psf_z"],
                                                      psf_yx=p["psf_yx"])

                        # updates table
                        for i in range(spots_subpixel.shape[0]):
                            z,x,y = spots_subpixel[i,:]
                            Table_entry = [str(uuid.uuid4()),
                                           roi,
                                           0,
                                           int(label.split('RT')[1]),
                                           i,
                                           z+zOffset,
                                           y,
                                           x,
                                           sharpness[i],
                                           roundness1[i],
                                           roundness2[i],
                                           npix[i],
                                           sky[i],
                                           peak[i],
                                           flux[i],
                                           mag[i],
                                           ]
                            outputTable.add_row(Table_entry)

                        # represents image in 3D with localizations
                        figures=list()
                        figures.append([self.plotsImage3D(image3D_aligned,localizations=[spots_subpixel,spots]),'_3DimageNlocalizations.png'])

                        # saves figures
                        outputFileNames = [self.dataFolder.outputFolders["segmentedObjects"]+os.sep+os.path.basename(fileName2Process)+x[1] for x in figures]

                        for fig, file in zip(figures,outputFileNames):
                            fig[0].savefig(file)

                        # saves Table with all shifts in every iteration to avoid loosing computed data
                        outputTable.write(
                            self.outputFileName,
                            format="ascii.ecsv",
                            overwrite=True,
                        )

                        del image3D_aligned, image3D, image3D0

        print("segmentSources3D procesing time: {}".format(datetime.now() - now))


    def segmentSources3D(self):
        """
        runs 3D fitting routine in rootFolder

        Returns
        -------
        None.

        """
        sessionName = "segmentSources3D"

        # processes folders and files
        self.log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        self.log1.report("folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.log1.fileNameMD, "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.label = self.param.param["acquisition"]["label"]
        self.outputFileName = self.dataFolder.outputFiles["segmentedObjects"] + "_3D_" + self.label + ".dat"

        self.log1.report("-------> Processing Folder: {}".format(self.currentFolder))
        self.log1.parallel = self.parallel

        self.segmentSources3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        self.log1.report("segmentedObjects run in {} finished".format(self.currentFolder), "info")

        return 0

