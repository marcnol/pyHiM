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

from astropy.table import Table,vstack

from imageProcessing.imageProcessing import (
    appliesXYshift3Dimages,
    preProcess3DImage,
    _segments3DvolumesByThresholding,
    display3D_assembled,
    _reinterpolatesFocalPlane,
    reinterpolateZ,
    _plotsImage3D,
    _segments3Dvolumes_StarDist,
    imageAdjust,
)
from fileProcessing.fileManagement import folders, writeString2File, printDict,try_get_client, getDictionaryValue
from fileProcessing.fileManagement import loadsAlignmentDictionary, retrieveNumberROIsFolder, printLog

from skimage import exposure

from apifish.detection.spot_modeling import fit_subpixel
from skimage.measure import regionprops


# =============================================================================
# CLASSES
# =============================================================================

class segmentSources3D:
    def __init__(self, param, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.parallel = parallel

        self.p=dict()

        # parameters from infoList.json
        self.p["referenceBarcode"] = self.param.param["alignImages"]["referenceFiducial"]
        self.p["brightest"] = self.param.param["segmentedObjects"]["brightest"]
        self.p["blockSizeXY"] = self.param.param["zProject"]["blockSize"]
        self.p["regExp"] =self.param.param["acquisition"]["fileNameRegExp"]
        self.p["zBinning"] = getDictionaryValue(self.param.param['acquisition'], "zBinning", default=1)

        self.p["zWindow"] = int(self.param.param["zProject"]["zwindows"]/self.p["zBinning"])
        self.p["pixelSizeXY"] = self.param.param["acquisition"]["pixelSizeXY"]
        self.p["pixelSizeZ"] = self.param.param["acquisition"]["pixelSizeZ"]

        self.p["parallelizePlanes"] = getDictionaryValue(self.param.param['acquisition'], "parallelizePlanes", default=1)


        # decides what segmentation method to use
        self.p["3Dmethod"]=getDictionaryValue(self.param.param["segmentedObjects"], "3Dmethod", default='thresholding')
        self.p["reducePlanes"]=getDictionaryValue(self.param.param["segmentedObjects"], "reducePlanes", default=True)

        # parameters used for 3D segmentation and deblending
        self.p["threshold_over_std"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_threshold_over_std", default=1)
        self.p["sigma"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_sigma", default=3)
        boxSize=getDictionaryValue(self.param.param["segmentedObjects"], "3D_boxSize", default=32)
        self.p["boxSize"] = (boxSize,boxSize)
        filter_size = getDictionaryValue(self.param.param["segmentedObjects"], "3D_filter_size", default=3)
        self.p["filter_size"]=(filter_size,filter_size)
        self.p["area_min"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_area_min", default=3)
        self.p["area_max"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_area_max", default=1000)
        self.p["nlevels"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_nlevels", default=64)
        self.p["contrast"]=getDictionaryValue(self.param.param["segmentedObjects"], "3D_contrast", default=0.001)

        # parameters for stardist
        self.p["stardist_basename"]=getDictionaryValue(self.param.param["segmentedObjects"], "stardist_basename", default='/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/')
        self.p["stardist_network"]=getDictionaryValue(self.param.param["segmentedObjects"], "stardist_network3D", default='stardist_18032021_single_loci')

        # parameters used for 3D gaussian fitting
        self.p["voxel_size_z"] = float(1000*self.p["pixelSizeZ"]*self.p["zBinning"])
        self.p["voxel_size_yx"] = float(1000*self.p["pixelSizeXY"])
        self.p["psf_z"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_psf_z", default=500)
        self.p["psf_yx"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_psf_yx", default=200)

        # range used for adjusting image levels during pre-precessing
        self.p["lower_threshold"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_lower_threshold", default=0.9)
        self.p["higher_threshold"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_higher_threshold", default=0.9999)

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

        if len(properties)>0:
            # selects nTolerance brightest spots and keeps only these for further processing
            try:
                   # compatibility with scikit_image versions <= 0.18
                   peak0=[x.max_intensity for x in properties]
            except AttributeError:
                   # compatibility with scikit_image versions >=0.19
                   peak0=[x.intensity_max for x in properties]

            peakList = peak0.copy()
            peakList.sort()
            
            if nTolerance == "None":
                last2keep = len(peakList)
            else:
                last2keep = np.min([nTolerance,len(peakList)])
                
            highestPeakValue  = peakList[-last2keep]
            selection = list(np.nonzero(peak0>highestPeakValue)[0])

            # attributes properties using the brightests spots selected
            try:
                   # compatibility with scikit_image versions <= 0.18
                   peak=[properties[x].max_intensity for x in selection]
                   centroids=[properties[x].weighted_centroid for x in selection]
                   sharpness=[float(properties[x].filled_area/properties[x].bbox_area) for x in selection]
                   roundness1=[properties[x].equivalent_diameter for x in selection]
            except AttributeError:
                   # compatibility with scikit_image versions >=0.19
                   peak=[properties[x].intensity_max for x in selection]
                   centroids=[properties[x].centroid_weighted for x in selection]
                   sharpness=[float(properties[x].area_filled/properties[x].area_bbox) for x in selection]
                   roundness1=[properties[x].equivalent_diameter_area for x in selection]

            roundness2=[properties[x].extent for x in selection]
            npix=[properties[x].area for x in selection]
            sky=[0.0 for x in selection]

            try:
                   # compatibility with scikit_image versions <= 0.18
                   peak=[properties[x].max_intensity for x in selection]
                   flux=[100*properties[x].max_intensity/threshold for x in selection] # peak intensity over t$
            except AttributeError:
                   # compatibility with scikit_image versions >=0.19
                   peak=[properties[x].intensity_max for x in selection]
                   flux=[100*properties[x].intensity_max/threshold for x in selection] # peak intensity$

            mag=[-2.5*np.log10(x) for x in flux] # -2.5 log10(flux)

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
        else:
            # creates output lists to return
            return [],[],[],[],[],[],[],[],[]

    def plotsImage3D(self,image3D,localizations=None,masks=None,normalize=None):
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
        window = self.p["windowDisplay"]
        fig1 = _plotsImage3D(image3D,localizations=localizations,masks=masks,normalize=normalize,window = window)
        return fig1

    def _segments3Dvolumes(self,image3D_aligned):
        p=self.p

        if 'stardist' in p["3Dmethod"]:
            binary, segmentedImage3D  = _segments3Dvolumes_StarDist(image3D_aligned,
                                        area_min = p["area_min"],
                                        area_max=p["area_max"],
                                        nlevels=p["nlevels"],
                                        contrast=p["contrast"],
                                        deblend3D = True,
                                        axis_norm=(0,1,2),
                                        model_dir=p["stardist_basename"],
                                        model_name=p["stardist_network"])
        else:
            binary, segmentedImage3D = _segments3DvolumesByThresholding(image3D_aligned,
                                        threshold_over_std=p["threshold_over_std"],
                                        sigma = p["sigma"],
                                        boxSize=p["boxSize"],
                                        filter_size=p["filter_size"],
                                        area_min = p["area_min"],
                                        area_max=p["area_max"],
                                        nlevels=p["nlevels"],
                                        contrast=p["contrast"],
                                        deblend3D = True,
                                        parallelExecution=self.innerParallelLoop)

        return binary, segmentedImage3D

    def segmentSources3D_file(self,fileName2Process):

        p = self.p
        # excludes the reference fiducial and processes files in the same ROI
        roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
        label = str(self.param.decodesFileParts(os.path.basename(fileName2Process))["cycle"])

        # creates Table that will hold results
        outputTable = self.createsOutputTable()

        # - load  and preprocesses 3D fiducial file
        printLog("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi,label))
        printLog("$ File:{}".format(os.path.basename(fileName2Process)))
        image3D0 = io.imread(fileName2Process).squeeze()

        # reinterpolates image in z if necessary
        image3D0 = reinterpolateZ(image3D0, range(0,image3D0.shape[0],p["zBinning"]),mode='remove')

        # restricts analysis to a sub volume containing sources
        if p["reducePlanes"]:
            focalPlaneMatrix, zRange, _= _reinterpolatesFocalPlane(image3D0,blockSizeXY = p["blockSizeXY"], window=p["zWindow"])
            zOffset = zRange[1][0]
            image3D = image3D0[zRange[1],:,:].copy()
            printLog("$ Focal plane found: {}, zRange = {}, imageSize = {}".format(zRange[0],zRange[1],image3D.shape))
        else:
            image3D = image3D0.copy()
            zOffset = 0
            printLog("$ zRange used = 0-{}".format(image3D.shape[0]))

        # preprocesses image by background substraction and level normalization
        if 'stardist' not in p["3Dmethod"]:
            image3D = preProcess3DImage(image3D,
                                        p["lower_threshold"],
                                        p["higher_threshold"],
                                        parallelExecution=self.innerParallelLoop)

        # drifts 3D stack in XY
        shift = None
        if self.dictShiftsAvailable and  label != p["referenceBarcode"]:
            # uses existing shift calculated by alignImages
            try:
                shift = self.dictShifts["ROI:" + roi][label]
                printLog("> Applying existing XY shift...")
            except KeyError:
                shift = None

        if shift is None and label != p["referenceBarcode"]:
            raise SystemExit(
                "> Existing with ERROR: Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                    "ROI:" + self.ROI, label))

        # applies XY shift to 3D stack
        if label != p["referenceBarcode"]:
            # printLog("$ Applies shift = {:.2f}".format(shift))
            printLog("$ Applies shift = [{:.2f} ,{:.2f}]".format(shift[0],shift[1]))
            image3D_aligned = appliesXYshift3Dimages(image3D, shift,parallelExecution=self.innerParallelLoop)
        else:
            printLog("$ Running reference fiducial cycle: no shift applied!")
            shift = np.array([0.,0.])
            image3D_aligned = image3D

        # segments 3D volumes
        binary, segmentedImage3D  = self._segments3Dvolumes(image3D_aligned)

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
            ) = self.getMaskProperties(segmentedImage3D, image3D_aligned,threshold = p["threshold_over_std"],nTolerance=p["brightest"])

        numberSources = len(peak)
        printLog("$ Number of sources detected by image segmentation: {}".format(numberSources))

        if numberSources >0:
            printLog("> Refits spots using gaussian 3D fittings...")

            printLog(" > Rescales image values after reinterpolation")
            image3D_aligned = exposure.rescale_intensity(image3D_aligned, out_range=(0, 1)) # removes negative intensity levels

            # calls bigfish to get 3D sub-pixel coordinates based on 3D gaussian fitting
            # compatibility with latest version of bigfish. To be removed if stable.
            try:
                # version 0.4 commit fa0df4f
                spots_subpixel = fit_subpixel(image3D_aligned,
                                              spots,
                                              voxel_size_z=p["voxel_size_z"],
                                              voxel_size_yx=p["voxel_size_yx"],
                                              psf_z=p["psf_z"],
                                              psf_yx=p["psf_yx"])
            except TypeError:                
                # version > 0.5
                spots_subpixel = fit_subpixel(image3D_aligned,
                                              spots,
                                              (p["voxel_size_z"],p["voxel_size_yx"],p["voxel_size_yx"]), # voxel size
                                              (p["psf_z"],p["psf_yx"],p["psf_yx"])) # spot radius

            printLog(" > Updating table and saving results")
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
            if 'stardist' in p["3Dmethod"]:
                image3D_aligned= imageAdjust(image3D_aligned, lower_threshold=p["lower_threshold"], higher_threshold=p["higher_threshold"])[0]

            figures.append([self.plotsImage3D(image3D_aligned,
                                              masks=segmentedImage3D,
                                              localizations=[spots_subpixel,spots]),'_3DimageNlocalizations.png'])

            # saves figures
            outputFileNames = [self.dataFolder.outputFolders["segmentedObjects"]+os.sep+os.path.basename(fileName2Process)+x[1] for x in figures]

            for fig, file in zip(figures,outputFileNames):
                fig[0].savefig(file)

        del image3D_aligned, image3D, image3D0

        return outputTable

    def segmentSources3DinFolder(self):
        """
        Fits sources in all files in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()

        # Reads list of parameters assigned upon Class initialization
        p=self.p
        printDict(p)

        # Finds images to process
        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)
        self.ROIList = retrieveNumberROIsFolder(self.currentFolder, p["regExp"], ext="tif")
        self.numberROIs = len(self.ROIList)
        printLog("\n$ Detected {} ROIs".format(self.numberROIs))
        printLog("$ Number of images to be processed: {}".format(len(self.param.fileList2Process)))

        # loads dicShifts with shifts for all ROIs and all labels
        self.dictShifts, self.dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder)

        # creates Table that will hold results
        outputTableGlobal = self.createsOutputTable()
        outputTables = list()

        if self.p["parallelizePlanes"]:
            client = None
        else:
            client = try_get_client()

        if self.numberROIs > 0:

            # loops over ROIs
            for ROI in self.ROIList:
                # loads reference fiducial image for this ROI
                self.ROI = ROI
                
                self.fileName2ProcessList = [x for x in self.param.fileList2Process\
                                        if self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI and \
                                            "RT" in self.param.decodesFileParts(os.path.basename(x))["cycle"]]

                Nfiles2Process=len(self.fileName2ProcessList)
                printLog("$ Found {} files in ROI [{}]".format(Nfiles2Process, ROI))
                printLog("$ [roi:cycle] {}".format(" | ".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                                + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"]) for x in self.fileName2ProcessList])))

                if client is None:
                    self.innerParallelLoop = True

                    for fileIndex, fileName2Process in enumerate(self.fileName2ProcessList): #self.param.fileList2Process):
                        printLog("\n\n>>>Iteration: {}/{}<<<".format(fileIndex,Nfiles2Process))

                        outputTables.append(self.segmentSources3D_file(fileName2Process))
                else:
                    self.innerParallelLoop = False
                    printLog("> Aligning {} files using {} workers...".format(Nfiles2Process,len(client.scheduler_info()['workers'])))

                    futures = [client.submit(self.segmentSources3D_file,
                                                     x) for x in self.fileName2ProcessList]

                    outputTables = client.gather(futures)
                    printLog(" > Retrieving {} results from cluster".format(len(outputTables)))

                # Merges Tables for different cycles and appends results Table to that of previous ROI
                outputTableGlobal = vstack([outputTableGlobal]+outputTables)

        printLog("$ segmentSources3D procesing time: {}".format(datetime.now() - now))

        # saves Table with all shifts in every iteration to avoid loosing computed data
        outputTableGlobal.write(
            self.outputFileName,
            format="ascii.ecsv",
            overwrite=True,
        )


    def segmentSources3D(self):
        """
        runs 3D fitting routine in rootFolder

        Returns
        -------
        None.

        """
        sessionName = "segmentSources3D"

        # processes folders and files
        printLog("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        printLog("$ folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.param.param["fileNameMD"], "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.label = self.param.param["acquisition"]["label"]
        self.outputFileName = self.dataFolder.outputFiles["segmentedObjects"] + "_3D_" + self.label + ".dat"

        printLog("> Processing Folder: {}".format(self.currentFolder))

        self.segmentSources3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        printLog("$ segmentedObjects run in {} finished".format(self.currentFolder))

        return 0

##########################################
# FUNCTIONS
##########################################
