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
)
from fileProcessing.fileManagement import folders, writeString2File, printDict,try_get_client
from fileProcessing.fileManagement import loadsAlignmentDictionary, retrieveNumberROIsFolder

from skimage import exposure

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
        # self.window = 3
        self.parallel = parallel

        self.p=dict()

        # parameters from infoList.json
        self.p["referenceBarcode"] = self.param.param["alignImages"]["referenceFiducial"]
        self.p["brightest"] = self.param.param["segmentedObjects"]["brightest"]
        self.p["blockSizeXY"] = self.param.param["zProject"]["blockSize"]
        self.p["regExp"] =self.param.param["acquisition"]["fileNameRegExp"]
        if 'zBinning' in self.param.param['acquisition']:
            self.p["zBinning"] = int(self.param.param['acquisition']['zBinning'])
        else:
            self.p["zBinning"] = 1
        self.p["zWindow"] = int(self.param.param["zProject"]["zwindows"]/self.p["zBinning"])
        self.p["pixelSizeXY"] = self.param.param["acquisition"]["pixelSizeXY"]
        self.p["pixelSizeZ"] = self.param.param["acquisition"]["pixelSizeZ"]

        if 'parallelizePlanes' in self.param.param['acquisition']:
            self.p["parallelizePlanes"] = self.param.param['acquisition']['parallelizePlanes']
        else:
            self.p["parallelizePlanes"]= 1

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
        self.p["voxel_size_z"] = float(1000*self.p["pixelSizeZ"]*self.p["zBinning"])
        self.p["voxel_size_yx"] = float(1000*self.p["pixelSizeXY"])
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

        if len(properties)>0:
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
            flux=[100*properties[x].max_intensity/threshold for x in selection] # peak intensity over the detection threshold
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

    def segmentSources3D_file(self,fileName2Process):

        p = self.p
        # excludes the reference fiducial and processes files in the same ROI
        roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
        label = str(self.param.decodesFileParts(os.path.basename(fileName2Process))["cycle"])

        # creates Table that will hold results
        outputTable = self.createsOutputTable()

        # - load  and preprocesses 3D fiducial file
        print("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi,label))
        print("$ File:{}".format(os.path.basename(fileName2Process)))
        image3D0 = io.imread(fileName2Process).squeeze()

        # reinterpolates image in z if necessary
        image3D0 = reinterpolateZ(image3D0, range(0,image3D0.shape[0],p["zBinning"]),mode='remove')

        # restricts analysis to a sub volume containing sources
        focalPlaneMatrix, zRange, _= _reinterpolatesFocalPlane(image3D0,blockSizeXY = p["blockSizeXY"], window=p["zWindow"])
        zOffset = zRange[1][0]
        image3D = image3D0[zRange[1],:,:].copy()

        print("$ Focal plane found: {}, zRange = {}, imageSize = {}".format(zRange[0],zRange[1],image3D.shape))

        # preprocesses image by background substraction and level normalization
        image3D = preProcess3DImage(image3D,
                                    p["lower_threshold"],
                                    p["higher_threshold"],
                                    parallelExecution=self.innerParallelLoop)

        # drifts 3D stack in XY
        if self.dictShiftsAvailable and  label != p["referenceBarcode"]:
            # uses existing shift calculated by alignImages
            try:
                shift = self.dictShifts["ROI:" + roi][label]
                print("> Applying existing XY shift...")
            except KeyError:
                shift = None
                raise SystemExit(
                    "# Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                        "ROI:" + self.ROI, label))

        # applies XY shift to 3D stack
        if label != p["referenceBarcode"]:
            # print("$ Applies shift = {:.2f}".format(shift))
            print("$ Applies shift = [{:.2f} ,{:.2f}]".format(shift[0],shift[1]))
            image3D_aligned = appliesXYshift3Dimages(image3D, shift,parallelExecution=self.innerParallelLoop)
        else:
            print("$ Running reference fiducial cycle: no shift applied!")
            shift = np.array([0.,0.])
            image3D_aligned = image3D

        # segments 3D volumes
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
        print("$ Number of sources detected by image segmentation: {}".format(numberSources))

        if numberSources >0:
            print("> Refits spots using gaussian 3D fittings...")

            print(" > Rescales image values after reinterpolation")
            image3D_aligned = exposure.rescale_intensity(image3D_aligned, out_range=(0, 1)) # removes negative backgrounds


            # calls bigfish to get 3D sub-pixel coordinates based on 3D gaussian fitting
            spots_subpixel = fit_subpixel(image3D_aligned,
                                          spots,
                                          voxel_size_z=p["voxel_size_z"],
                                          voxel_size_yx=p["voxel_size_yx"],
                                          psf_z=p["psf_z"],
                                          psf_yx=p["psf_yx"])

            print(" > Updating table and saving results")
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
        self.log1.info("$ Detected {} ROIs".format(self.numberROIs))
        self.log1.info("$ Number of images to be processed: {}".format(len(self.param.fileList2Process)))

        # loads dicShifts with shifts for all ROIs and all labels
        self.dictShifts, self.dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder, self.log1)

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
                # ROI = ROIList[fileNameReference]
                self.fileName2ProcessList = [x for x in self.param.fileList2Process\
                                        if self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI]
                Nfiles2Process=len(self.fileName2ProcessList)
                print("$ Found {} files in ROI [{}]".format(Nfiles2Process, ROI))
                print("$ [roi:cycle] {}".format(" | ".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                                + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"]) for x in self.fileName2ProcessList])))

                if client is None:
                    self.innerParallelLoop = True

                    for fileIndex, fileName2Process in enumerate(self.fileName2ProcessList): #self.param.fileList2Process):
                        print("\n\n>>>Iteration: {}/{}<<<".format(fileIndex,Nfiles2Process))

                        outputTables.append(self.segmentSources3D_file(fileName2Process))
                else:
                    self.innerParallelLoop = False
                    print("> Aligning {} files using {} workers...".format(Nfiles2Process,len(client.scheduler_info()['workers'])))

                    futures = [client.submit(self.segmentSources3D_file,
                                                     x) for x in self.fileName2ProcessList]

                    outputTables = client.gather(futures)
                    print(" > Retrieving {} results from cluster".format(len(outputTables)))

                # Merges Tables for different cycles and appends results Table to that of previous ROI
                outputTableGlobal = vstack([outputTableGlobal]+outputTables)

        print("$ segmentSources3D procesing time: {}".format(datetime.now() - now))

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
        self.log1.addSimpleText("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        self.log1.addSimpleText("$ folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.log1.fileNameMD, "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.label = self.param.param["acquisition"]["label"]
        self.outputFileName = self.dataFolder.outputFiles["segmentedObjects"] + "_3D_" + self.label + ".dat"

        self.log1.addSimpleText("> Processing Folder: {}".format(self.currentFolder))
        self.log1.parallel = self.parallel

        self.segmentSources3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        self.log1.report("$ segmentedObjects run in {} finished".format(self.currentFolder), "info")

        return 0

