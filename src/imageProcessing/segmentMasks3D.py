#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 14:31:47 2021

@author: marcnol
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Purpose: Segments Masks in 3D using AI (stardist)

steps:
    - iterate over ROIs
    - load 3D file for cycle <i>
    - load global alignment for this cycle
    - re-align 3D image using XY alignment
    - segment objects to get labeled masks

    - display results:
    - output results

"""

# =============================================================================
# IMPORTS
# =============================================================================

import glob, os
import numpy as np
from datetime import datetime
from skimage import io

from imageProcessing.imageProcessing import (
    appliesXYshift3Dimages,
    reinterpolateZ,
    plotRawImagesAndLabels,
    _segments3DMasks,
)
from fileProcessing.fileManagement import folders, writeString2File, printDict, getDictionaryValue
from fileProcessing.fileManagement import loadsAlignmentDictionary, retrieveNumberROIsFolder, printLog

from skimage.measure import regionprops

# =============================================================================
# CLASSES
# =============================================================================

class segmentMasks3D:
    def __init__(self, param, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.parallel = parallel

        self.p = dict()

        # parameters from infoList.json
        self.p["referenceBarcode"] = self.param.param["alignImages"]["referenceFiducial"]
        self.p["brightest"] = self.param.param["segmentedObjects"]["brightest"]
        self.p["blockSizeXY"] = self.param.param["zProject"]["blockSize"]
        self.p["regExp"] = self.param.param["acquisition"]["fileNameRegExp"]
        self.p["zBinning"] = getDictionaryValue(self.param.param["acquisition"], "zBinning", default=1)

        self.p["zWindow"] = int(self.param.param["zProject"]["zwindows"] / self.p["zBinning"])
        self.p["pixelSizeXY"] = self.param.param["acquisition"]["pixelSizeXY"]
        self.p["pixelSizeZ"] = self.param.param["acquisition"]["pixelSizeZ"]

        self.p["parallelizePlanes"] = getDictionaryValue(
            self.param.param["acquisition"], "parallelizePlanes", default=1
        )

        # decides what segmentation method to use
        self.p["3Dmethod"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "3Dmethod", default="thresholding"
        )

        # parameters used for 3D segmentation and deblending
        self.p["threshold_over_std"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "3D_threshold_over_std", default=1
        )
        self.p["sigma"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_sigma", default=3)
        boxSize = getDictionaryValue(self.param.param["segmentedObjects"], "3D_boxSize", default=32)
        self.p["boxSize"] = (boxSize, boxSize)
        filter_size = getDictionaryValue(self.param.param["segmentedObjects"], "3D_filter_size", default=3)
        self.p["filter_size"] = (filter_size, filter_size)
        self.p["area_min"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_area_min", default=3)
        self.p["area_max"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_area_max", default=1000)
        self.p["nlevels"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_nlevels", default=64)
        self.p["contrast"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_contrast", default=0.001)

        # parameters for stardist
        self.p["stardist_basename"] = getDictionaryValue(
            self.param.param["segmentedObjects"],
            "stardist_basename",
            default="/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
        ).rstrip('/')

        self.p["stardist_network"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "stardist_network", default="stardist_18032021_single_loci"
        )
        
        self.p["stardist_network3D"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "stardist_network3D", default="stardist_20210625_deconvolved"
        )

        # parameters used for 3D gaussian fitting
        self.p["voxel_size_z"] = float(1000 * self.p["pixelSizeZ"] * self.p["zBinning"])
        self.p["voxel_size_yx"] = float(1000 * self.p["pixelSizeXY"])
        self.p["psf_z"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_psf_z", default=500)
        self.p["psf_yx"] = getDictionaryValue(self.param.param["segmentedObjects"], "3D_psf_yx", default=200)

        # range used for adjusting image levels during pre-precessing
        self.p["lower_threshold"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "3D_lower_threshold", default=0.9
        )
        self.p["higher_threshold"] = getDictionaryValue(
            self.param.param["segmentedObjects"], "3D_higher_threshold", default=0.9999
        )

        # parameters used for plotting 3D image
        # sets the number of planes around the center of the image used to represent localizations in XZ and ZY
        self.p["windowDisplay"] = 10

    def getMaskProperties(self, segmentedImage3D, image3D_aligned, threshold=10, nTolerance=1000):
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

        if len(properties) > 0:
            # selects nTolerance brightest spots and keeps only these for further processing
            peak0 = [x.max_intensity for x in properties]
            peakList = peak0.copy()
            peakList.sort()
            last2keep = np.min([nTolerance, len(peakList)])
            highestPeakValue = peakList[-last2keep]
            selection = list(np.nonzero(peak0 > highestPeakValue)[0])

            # attributes properties using the brightests spots selected
            peak = [properties[x].max_intensity for x in selection]
            centroids = [properties[x].weighted_centroid for x in selection]
            sharpness = [float(properties[x].filled_area / properties[x].bbox_area) for x in selection]
            roundness1 = [properties[x].equivalent_diameter for x in selection]
            roundness2 = [properties[x].extent for x in selection]
            npix = [properties[x].area for x in selection]
            sky = [0.0 for x in selection]
            peak = [properties[x].max_intensity for x in selection]
            flux = [
                100 * properties[x].max_intensity / threshold for x in selection
            ]  # peak intensity over the detection threshold
            mag = [-2.5 * np.log10(x) for x in flux]  # -2.5 log10(flux)

            # converts centroids to spot coordinates for bigfish to run 3D gaussian fits
            z = [x[0] for x in centroids]
            y = [x[1] for x in centroids]
            x = [x[2] for x in centroids]

            spots = np.zeros((len(z), 3))
            spots[:, 0] = z
            spots[:, 1] = y
            spots[:, 2] = x
            spots = spots.astype("int64")

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
            return [], [], [], [], [], [], [], [], []

    def plotsImage3D(self, image3D, masks, normalize=False):
        """
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

        """
        window = self.p["windowDisplay"]

        fig1 = plotRawImagesAndLabels(image3D, masks, normalize=normalize, window=window)

        return fig1

    def _segments3Dvolumes(self, image3D_aligned):
        p = self.p

        binary, segmentedImage3D = _segments3DMasks(
            image3D_aligned,
            axis_norm=(0, 1, 2),
            pmin=1,
            pmax=99.8,
            model_dir=p["stardist_basename"],
            model_name=p["stardist_network3D"],
        )

        return binary, segmentedImage3D

    def segmentMasks3D_file(self, fileName2Process):

        p = self.p
        # excludes the reference fiducial and processes files in the same ROI
        roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
        label = str(self.param.decodesFileParts(os.path.basename(fileName2Process))["cycle"])

        # load  and preprocesses 3D fiducial file
        printLog("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi,label))
        printLog("$ File:{}".format(os.path.basename(fileName2Process)))
        image3D0 = io.imread(fileName2Process).squeeze()

        # reinterpolates image in z if necessary
        image3D = reinterpolateZ(image3D0, range(0, image3D0.shape[0], p["zBinning"]), mode="remove")

        # drifts 3D stack in XY
        if self.dictShiftsAvailable:
            # uses existing shift calculated by alignImages
            try:
                shift = self.dictShifts["ROI:" + roi][label]
                printLog("> Applying existing XY shift...")
            except KeyError:
                shift = None
                raise SystemExit(
                    "# Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                        "ROI:" + self.ROI, label
                    )
                )

        # applies XY shift to 3D stack
        if label != p["referenceBarcode"]:
            printLog("$ Applies shift = [{:.2f} ,{:.2f}]".format(shift[0], shift[1]))
            image3D_aligned = appliesXYshift3Dimages(image3D, shift, parallelExecution=self.innerParallelLoop)
        else:
            printLog("$ Running reference fiducial cycle: no shift applied!")
            shift = np.array([0.0, 0.0])
            image3D_aligned = image3D

        # segments 3D volumes
        binary, segmentedImage3D = self._segments3Dvolumes(image3D_aligned)

        numberMasks = np.max(segmentedImage3D)
        printLog("$ Number of masks detected: {}".format(numberMasks))

        if numberMasks > 0:
            outputExtension = "_3Dmasks"
            NPY_labeled_image_fileName = (
                self.dataFolder.outputFolders["segmentedObjects"] + os.sep + os.path.basename(fileName2Process)
            )
            NPY_labeled_image_fileName = NPY_labeled_image_fileName.split(".")[0] + "." + outputExtension + ".npy"
            printLog(" > Saving output labeled image: {}".format(NPY_labeled_image_fileName))
            np.save(NPY_labeled_image_fileName, segmentedImage3D)

            # represents image in 3D with localizations
            printLog("> plotting outputs...")

            figures = list()
            figures.append([self.plotsImage3D(image3D_aligned, segmentedImage3D,), outputExtension + ".png"])

            # saves figures
            outputFileNames = [
                self.dataFolder.outputFolders["segmentedObjects"] + os.sep + os.path.basename(fileName2Process) + x[1]
                for x in figures
            ]

            for fig, file in zip(figures, outputFileNames):
                fig[0].savefig(file)

        del image3D_aligned, image3D, image3D0

    def segmentMasks3DinFolder(self):
        """
        Segments 3D Masks in all files in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()

        # Reads list of parameters assigned upon Class initialization
        p = self.p
        printDict(p)

        # Finds images to process
        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)
        self.ROIList = retrieveNumberROIsFolder(self.currentFolder, p["regExp"], ext="tif")
        self.numberROIs = len(self.ROIList)
        printLog("$ Detected {} ROIs".format(self.numberROIs))
        printLog("$ Images to be processed: {}".format(self.param.fileList2Process))
        printLog("$ Number of images to be processed: {}".format(len(self.param.fileList2Process)))

        # loads dicShifts with shifts for all ROIs and all labels
        self.dictShifts, self.dictShiftsAvailable = loadsAlignmentDictionary(self.dataFolder)

        if self.numberROIs > 0:

            # loops over ROIs
            for ROI in self.ROIList:

                # loads reference fiducial image for this ROI
                self.fileName2ProcessList = [
                    x
                    for x in self.param.fileList2Process
                    if self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI
                    and (
                        "DAPI" in self.param.decodesFileParts(os.path.basename(x))["cycle"]
                        or "mask" in self.param.decodesFileParts(os.path.basename(x))["cycle"]
                    )
                ]
                Nfiles2Process = len(self.fileName2ProcessList)
                printLog("$ Found {} files in ROI [{}]".format(Nfiles2Process, ROI))
                printLog(
                    "$ [roi:cycle] {}".format(
                        " | ".join(
                            [
                                str(self.param.decodesFileParts(os.path.basename(x))["roi"])
                                + ":"
                                + str(self.param.decodesFileParts(os.path.basename(x))["cycle"])
                                for x in self.fileName2ProcessList
                            ]
                        )
                    )
                )

                self.innerParallelLoop = True
                # processes files in this ROI
                for fileIndex, fileName2Process in enumerate(self.fileName2ProcessList):
                    printLog("\n\n>>>Iteration: {}/{}<<<".format(fileIndex, Nfiles2Process))
                    self.segmentMasks3D_file(fileName2Process)

        printLog("$ segmentMasks3D procesing time: {}".format(datetime.now() - now))

    def segmentMasks3D(self):
        """
        segments 3D masks in rootFolder

        Returns
        -------
        None.

        """
        sessionName = "segmentMasks3D"

        # processes folders and files
        printLog("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        printLog("$ folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.param.param["fileNameMD"], "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.label = self.param.param["acquisition"]["label"]

        printLog("> Processing Folder: {}".format(self.currentFolder))

        self.segmentMasks3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        printLog("$ segmentedObjects run in {} finished".format(self.currentFolder))

        return 0


