# -*- coding: utf-8 -*-
"""
Spyder Editor


Purpose: Corrects drift in 3D

Often there is drift in the z-position from cycle to cycle.

The drift correction routines take care of the corrections in XY but not in Z.

steps:
    - iterate over ROIs
    - load 3D fiducial file for reference fiducial
    - iterate over cycles <i>
    - load 3D fiducial file for fiducial barcode <i>
    - re-align 3D fiducial image using XY alignment
    - perform block alignment in 3D by cross=correlating blocks in 3D.

    ** determine if we keep or not based on conditions to be tested **

    - store in database.
    - display results:
        - drift correction maps in X-Y-Z
        - corrected blocks in XY, ZX, ZY

During buildMatrix, if available, this database is loaded.
    - check database exist and load it
    - correct z-coordinate of the barcode provided the correction given in the dict




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
)
from fileProcessing.fileManagement import folders, writeString2File, loadJSON
from fileProcessing.fileManagement import daskCluster, RT2fileName

from photutils import Background2D, MedianBackground
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats

from astropy.stats import SigmaClip
from imageProcessing.segmentMasks import _showsImageSources
from tqdm import trange, tqdm
from skimage import exposure
from skimage.registration import phase_cross_correlation


# =============================================================================
# CLASSES
# =============================================================================


class drift3D:
    def __init__(self, param, log1, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.log1 = log1
        self.window = 3
        self.parallel = parallel

    def findsFile2Process(self, nBarcode, nROI):
        Barcode = "RT" + str(nBarcode)
        ROI = str(nROI) + "_ROI"
        channelbarcode = self.param.setsChannel("barcode_channel", "ch01")

        filesFolder = glob.glob(self.dataFolder.masterFolder + os.sep + "*.tif")
        imageFile = [x for x in filesFolder if ROI in x and Barcode in x and channelbarcode in x]

        return imageFile



    def alignFiducials3DinFolder(self):
        """
        Refits all the barcode files found in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()

        referenceBarcode = self.param.param["alignImages"]["referenceFiducial"]
        self.log1.info("\nReference Barcode: {}".format(referenceBarcode))
        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)

        fileNameReferenceList, ROIList = RT2fileName(self.param, referenceBarcode)

        numberROIs = len(ROIList)
        self.log1.info("\nDetected {} ROIs".format(numberROIs))

        if numberROIs > 0:

            # loops over ROIs
            for fileNameReference in fileNameReferenceList:

                # loads reference fiducial image for this ROI
                ROI = ROIList[fileNameReference]
                self.log1.report("Loading reference 3D image: {}".format(fileNameReference))
                imageRef0 = io.imread(fileNameReference).squeeze()
                imageRef = preProcess3DImage(imageRef0, self.lower_threshold, self.higher_threshold)

                fileName2ProcessList = [x for x in self.param.fileList2Process\
                                        if (x not in fileNameReference) and self.param.decodesFileParts(os.path.basename(x))["roi"] == ROI]

                print("Found {} files in ROI: {}".format(len(fileName2ProcessList), ROI))
                print("[roi:cycle] {}".format("|".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                                + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"]) for x in fileName2ProcessList])))

                for fileName2Process in self.param.fileList2Process:
                    # excludes the reference fiducial and processes files in the same ROI
                    roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]

                    if (fileName2Process not in fileNameReference) and roi == ROI:


                        # - load  and preprocesses 3D fiducial file
                        print("Processing cycle {}".format(os.path.basename(fileName2Process)))
                        image3D0 = io.imread(fileName2Process).squeeze()
                        image3D = preProcess3DImage(image3D0, self.lower_threshold, self.higher_threshold)

                        # shows original images and background substracted
                        images0 = [imageRef0,image3D0] # list with unprocessed 3D stacks
                        images = [imageRef,image3D] # list with processed 3D stacks
                        allimages = images0 + images
                        fig1 = plots4images(allimages, titles=['reference','cycle <i>','processed reference','processed cycle <i>'])

                        # calculates XY shift  and applies it
                        images_2D = [np.sum(x, axis=0) for x in images]

                        print("Calculating shifts...")
                        shift, error, diffphase = phase_cross_correlation(images_2D[0], images_2D[1], upsample_factor=self.upsample_factor)
                        print("shifts XY = {}".format(shift))

                        images_2D.append(shiftImage(images_2D[1], shift))

                        # reinterpolate second file in XY using dictionnary to get rough alignment
                        images.append(appliesXYshift3Dimages(images[1], shift))

                        # 3D image alignment by block

                        shiftMatrices, block_ref, block_target = imageBlockAlignment3D(images, blockSizeXY=self.blockSizeXY, upsample_factor=self.upsample_factor)

                        # [plots shift matrices]
                        fig2 = plots3DshiftMatrices(shiftMatrices, fontsize=8)

                        # combines blocks into a single matrix for display instead of plotting a matrix of subplots each with a block
                        outputs = []
                        for axis in self.axes2Plot:
                            outputs.append(combinesBlocksImageByReprojection(block_ref, block_target, shiftMatrices, axis1=axis))

                        fig3 = plt.figure(constrained_layout=False)
                        fig3.set_size_inches((20 * 2, 20))
                        gs = fig3.add_gridspec(2, 2)
                        ax = [fig3.add_subplot(gs[:, 0]), fig3.add_subplot(gs[0, 1]), fig3.add_subplot(gs[1, 1])]

                        titles = ["Z-projection", "X-projection", "Y-projection"]

                        for axis, output, i in zip(ax, outputs, self.axes2Plot):
                            axis.imshow(output)
                            axis.set_title(titles[i])

                        fig3.tight_layout()

                        # saves figures
                        figTitles = ['_bkgSubstracted.png','_shiftMatrices.png','_3Dalignments.png']
                        outputFileNames = ['/home/marcnol/Documents'+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]

                        figs=[fig1,fig2,fig3]
                        for fig, file in zip(figs,outputFileNames):
                            fig.savefig(file)

                        # saves results to database in dict format and to a Table

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


# =============================================================================
#   FUNCTIONS
# =============================================================================
def preProcess3DImage(x,lower_threshold, higher_threshold):

    # images0= [x/x.max() for x in images0]
    image = exposure.rescale_intensity(x, out_range=(0, 1))

    print("Removing inhomogeneous background...")
    image = _removesInhomogeneousBackground(image)

    print("Rescaling grey levels...")
    image = imageAdjust(image, lower_threshold=lower_threshold, higher_threshold=higher_threshold)[0]

    return image
