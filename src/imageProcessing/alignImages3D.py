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
    makesShiftMatrixHiRes,
    preProcess3DImage,
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

    def createsOutputTable(self):
        return Table(
                names=(
                    "reference file",
                    "aligned file",
                    "blockXY",
                    "ROI #",
                    "label",
                    "block_i",
                    "block_j",
                    "shift_z",
                    "shift_x",
                    "shift_y",
                    "quality_xy",
                    "quality_zy",
                    "quality_zx",
                ),
                dtype=("S2", "S2", "int", "int", "S2", "int", "int","f4", "f4", "f4", "f4", "f4", "f4"),
            )
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

        # loads dicShifts with shifts for all ROIs and all labels
        dictShifts, dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder, self.log1)
        # dictFileName = os.path.splitext(self.dataFolder.outputFiles["dictShifts"])[0] + ".json"

        # # dictFileName = dataFolder.outputFiles["dictShifts"] + ".json"
        # dictShifts = loadJSON(dictFileName)
        # if len(dictShifts) == 0:
        #     self.log1.report("File with dictionary not found!: {}".format(dictFileName))
        #     dictShiftsAvailable = False
        # else:
        #     self.log1.report("Dictionary File loaded: {}".format(dictFileName))
        #     dictShiftsAvailable = True

        # creates Table that will hold results
        alignmentResultsTable=self.createsOutputTable()

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
                    label = os.path.basename(fileName2Process).split("_")[2]  # to FIX

                    if (fileName2Process not in fileNameReference) and roi == ROI:

                        # - load  and preprocesses 3D fiducial file
                        print("\n\nProcessing cycle {}".format(os.path.basename(fileName2Process)))
                        image3D0 = io.imread(fileName2Process).squeeze()
                        image3D = preProcess3DImage(image3D0, self.lower_threshold, self.higher_threshold)

                        # shows original images and background substracted
                        images0 = [imageRef0,image3D0] # list with unprocessed 3D stacks
                        images = [imageRef,image3D] # list with processed 3D stacks
                        allimages = images0 + images
                        fig1 = plots4images(allimages, titles=['reference','cycle <i>','processed reference','processed cycle <i>'])

                        # drifts 3D stack in XY
                        if dictShiftsAvailable:
                            # uses existing shift calculated by alignImages
                            try:
                                shift = dictShifts["ROI:" + roi][label]
                                print("Applying existing XY shift...")
                            except KeyError:
                                shift = None
                                self.log1.report(
                                    "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                                        "ROI:" + ROI, label
                                    ),
                                    "ERROR",
                                )
                        if not dictShiftsAvailable or shift == None:
                            # if dictionary of shift or key for this cycle was not found, then it will recalculate XY shift
                            images_2D = [np.sum(x, axis=0) for x in images]

                            print("Calculating XY shift...")
                            shift, error, diffphase = phase_cross_correlation(images_2D[0], images_2D[1], upsample_factor=self.upsample_factor)

                        # applies XY shift to 3D stack
                        print("shifts XY = {}".format(shift))

                        # images_2D.append(shiftImage(images_2D[1], shift))

                        # reinterpolate second file in XY using dictionnary to get rough alignment
                        images.append(appliesXYshift3Dimages(images[1], shift))

                        # 3D image alignment by block
                        shiftMatrices, block_ref, block_target = imageBlockAlignment3D(images, blockSizeXY=self.blockSizeXY,\
                                                                                       upsample_factor=self.upsample_factor)

                        # [plots shift matrices]
                        fig2 = plots3DshiftMatrices(shiftMatrices, fontsize=8)

                        # combines blocks into a single matrix for display instead of plotting a matrix of subplots each with a block
                        outputs = []
                        for axis in self.axes2Plot:
                            outputs.append(combinesBlocksImageByReprojection(block_ref, block_target, shiftMatrices=shiftMatrices, axis1=axis))

                        SSIM_matrices = [x[1] for x in outputs]
                        MSE_matrices = [x[2] for x in outputs]
                        NRMSE_matrices = [x[3] for x in outputs]

                        fig3 = plt.figure(constrained_layout=False)
                        fig3.set_size_inches((20 * 2, 20))
                        gs = fig3.add_gridspec(2, 2)
                        ax = [fig3.add_subplot(gs[:, 0]), fig3.add_subplot(gs[0, 1]), fig3.add_subplot(gs[1, 1])]

                        titles = ["Z-projection", "X-projection", "Y-projection"]

                        for axis, output, i in zip(ax, outputs, self.axes2Plot):
                            axis.imshow(output[0])
                            axis.set_title(titles[i])

                        fig3.tight_layout()

                        fig4 = plots3DshiftMatrices(SSIM_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        fig4.suptitle("SSIM block matrices")

                        fig5 = plots3DshiftMatrices(MSE_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        fig5.suptitle("mean square root block matrices")

                        fig6 = plots3DshiftMatrices(NRMSE_matrices, fontsize=6, log=False,valfmt="{x:.2f}")
                        fig6.suptitle("normalized root mean square block matrices")

                        # saves figures
                        figTitles = ['_bkgSubstracted.png','_shiftMatrices.png',
                                     '_3Dalignments.png','_SSIMblocks.png',
                                     '_MSEblocks.png','_NRMSEblocks.png']
                        outputFileNames = ['/home/marcnol/Documents'+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]
                        outputFileNames = [self.dataFolder.outputFolders["alignImages"]+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]

                        figs=[fig1,fig2,fig3,fig4,fig5,fig6]
                        for fig, file in zip(figs,outputFileNames):
                            fig.savefig(file)

                        # dict with shiftMatrix and NRMSEmatrix: https://en.wikipedia.org/wiki/Root-mean-square_deviation
                        # These matrices can be used to apply and assess zxy corrections for any pixel in the 3D image
                        # reference file,aligned file,ROI,label,block_i,block_j,shift_z,shift_x,shift_y,quality_xy,quality_zy,quality_zx
                        numBlocks,blockXY = block_ref.shape[0], block_ref.shape[-1]
                        for i in range(numBlocks):
                            for j in range(numBlocks):
                                Table_entry = [os.path.basename(fileNameReference),
                                               os.path.basename(fileName2Process),
                                               int(blockXY),
                                               int(roi),
                                               label,
                                               i,
                                               j,
                                               shiftMatrices[0][i,j],
                                               shiftMatrices[1][i,j],
                                               shiftMatrices[2][i,j],
                                               NRMSE_matrices[0][i,j],
                                               NRMSE_matrices[1][i,j],
                                               NRMSE_matrices[2][i,j],
                                               ]
                                alignmentResultsTable.add_row(Table_entry)

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


# =============================================================================
#   FUNCTIONS
# =============================================================================
