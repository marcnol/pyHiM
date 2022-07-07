# -*- coding: utf-8 -*-
"""
Spyder Editor


Purpose: Corrects drift in 3D

The drift correction routines in alignImages.py take care of the corrections in XY but not in Z.

Drift in the z-position from cycle to cycle is small, typically 100-200 nm, but needs to be checked
and corrected. In addition, sample deformation can lead to inhomogeneous drift that cannot be corrected by a rigid body translation.

This routine solves both issues by correcting drift in 3D by block decomposition.

steps:
    - iterate over ROIs
    - load 3D fiducial file for reference fiducial
    - iterate over cycles <i>
    - load 3D fiducial file for fiducial barcode <i>
    - re-align 3D fiducial image using XY alignment
    - perform block alignment in 3D by cross=correlating blocks in 3D.

    - store in database.
    - display results:
        - drift correction maps in X-Y-Z
        - corrected blocks in XY, ZX, ZY

During buildMatrix this database is loaded, if available.
    - check database exist and load it
    - correct z-coordinate of the barcode provided the correction given in the dict

"""
# =============================================================================
# IMPORTS
# =============================================================================

import glob, os, time
import matplotlib.pylab as plt
import numpy as np
from datetime import datetime
from skimage import io
import resource
# import pympler
from pympler import tracker

from astropy.table import Table, vstack

from imageProcessing.imageProcessing import (
    appliesXYshift3Dimages,
    imageBlockAlignment3D,
    plots3DshiftMatrices,
    combinesBlocksImageByReprojection,
    plots4images,
    preProcess3DImage,
    reinterpolateZ,
)
from fileProcessing.fileManagement import folders, writeString2File, getDictionaryValue
from fileProcessing.fileManagement import RT2fileName, loadsAlignmentDictionary
from fileProcessing.fileManagement import try_get_client, printDict, printLog

from skimage.registration import phase_cross_correlation


# =============================================================================
# CLASSES
# =============================================================================

class drift3D:
    def __init__(self, param, session1, parallel=False):
        self.param = param
        self.session1 = session1
        self.window = 3
        self.parallel = parallel
        self.p = dict()

        self.p["blockSizeXY"] = 128
        self.p["upsample_factor"]=100
        self.p["lower_threshold"] = 0.9
        self.p["higher_threshold"]=0.9999

        self.p["blockSizeXY"]=getDictionaryValue(self.param.param["alignImages"], "blockSizeXY", default=128)
        self.p["lower_threshold"]=getDictionaryValue(self.param.param["alignImages"], "3D_lower_threshold", default=0.9)
        self.p["higher_threshold"]=getDictionaryValue(self.param.param["alignImages"], "3D_higher_threshold", default=0.9999)

        self.p["axes2Plot"] = range(3)
        self.p["referenceBarcode"] = self.param.param["alignImages"]["referenceFiducial"]

        if 'zBinning' in self.param.param['acquisition']:
            self.p["zBinning"] = self.param.param['acquisition']['zBinning']
        else:
            self.p["zBinning"]= 1

        if 'parallelizePlanes' in self.param.param['acquisition']:
            self.p["parallelizePlanes"] = self.param.param['acquisition']['parallelizePlanes']
        else:
            self.p["parallelizePlanes"]= 1

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

    def alignFiducials3Dfile(self,fileName2Process):
        """
        Aligns <fileName2Process> fiducial against reference

        Returns
        -------
        None.

        """

        p=self.p
        alignmentResultsTable=self.createsOutputTable()

        # excludes the reference fiducial and processes files in the same ROI
        roi = self.param.decodesFileParts(os.path.basename(fileName2Process))["roi"]
        label = str(self.param.decodesFileParts(os.path.basename(fileName2Process))["cycle"])
        innerParallelLoop = self.innerParallelLoop
        imageRef = self.imageRef
        imageRef0 = self.imageRef0
        dictShiftsAvailable = self.dictShiftsAvailable
        dictShifts = self.dictShifts
        outputFolder = self.dataFolder.outputFolders["alignImages"]

        return _alignFiducials3Dfile(fileName2Process,
                                     alignmentResultsTable,
                                     p,
                                     roi,
                                     label,
                                     innerParallelLoop,
                                     imageRef,
                                     imageRef0,
                                     dictShifts,
                                     dictShiftsAvailable,
                                     outputFolder)


    def loadsReferenceFiducial(self,fileNameReference):
        """
        Loads Reference fiducial image and reports on the number of cycles to process for this reference

        Returns
        -------
        None.

        """

        self.p["fileNameReference"] = fileNameReference
        self.p["ROI"]= self.ROIList[fileNameReference]
        printLog("Loading reference 3D image: {}".format(fileNameReference))

        self.imageRef0, self.imageRef= loadNpreprocessImage(fileNameReference,
                                                            self.p["zBinning"],
                                                            self.p["lower_threshold"],
                                                            self.p["higher_threshold"],
                                                            parallelExecution=False)

        self.imageRef0 = np.sum(self.imageRef0, axis=0) # replaces 3D with a 2D projection

        self.fileName2ProcessList = [x for x in self.param.fileList2Process\
                                if (x not in fileNameReference) and \
                                    self.param.decodesFileParts(os.path.basename(x))["roi"] == self.p["ROI"]]

        printLog("$ Found {} files in ROI: {}".format(len(self.fileName2ProcessList), self.p["ROI"]))
        printLog("$ [roi:cycle] {}".format("|".join([str(self.param.decodesFileParts(os.path.basename(x))["roi"])\
                        + ":" + str(self.param.decodesFileParts(os.path.basename(x))["cycle"])\
                            for x in self.fileName2ProcessList])))


    def loadsDictShifts(self):
        """
        Lods dictionary of XY shifts

        Returns
        -------
        None.

        """
        printLog("\nReference Barcode: {}".format(self.p["referenceBarcode"]))
        self.fileNameReferenceList, self.ROIList = RT2fileName(self.param, self.p["referenceBarcode"])

        self.numberROIs = len(self.ROIList)
        printLog("\nDetected {} ROIs".format(self.numberROIs))

        # loads dicShifts with shifts for all ROIs and all labels
        self.dictShifts, self.dictShiftsAvailable  = loadsAlignmentDictionary(self.dataFolder)


    def alignFiducials3DinFolder(self):
        """
        Refits all the barcode files found in rootFolder

        Returns
        -------
        None.

        """
        now = datetime.now()
        printDict(self.p)
        tr = tracker.SummaryTracker()

        # gets files to process
        filesFolder = glob.glob(self.currentFolder + os.sep + "*.tif")
        self.param.files2Process(filesFolder)

        # loads dictinary of shifts
        self.loadsDictShifts()

        # creates Table that will hold results
        alignmentResultsTableGlobal=self.createsOutputTable()
        alignmentResultsTables = list()

        if self.p["parallelizePlanes"]:
            client = None
        else:
            client = try_get_client()

        if self.numberROIs > 0:

            # loops over ROIs
            for fileNameReference in self.fileNameReferenceList:

                # loads reference fiducial image for this ROI
                self.loadsReferenceFiducial(fileNameReference)
                numberFiles = len(self.fileName2ProcessList)
                tr.print_diff()

                # for fileIndex, fileName2Process in enumerate(self.param.fileList2Process):
                if client is None:
                    self.innerParallelLoop = True
                    for fileIndex, fileName2Process in enumerate(self.fileName2ProcessList):

                        printLog("\n\n>>>Iteration: {}/{}<<<".format(fileIndex,numberFiles))

                        alignmentResultsTables.append(self.alignFiducials3Dfile(fileName2Process))

                        tr.print_diff()

                else:
                    self.innerParallelLoop = False
                    printLog("> Aligning {} files using {} workers...".format(numberFiles,len(client.scheduler_info()['workers'])))

                    futures = [client.submit(self.alignFiducials3Dfile,
                                                     x) for x in self.fileName2ProcessList]

                    alignmentResultsTables = client.gather(futures)
                    printLog(" > Retrieving {} results from cluster".format(len(alignmentResultsTables)))

                    # del futures

                # Merges Tables for different cycles and appends results Table to that of previous ROI
                alignmentResultsTableGlobal = vstack([alignmentResultsTableGlobal]+alignmentResultsTables)

        # saves Table with all shifts
        output_filename = self.dataFolder.outputFiles["alignImages"].split(".")[0] + "_block3D.dat"
        alignmentResultsTableGlobal.write(
            output_filename,
            format="ascii.ecsv",
            overwrite=True,
        )

        printLog("$ alignImages3D procesing time: {}".format(datetime.now() - now))
        printLog(f"$ alignImages3D output Table saved in: {output_filename}")
        
    def alignFiducials3D(self):
        """
        runs refitting routine in rootFolder

        Returns
        -------
        None.

        """
        sessionName = "alignImages3D"

        # processes folders and files
        printLog("\n===================={}====================\n".format(sessionName))
        self.dataFolder = folders(self.param.param["rootFolder"])
        printLog("folders read: {}".format(len(self.dataFolder.listFolders)))
        writeString2File(self.param.param["fileNameMD"], "## {}\n".format(sessionName), "a")

        # creates output folders and filenames
        self.currentFolder = self.dataFolder.listFolders[0]

        self.dataFolder.createsFolders(self.currentFolder, self.param)
        self.outputFileName = self.dataFolder.outputFiles["alignImages"]

        printLog("-------> Processing Folder: {}".format(self.currentFolder))
        # self.log1.parallel = self.parallel

        self.alignFiducials3DinFolder()

        self.session1.add(self.currentFolder, sessionName)

        printLog("HiM matrix in {} processed".format(self.currentFolder))

        return 0


# =============================================================================
#   FUNCTIONS
# =============================================================================

def loadNpreprocessImage(fileName2Process, zBinning, lower_threshold, higher_threshold,parallelExecution=True):

    printLog("$ File:{}".format(os.path.basename(fileName2Process)))

    image3D0 = io.imread(fileName2Process).squeeze()

    # reinterpolates image in z if necessary
    image3D0 = reinterpolateZ(image3D0, range(0,image3D0.shape[0],zBinning),mode='remove')

    image3D = preProcess3DImage(image3D0, lower_threshold, higher_threshold,parallelExecution=parallelExecution)

    return image3D0, image3D

def _alignFiducials3Dfile(fileName2Process,
                                 alignmentResultsTable,
                                 p,
                                 roi,
                                 label,
                                 innerParallelLoop,
                                 imageRef,
                                 imageRef0,
                                 dictShifts,
                                 dictShiftsAvailable,
                                 outputFolder):

    baseMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000

    # - load  and preprocesses 3D fiducial file
    printLog("\n\n>>>Processing roi:[{}] cycle:[{}]<<<".format(roi,label))
    image3D0, image3D = loadNpreprocessImage(fileName2Process,
                                             p["zBinning"],
                                             p["lower_threshold"],
                                             p["higher_threshold"],
                                             parallelExecution=innerParallelLoop)

    # shows original images and background substracted
    image3D0 = np.sum(image3D0,axis=0) # replaces by a 2D projection
    images = [imageRef,image3D]
    images2D=[np.sum(x,axis=0) for x in images]
    fig1 = plots4images([imageRef0,image3D0]+images2D, titles=['reference','cycle <i>','processed reference','processed cycle <i>'])

    del image3D0

    # drifts 3D stack in XY
    # ---------------------
    if dictShiftsAvailable:
        # uses existing shift calculated by alignImages
        try:
            shift = dictShifts["ROI:" + roi][label]
            printLog("> Applying existing XY shift...")
        except KeyError:
            shift = None
            printLog(
                "Could not find dictionary with alignment parameters for this ROI: {}, label: {}".format(
                    "ROI:" + roi, label
                ),
                status = "WARN",
            )
    if not dictShiftsAvailable or shift == None:
        # if dictionary of shift or key for this cycle was not found, then it will recalculate XY shift
        images_2D = [np.sum(x, axis=0) for x in images]

        printLog("> Calculating XY shift...")
        shift, error, diffphase = phase_cross_correlation(images_2D[0], images_2D[1], upsample_factor=p["upsample_factor"])

    # applies XY shift to 3D stack
    # ----------------------------
    printLog("$ shifts XY = {}".format(shift))

    # reinterpolate second file in XY using dictionnary to get rough alignment
    images.append(appliesXYshift3Dimages(image3D, shift,parallelExecution=innerParallelLoop))

    del images[1], image3D # removes unshifted image to save memory

    # 3D image alignment by block
    # ---------------------------
    printLog("> Block-aligning images in 3D...")
    shiftMatrices, block_ref, block_target = imageBlockAlignment3D(images, blockSizeXY=p["blockSizeXY"],\
                                                                   upsample_factor=p["upsample_factor"])
    del images # deletes image list to save memory

    # [plots shift matrices]
    fig2 = plots3DshiftMatrices(shiftMatrices, fontsize=8)

    # combines blocks into a single matrix for display instead of plotting a matrix of subplots each with a block
    outputs = []
    for axis in p["axes2Plot"]:
        outputs.append(combinesBlocksImageByReprojection(block_ref, block_target, shiftMatrices=shiftMatrices, axis1=axis))

    SSIM_matrices = [x[1] for x in outputs]
    MSE_matrices = [x[2] for x in outputs]
    NRMSE_matrices = [x[3] for x in outputs]

    fig3 = plt.figure(constrained_layout=False)
    fig3.set_size_inches((20 * 2, 20))
    gs = fig3.add_gridspec(2, 2)
    ax = [fig3.add_subplot(gs[:, 0]), fig3.add_subplot(gs[0, 1]), fig3.add_subplot(gs[1, 1])]

    titles = ["Z-projection", "X-projection", "Y-projection"]

    for axis, output, i in zip(ax, outputs, p["axes2Plot"]):
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
    # -------------
    figTitles = ['_bkgSubstracted.png','_shiftMatrices.png',
                 '_3Dalignments.png','_SSIMblocks.png',
                 '_MSEblocks.png','_NRMSEblocks.png']
    outputFileNames = [outputFolder+os.sep+os.path.basename(fileName2Process)+x for x in figTitles]

    figs=[fig1,fig2,fig3,fig4,fig5,fig6]
    for fig, file in zip(figs,outputFileNames):
        fig.savefig(file)
        plt.close(fig)

    # Saves results
    # -------------
    # dict with shiftMatrix and NRMSEmatrix: https://en.wikipedia.org/wiki/Root-mean-square_deviation
    # These matrices can be used to apply and assess zxy corrections for any pixel in the 3D image
    # reference file,aligned file,ROI,label,block_i,block_j,shift_z,shift_x,shift_y,quality_xy,quality_zy,quality_zx
    numBlocks,blockXY = block_ref.shape[0], block_ref.shape[-1]
    for i in range(numBlocks):
        for j in range(numBlocks):
            Table_entry = [os.path.basename(p["fileNameReference"]),
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

    printLog("$ delta memory used: {} Mb".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000-baseMemory))
    printLog("Erasing {} variables\n".format(len(dir())-1))
    for var in dir():
        if var != "alignmentResultsTable":
            del var

    printLog("$ Memory used after variables deleted: {} Mb".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000-baseMemory))

    printLog("Variables still alive: {}".format(dir()))

    return alignmentResultsTable