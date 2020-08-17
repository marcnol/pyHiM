#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 14:59:43 2020

@author: marcnol

File containing all functions responsible for segmentation of masks for Hi-M,
including DNA masks, barcodes, and fiducials

At the moment, fittings of the 2D positions of barcodes is also performed just 
after image segmentation.

"""


# =============================================================================
# IMPORTS
# =============================================================================

# ---- stardist
from __future__ import print_function, unicode_literals, absolute_import, division

import glob, os
import matplotlib.pylab as plt
import numpy as np
import uuid

from astropy.stats import sigma_clipped_stats, SigmaClip, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, vstack, Column
from photutils import DAOStarFinder, CircularAperture, detect_sources
from photutils import detect_threshold, deblend_sources
from photutils import Background2D, MedianBackground
from photutils.segmentation.core import SegmentationImage

from imageProcessing.imageProcessing import Image, saveImage2Dcmd
from fileProcessing.fileManagement import (
    folders, writeString2File)

# ---- stardist
import matplotlib
matplotlib.rcParams["image.interpolation"] = None
from csbdeep.utils import normalize

from stardist import random_label_cmap
from stardist.models import StarDist2D

np.random.seed(6)
lbl_cmap = random_label_cmap()

# to remove in a future version
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================


def showsImageSources(im, im1_bkg_substracted, log1, sources, outputFileName):
    # show results
    fig = plt.figure()
    fig.set_size_inches((30, 30))

    positions = np.transpose(
        (sources["xcentroid"] + 0.5, sources["ycentroid"] + 0.5)
    )  # for some reason sources are always displays 1/2 px from center of spot

    apertures = CircularAperture(positions, r=4.0)
    norm = simple_norm(im, "sqrt", percent=99.99)
    # norm = ImageNormalize(stretch=SqrtStretch())
    # plt.imshow(im1_bkg_substracted, clim=(0, 1), cmap="Greys", origin="lower", norm=norm)
    plt.imshow(im1_bkg_substracted, cmap="Greys", origin="lower", norm=norm)
    apertures.plot(color="blue", lw=1.5, alpha=0.35)
    plt.xlim(0, im.shape[1] - 1)
    plt.ylim(0, im.shape[0] - 1)
    plt.savefig(outputFileName + "_segmentedSources.png")
    plt.axis("off")
    plt.close()
    writeString2File(
        log1.fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_segmentedSources.png"),
        "a",
    )
    
    # np.save('/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug/segmentedObjects/im1_bkg_substracted.npy',im1_bkg_substracted)


def showsImageMasks(im, log1, segm_deblend, outputFileName):

    norm = ImageNormalize(stretch=SqrtStretch())
    # cmap = segm_deblend.make_cmap(random_state=12345)
    cmap = lbl_cmap

    fig = plt.figure()
    fig.set_size_inches((30, 30))
    plt.imshow(im, cmap="Greys_r", origin="lower", norm=norm)
    plt.imshow(segm_deblend, origin="lower", cmap=cmap, alpha=0.5)
    plt.savefig(outputFileName + "_segmentedMasks.png")
    plt.close()
    writeString2File(
        log1.fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_segmentedMasks.png"),
        "a",
    )


def segmentSourceInhomogBackground(im, param):
    """
    Segments barcodes by estimating inhomogeneous background
    Parameters
    ----------
    im : NPY 2D
        image to be segmented
    param : Parameters
        parameters object.

    Returns
    -------
    table : `astropy.table.Table` or `None`
    A table of found stars with the following parameters:
    
    * ``id``: unique object identification number.
    * ``xcentroid, ycentroid``: object centroid.
    * ``sharpness``: object sharpness.
    * ``roundness1``: object roundness based on symmetry.
    * ``roundness2``: object roundness based on marginal Gaussian
      fits.
    * ``npix``: the total number of pixels in the Gaussian kernel
      array.
    * ``sky``: the input ``sky`` parameter.
    * ``peak``: the peak, sky-subtracted, pixel value of the object.
    * ``flux``: the object flux calculated as the peak density in
      the convolved image divided by the detection threshold.  This
      derivation matches that of `DAOFIND`_ if ``sky`` is 0.0.
    * ``mag``: the object instrumental magnitude calculated as
      ``-2.5 * log10(flux)``.  The derivation matches that of
      `DAOFIND`_ if ``sky`` is 0.0.
    
    `None` is returned if no stars are found.
    
    img_bkc_substracted: 2D NPY array with background substracted image
    """

    threshold_over_std = param.param["segmentedObjects"]["threshold_over_std"]
    fwhm = param.param["segmentedObjects"]["fwhm"]
    brightest = param.param["segmentedObjects"]["brightest"]  # keeps brightest sources

    # estimates inhomogeneous background
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)

    im1_bkg_substracted = im - bkg.background
    mean, median, std = sigma_clipped_stats(im1_bkg_substracted, sigma=3.0)

    # estimates sources
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_over_std * std, brightest=brightest, exclude_border=True,)
    sources = daofind(im1_bkg_substracted)

    return sources, im1_bkg_substracted


def segmentSourceFlatBackground(im, param):
    """
    Segments barcodes using flat background
    Parameters
    ----------
    im : NPY 2D
        image to be segmented
    param : Parameters
        parameters object.

    Returns
    -------
    table : `~astropy.table.Table` or `None`
    A table of found stars with the following parameters:
    
    * ``id``: unique object identification number.
    * ``xcentroid, ycentroid``: object centroid.
    * ``sharpness``: object sharpness.
    * ``roundness1``: object roundness based on symmetry.
    * ``roundness2``: object roundness based on marginal Gaussian
      fits.
    * ``npix``: the total number of pixels in the Gaussian kernel
      array.
    * ``sky``: the input ``sky`` parameter.
    * ``peak``: the peak, sky-subtracted, pixel value of the object.
    * ``flux``: the object flux calculated as the peak density in
      the convolved image divided by the detection threshold.  This
      derivation matches that of `DAOFIND`_ if ``sky`` is 0.0.
    * ``mag``: the object instrumental magnitude calculated as
      ``-2.5 * log10(flux)``.  The derivation matches that of
      `DAOFIND`_ if ``sky`` is 0.0.
    
    `None` is returned if no stars are found.

    img_bkc_substracted: 2D NPY array with background substracted image
    """

    threshold_over_std = param.param["segmentedObjects"]["threshold_over_std"]
    fwhm = param.param["segmentedObjects"]["fwhm"]

    # removes background
    mean, median, std = sigma_clipped_stats(im, param.param["segmentedObjects"]["background_sigma"])
    im1_bkg_substracted = im - median

    # estimates sources
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_over_std * std, exclude_border=True)
    sources = daofind(im - median)

    return sources, im1_bkg_substracted


def segmentMaskInhomogBackground(im, param):
    """
    Function used for segmenting DAPI masks with the ASTROPY library that uses image processing    

    Parameters
    ----------
    im : 2D np array
        image to be segmented.
    param : Parameters class
        parameters.

    Returns
    -------
    segm_deblend: 2D np array where each pixel contains the label of the mask segmented. Background: 0

    """
    # removes background
    threshold = detect_threshold(im, nsigma=2.0)
    sigma_clip = SigmaClip(sigma=param.param["segmentedObjects"]["background_sigma"])

    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
    threshold = bkg.background + (
        param.param["segmentedObjects"]["threshold_over_std"] * bkg.background_rms
    )  # background-only error image, typically 1.0

    sigma = param.param["segmentedObjects"]["fwhm"] * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    # estimates masks and deblends
    segm = detect_sources(im, threshold, npixels=param.param["segmentedObjects"]["area_min"], filter_kernel=kernel,)

    # removes masks too close to border
    segm.remove_border_labels(border_width=10)  # parameter to add to infoList

    segm_deblend = deblend_sources(
        im,
        segm,
        npixels=param.param["segmentedObjects"]["area_min"],  # typically 50 for DAPI
        filter_kernel=kernel,
        nlevels=32,
        contrast=0.001,  # try 0.2 or 0.3
        relabel=True,
    )

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        # print('label {}, with area {}'.format(label,area))
        if area < param.param["segmentedObjects"]["area_min"] or area > param.param["segmentedObjects"]["area_max"]:
            segm_deblend.remove_label(label=label)
            # print('label {} removed'.format(label))

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend


def segmentMaskStardist(im, param):
    """
    Function used for segmenting DAPI masks with the STARDIST package that uses Deep Convolutional Networks

    Parameters
    ----------
    im : 2D np array
        image to be segmented.
    param : Parameters class
        parameters.

    Returns
    -------
    segm_deblend: 2D np array where each pixel contains the label of the mask segmented. Background: 0

    """
    # removes background
    threshold = detect_threshold(im, nsigma=2.0)
    sigma_clip = SigmaClip(sigma=param.param["segmentedObjects"]["background_sigma"])

    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
    threshold = bkg.background + (
        param.param["segmentedObjects"]["threshold_over_std"] * bkg.background_rms
    )  # background-only error image, typically 1.0

    sigma = param.param["segmentedObjects"]["fwhm"] * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    # outputFileName = dataFolder.outputFolders["segmentedObjects"] + os.sep + rootFileName

    n_channel = 1 if im.ndim == 2 else im.shape[-1]
    axis_norm = (0, 1)  # normalize channels independently

    if n_channel > 1:
        print(
            "Normalizing image channels %s." % ("jointly" if axis_norm is None or 2 in axis_norm else "independently")
        )

    model = StarDist2D(
        None,
        name=param.param["segmentedObjects"]["stardist_network"],
        basedir=param.param["segmentedObjects"]["stardist_basename"],
    )

    img = normalize(im, 1, 99.8, axis=axis_norm)
    labeled, details = model.predict_instances(img)

    if True:
        plt.figure(figsize=(8, 8))
        plt.imshow(img, clim=(0, 1), cmap="gray")
        plt.imshow(labeled, cmap=lbl_cmap, alpha=0.5)
        plt.axis("off")

    # estimates masks and deblends
    segm = SegmentationImage(labeled)
    # threshold = 0.5
    # segm = detect_sources(labeled, threshold, npixels=param.param["segmentedObjects"]["area_min"], filter_kernel=kernel,)

    # removes masks too close to border
    segm.remove_border_labels(border_width=10)  # parameter to add to infoList
    segm_deblend = segm

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        # print('label {}, with area {}'.format(label,area))
        if area < param.param["segmentedObjects"]["area_min"] or area > param.param["segmentedObjects"]["area_max"]:
            segm_deblend.remove_label(label=label)
            # print('label {} removed'.format(label))

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend, labeled


def makesSegmentations(fileName, param, log1, session1, dataFolder):

    rootFileName = os.path.basename(fileName).split(".")[0]
    outputFileName = dataFolder.outputFolders["segmentedObjects"] + os.sep + rootFileName
    fileName_2d_aligned = dataFolder.outputFolders["alignImages"] + os.sep + rootFileName + "_2d_registered.npy"

    print("searching for {}".format(fileName_2d_aligned))
    if (
        # (fileName in session1.data)
        param.param["segmentedObjects"]["operation"] == "overwrite"
        and os.path.exists(fileName_2d_aligned)
    ):  # file exists

        ROI = os.path.basename(fileName).split("_")[param.param["acquisition"]["positionROIinformation"]]
        label = param.param["acquisition"]["label"]

        # loading registered 2D projection
        Im = Image()
        Im.loadImage2D(
            fileName, log1, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
        )
        im = Im.data_2D
        log1.report("[{}] Loaded 2D registered file: {}".format(label, os.path.basename(fileName)))

        if label == "barcode" and len([i for i in rootFileName.split("_") if "RT" in i]) > 0:
            if param.param["segmentedObjects"]["background_method"] == "flat":
                output = segmentSourceFlatBackground(im, param)
            elif param.param["segmentedObjects"]["background_method"] == "inhomogeneous":
                output, im1_bkg_substracted = segmentSourceInhomogBackground(im, param)
            else:
                log1.report(
                    "segmentedObjects/background_method not specified in json file", "ERROR",
                )
                return Table()

            # show results
            showsImageSources(im, im1_bkg_substracted, log1, output, outputFileName)

            # [ formats results Table for output by adding buid, barcodeID, CellID and ROI]

            # buid
            buid = []
            for i in range(len(output)):
                buid.append(str(uuid.uuid4()))
            colBuid = Column(buid, name="Buid", dtype=str)

            # barcodeID, cellID and ROI
            barcodeID = os.path.basename(fileName).split("_")[2].split("RT")[1]
            colROI = Column(int(ROI) * np.ones(len(output)), name="ROI #", dtype=int)
            colBarcode = Column(int(barcodeID) * np.ones(len(output)), name="Barcode #", dtype=int)
            colCellID = Column(np.zeros(len(output)), name="CellID #", dtype=int)

            # adds to table
            output.add_column(colBarcode, index=0)
            output.add_column(colROI, index=0)
            output.add_column(colBuid, index=0)
            output.add_column(colCellID, index=2)

            # changes format of table
            # for col in output.colnames:
            #    output[col].info.format = '%.8g'  # for consistent table output

        elif label == "DAPI" and rootFileName.split("_")[2] == "DAPI":
            if param.param["segmentedObjects"]["background_method"] == "flat":
                output = segmentMaskInhomogBackground(im, param)
            elif param.param["segmentedObjects"]["background_method"] == "inhomogeneous":
                output = segmentMaskInhomogBackground(im, param)
            elif param.param["segmentedObjects"]["background_method"] == "stardist":
                output, labeled = segmentMaskStardist(im, param)
            else:
                log1.report(
                    "segmentedObjects/background_method not specified in json file", "ERROR",
                )
                output = np.zeros(1)
                return output

            # show results
            if "labeled" in locals():
                outputFileNameStarDist = (
                    dataFolder.outputFolders["segmentedObjects"] + os.sep + rootFileName + "_stardist"
                )
                showsImageMasks(im, log1, labeled, outputFileNameStarDist)
            showsImageMasks(im, log1, output, outputFileName)

            # saves output 2d zProjection as matrix
            Im.saveImage2D(log1, dataFolder.outputFolders["zProject"])
            saveImage2Dcmd(output, outputFileName + "_Masks", log1)
        else:
            output = []
        del Im

        return output
    else:
        log1.report(
            "2D aligned file does not exist:{}\n{}\n{}\n{}".format(
                fileName_2d_aligned,
                fileName in session1.data.keys(),
                param.param["segmentedObjects"]["operation"] == "overwrite",
                os.path.exists(fileName_2d_aligned),
            ),
            "Error",
        )
        return []


def segmentMasks(param, log1, session1,fileName=None):
    sessionName = "segmentMasks"

    # processes folders and files
    dataFolder = folders(param.param["rootFolder"])
    log1.addSimpleText(
        "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
    )
    log1.report("folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        log1.fileNameMD, "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
    )
    barcodesCoordinates = Table()

    for currentFolder in dataFolder.listFolders:
        # currentFolder=dataFolder.listFolders[0]
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        log1.report("-------> Processing Folder: {}".format(currentFolder))
        log1.report("About to read {} files\n".format(len(param.fileList2Process)))

        for fileName2Process in param.fileList2Process:
            if fileName==None or (fileName!=None and os.path.basename(fileName)==os.path.basename(fileName2Process)):
                label = param.param["acquisition"]["label"]
                if label != "fiducial":
                    output = makesSegmentations(fileName2Process, param, log1, session1, dataFolder)
                    if label == "barcode":
                        outputFile = dataFolder.outputFiles["segmentedObjects"] + "_" + label + ".dat"
                        barcodesCoordinates = vstack([barcodesCoordinates, output])
                        barcodesCoordinates.write(outputFile, format="ascii.ecsv", overwrite=True)
                        log1.report("File {} written to file.".format(outputFile), "info")
                    session1.add(fileName2Process, sessionName)


# # =============================================================================
# # MAIN
# # =============================================================================
# if __name__ == "__main__":
#     begin_time = datetime.now()

#     parser = argparse.ArgumentParser()
#     parser.add_argument("-F", "--rootFolder", help="Folder with images")
#     parser.add_argument("-x", "--fileName", help="fileName to analyze")

#     args = parser.parse_args()

#     print("\n--------------------------------------------------------------------------")

#     if args.rootFolder:
#         rootFolder = args.rootFolder
#     else:
#         rootFolder = os.getcwd()
    
#     if args.fileName:
#         fileName = args.fileName
#     else:
#         fileName = None
        
#     print("parameters> rootFolder: {}".format(rootFolder))
#     now = datetime.now()

#     labels2Process = [
#         {"label": "fiducial", "parameterFile": "infoList_fiducial.json"},
#         {"label": "barcode", "parameterFile": "infoList_barcode.json"},
#         {"label": "DAPI", "parameterFile": "infoList_DAPI.json"},
#         {"label": "RNA", "parameterFile": "infoList_RNA.json"},
#     ]

#     # session
#     sessionName = "segmentMasks"
#     session1 = session(rootFolder, sessionName)
#     # setup logs
#     log1 = log(rootFolder)
#     log1.addSimpleText("\n^^^^^^^^^^^^^^^^^^^^^^^^^^{}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n".format(sessionName))
#     log1.report("Hi-M analysis MD: {}".format(log1.fileNameMD))
#     writeString2File(
#         log1.fileNameMD, "# Hi-M analysis {}".format(now.strftime("%Y/%m/%d %H:%M:%S")), "w",
#     )  # initialises MD file

#     for ilabel in range(len(labels2Process)):
#         label = labels2Process[ilabel]["label"]
#         labelParameterFile = labels2Process[ilabel]["parameterFile"]
#         log1.addSimpleText("**Analyzing label: {}**".format(label))

#         # sets parameters
#         param = Parameters(rootFolder, labelParameterFile)

#         # [applies registration to DAPI and barcodes]
#         if label != "fiducial" and param.param["acquisition"]["label"] != "fiducial":
#             # [segments DAPI and spot masks]
#             if label != "RNA" and param.param["acquisition"]["label"] != "RNA":
#                 segmentMasks(param, log1, session1,fileName)

#         print("\n")
#         del param
#     # exits
#     session1.save(log1)
#     log1.addSimpleText("\n===================={}====================\n".format("Normal termination"))

#     del log1, session1
#     print("Elapsed time: {}".format(datetime.now() - begin_time))

# # 