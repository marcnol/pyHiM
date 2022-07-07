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

import glob, os, time
import matplotlib.pylab as plt
from matplotlib.path import Path
from scipy.ndimage import gaussian_filter
from scipy.spatial import Voronoi
from skimage.measure import regionprops

import numpy as np
import uuid
from dask.distributed import get_client

from astropy.stats import sigma_clipped_stats, SigmaClip, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, vstack, Column
from photutils import DAOStarFinder, detect_sources
from photutils import detect_threshold, deblend_sources
from photutils import Background2D, MedianBackground
from photutils.segmentation.core import SegmentationImage

from imageProcessing.imageProcessing import Image, saveImage2Dcmd
from fileProcessing.fileManagement import folders, writeString2File, try_get_client, printLog

# ---- stardist
import matplotlib

matplotlib.rcParams["image.interpolation"] = None

# to remove in a future version
import warnings

warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================


def _showsImageSources(im, im1_bkg_substracted, x, y, flux, percent=99.5, vmin=0, vmax=2000):

    fig, ax = plt.subplots()
    fig.set_size_inches((50, 50))

    norm = simple_norm(im, "sqrt", percent=percent)
    ax.imshow(im1_bkg_substracted, cmap="Greys", origin="lower", norm=norm)

    p1 = ax.scatter(x, y, c=flux, s=50, facecolors="none", cmap="jet", marker="x", vmin=vmin, vmax=vmax)
    fig.colorbar(p1, ax=ax, fraction=0.046, pad=0.04)

    ax.set_xlim(0, im.shape[1] - 1)
    ax.set_ylim(0, im.shape[0] - 1)

    return fig


def showsImageSources(im, im1_bkg_substracted, sources, fileNameMD,outputFileName):

    percent = 99.5
    flux = sources["flux"]
    x = sources["xcentroid"] + 0.5
    y = sources["ycentroid"] + 0.5

    fig = _showsImageSources(im, im1_bkg_substracted, x, y, flux, percent=percent)
    fig.savefig(outputFileName + "_segmentedSources.png")
    plt.close(fig)

    writeString2File(
        fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_segmentedSources.png"),
        "a",
    )


def showsImageMasks(im, segm_deblend, fileNameMD, outputFileName):
    from stardist import random_label_cmap
    lbl_cmap = random_label_cmap()

    norm = ImageNormalize(stretch=SqrtStretch())
    cmap = lbl_cmap

    fig = plt.figure()
    fig.set_size_inches((30, 30))
    plt.imshow(im, cmap="Greys_r", origin="lower", norm=norm)
    plt.imshow(segm_deblend, origin="lower", cmap=cmap, alpha=0.5)
    plt.savefig(outputFileName + "_segmentedMasks.png")
    plt.close()
    writeString2File(
        fileNameMD,
        "{}\n ![]({})\n".format(os.path.basename(outputFileName), outputFileName + "_segmentedMasks.png"),
        "a",
    )


def _segmentSourceInhomogBackground(im, threshold_over_std, fwhm, brightest, sigma_clip):
    """
    Function that segments barcodes by estimating inhomogeneous background
    Parameters
    ----------
    im : NPY 2D
        image to be segmented

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

    # estimates and removes inhomogeneous background
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)

    im1_bkg_substracted = im - bkg.background
    mean, median, std = sigma_clipped_stats(im1_bkg_substracted, sigma=3.0)

    # estimates sources
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_over_std * std, brightest=brightest, exclude_border=True,)
    sources = daofind(im1_bkg_substracted)

    return sources, im1_bkg_substracted


def segmentSourceInhomogBackground(im, param):
    """
    Wrapper for function that segments barcodes by estimating inhomogeneous background
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

    # sigma_clip = SigmaClip(sigma=3.0)
    sigma_clip = SigmaClip(sigma=param.param["segmentedObjects"]["background_sigma"])

    sources, im1_bkg_substracted = _segmentSourceInhomogBackground(im, threshold_over_std, fwhm, brightest, sigma_clip)
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


def tessellate_masks(segm_deblend):
    """
    * takes a labeled mask (background 0, nuclei labeled 1, 2, ...)
    * calls get_tessellation(xy, img_shape)
    * returns the tesselated mask and the voronoi data structure

    Parameters
    ----------
    segm_deblend : TYPE
        Labeled mask.

    Returns
    -------
    voronoiData : TYPE
        DESCRIPTION.
    mask_voronoi : TYPE
        DESCRIPTION.

    """
    start_time = time.time()

    # get centroids
    mask_labeled = segm_deblend.data
    mask_binary = mask_labeled.copy()
    mask_binary[mask_binary > 0] = 1

    regions = regionprops(mask_labeled)

    numMasks = np.max(mask_labeled)
    centroid = np.zeros((numMasks + 1, 2))  # +1 as labels run from 0 to max

    for props in regions:
        y0, x0 = props.centroid
        label = props.label
        centroid[label, :] = x0, y0

    # tesselation
    # remove first centroid (this is the background label)
    xy = centroid[1:, :]
    voronoiData = get_tessellation(xy, mask_labeled.shape)

    # add some clipping to the tessellation
    # gaussian blur and thresholding; magic numbers!
    mask_blurred = gaussian_filter(mask_binary.astype("float64"), sigma=20)
    mask_blurred = mask_blurred > 0.01

    # convert tessellation to mask
    np.random.seed(42)

    mask_voronoi = np.zeros(mask_labeled.shape, dtype="int64")

    nx, ny = mask_labeled.shape

    # Create vertex coordinates for each grid cell...
    # (<0,0> is at the top left of the grid in this system)
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x, y = x.flatten(), y.flatten()

    points = np.vstack((x, y)).T

    # currently takes 1min for approx 600 polygons
    for label in range(0, numMasks):  # label is shifted by -1 now
        maskID = label + 1

        idx_vor_region = voronoiData.point_region[label]
        idx_vor_vertices = voronoiData.regions[idx_vor_region]  # list of indices of the Voronoi vertices

        vertices = np.full((len(idx_vor_vertices), 2), np.NaN)
        drop_vert = False
        for i in range(len(idx_vor_vertices)):
            idx = idx_vor_vertices[i]
            if idx == -1:  # this means a "virtual point" at infinity as the vertex is not closed
                drop_vert = True
                printLog('$ Detected "virtual point" at infinity. Skipping this mask.')
                break
            vertices[i, :] = voronoiData.vertices[idx]

        if drop_vert:  # region is not bounded
            continue

        poly_path = Path(vertices)
        mask = poly_path.contains_points(points)
        mask = mask.reshape((ny, nx))
        mask_voronoi[mask & mask_blurred] = maskID

    # printLog("--- Took {:.2f}s seconds ---".format(time.time() - start_time))
    printLog("$ Tessellation took {:.2f}s seconds.".format(time.time() - start_time))

    return voronoiData, mask_voronoi


def get_tessellation(xy, img_shape):
    """
    * runs the actual tesselation based on the xy position of the markers in
    an image of given shape

    # follow this tutorial
    # https://hpaulkeeler.com/voronoi-dirichlet-tessellations/
    # https://github.com/hpaulkeeler/posts/blob/master/PoissonVoronoi/PoissonVoronoi.py

    # changes:
    # added dummy points outside of the image corners (in quite some distance)
    # they are supposed "catch" all the vertices that end up at infinity
    # follows an answer given here
    # https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647

    # Attributes
    #    points ndarray of double, shape (npoints, ndim)
    #        Coordinates of input points.
    #
    #    vertices ndarray of double, shape (nvertices, ndim)
    #        Coordinates of the Voronoi vertices.
    #
    #    ridge_points ndarray of ints, shape (nridges, 2)
    #        Indices of the points between which each Voronoi ridge lies.
    #
    #    ridge_vertices list of list of ints, shape (nridges, *)
    #        Indices of the Voronoi vertices forming each Voronoi ridge.
    #
    #    regions list of list of ints, shape (nregions, *)
    #        Indices of the Voronoi vertices forming each Voronoi region. -1 indicates vertex outside the Voronoi diagram.
    #
    #    point_region list of ints, shape (npoints)
    #        Index of the Voronoi region for each input point. If qhull option “Qc” was not specified, the list will contain -1 for points that are not associated with a Voronoi region.
    #
    #    furthest_site
    #        True if this was a furthest site triangulation and False if not.
    #        New in version 1.4.0.

    Parameters
    ----------
    xy : TYPE
        DESCRIPTION.
    img_shape : TYPE
        DESCRIPTION.

    Returns
    -------
    voronoiData : TYPE
        DESCRIPTION.


    """

    x_center, y_center = np.array(img_shape) / 2
    x_max, y_max = np.array(img_shape)

    corner1 = [x_center - 100 * x_max, y_center - 100 * y_max]
    corner2 = [x_center + 100 * x_max, y_center - 100 * y_max]
    corner3 = [x_center - 100 * x_max, y_center + 100 * y_max]
    corner4 = [x_center + 100 * x_max, y_center + 100 * y_max]

    xy = np.append(xy, [corner1, corner2, corner3, corner4], axis=0)

    # perform Voroin tesseslation
    voronoiData = Voronoi(xy)

    return voronoiData


def segmentMaskInhomogBackground(im, param):
    """
    Function used for segmenting masks with the ASTROPY library that uses image processing

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
        npixels=param.param["segmentedObjects"]["area_min"],  # typically 50 for masks
        filter_kernel=kernel,
        nlevels=32,
        contrast=0.001,  # try 0.2 or 0.3
        relabel=True,
    )

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        # printLog('label {}, with area {}'.format(label,area))
        if area < param.param["segmentedObjects"]["area_min"] or area > param.param["segmentedObjects"]["area_max"]:
            segm_deblend.remove_label(label=label)
            # printLog('label {} removed'.format(label))

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend


def segmentMaskStardist(im, param):
    """
    Function used for segmenting masks with the STARDIST package that uses Deep Convolutional Networks

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

    from csbdeep.utils import normalize
    from stardist.models import StarDist2D
    np.random.seed(6)

    # removes background
    # threshold = detect_threshold(im, nsigma=2.0)
    # sigma_clip = SigmaClip(sigma=param.param["segmentedObjects"]["background_sigma"])

    # bkg_estimator = MedianBackground()
    # bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
    # threshold = bkg.background + (
    #     param.param["segmentedObjects"]["threshold_over_std"] * bkg.background_rms
    # )  # background-only error image, typically 1.0

    sigma = param.param["segmentedObjects"]["fwhm"] * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    n_channel = 1 if im.ndim == 2 else im.shape[-1]
    axis_norm = (0, 1)  # normalize channels independently

    if n_channel > 1:
        printLog(
            "> Normalizing image channels %s." % ("jointly" if axis_norm is None or 2 in axis_norm else "independently")
        )

    model = StarDist2D(
        None,
        name=param.param["segmentedObjects"]["stardist_network"],
        basedir=param.param["segmentedObjects"]["stardist_basename"],
    )

    img = normalize(im, 1, 99.8, axis=axis_norm)
    labeled, details = model.predict_instances(img)

    # estimates masks and deblends
    segm = SegmentationImage(labeled)

    # removes masks too close to border
    segm.remove_border_labels(border_width=10)  # parameter to add to infoList
    segm_deblend = segm

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        if area < param.param["segmentedObjects"]["area_min"] or area > param.param["segmentedObjects"]["area_max"]:
            segm_deblend.remove_label(label=label)

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend, labeled


def makesSegmentations(fileName, param, session1, dataFolder):

    rootFileName = os.path.basename(fileName).split(".")[0]
    outputFileName = dataFolder.outputFolders["segmentedObjects"] + os.sep + rootFileName
    fileName_2d_aligned = dataFolder.outputFolders["alignImages"] + os.sep + rootFileName + "_2d_registered.npy"

    printLog("> searching for {}".format(fileName_2d_aligned))
    if os.path.exists(fileName_2d_aligned):  # file exists

        ROI = os.path.basename(fileName).split("_")[param.param["acquisition"]["positionROIinformation"]]
        label = param.param["acquisition"]["label"]

        # loading registered 2D projection
        Im = Image(param)
        Im.loadImage2D(
            fileName, dataFolder.outputFolders["alignImages"], tag="_2d_registered",
        )
        im = Im.data_2D
        printLog("> [{}] Loaded 2D registered file: {}".format(label, os.path.basename(fileName)))

        ##########################################
        #               Segments barcodes
        ##########################################

        if label == "barcode" and len([i for i in rootFileName.split("_") if "RT" in i]) > 0:
            segmentationMathod = param.param["segmentedObjects"]["background_method"]
            printLog(f"\n$ Segmenting barcodes using method: {segmentationMathod }")
            if  segmentationMathod == "flat":
                output = segmentSourceFlatBackground(im, param)
            elif segmentationMathod == "inhomogeneous":
                output, im1_bkg_substracted = segmentSourceInhomogBackground(im, param)
            else:
                printLog(
                    f"# Method <{segmentationMathod}> not available for barcode segmentation!")
                return Table()

            # show results
            showsImageSources(im, im1_bkg_substracted, output, param.param["fileNameMD"], outputFileName)

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
            zcoord = Column(np.nan*np.zeros(len(output)), name="zcentroid", dtype=float)

            # adds to table
            output.add_column(colBarcode, index=0)
            output.add_column(colROI, index=0)
            output.add_column(colBuid, index=0)
            output.add_column(colCellID, index=2)
            output.add_column(zcoord, index=5)

            # changes format of table
            # for col in output.colnames:
            #    output[col].info.format = '%.8g'  # for consistent table output

        #######################################
        #           Segments Masks
        #######################################
        elif (label == "DAPI" or label == "mask"): # and rootFileName.split("_")[2] == "DAPI":
            if param.param["segmentedObjects"]["background_method"] == "flat":
                output = segmentMaskInhomogBackground(im, param)
            elif param.param["segmentedObjects"]["background_method"] == "inhomogeneous":
                output = segmentMaskInhomogBackground(im, param)
            elif param.param["segmentedObjects"]["background_method"] == "stardist":
                output, labeled = segmentMaskStardist(im, param)
            else:
                printLog(
                    "# segmentedObjects/background_method not specified in json file")
                output = np.zeros(1)
                return output

            if "tesselation" in param.param["segmentedObjects"].keys():
                if param.param["segmentedObjects"]["tesselation"]:
                    _, output = tessellate_masks(output)

            # show results
            if "labeled" in locals():
                outputFileNameStarDist = (
                    dataFolder.outputFolders["segmentedObjects"] + os.sep + rootFileName + "_stardist"
                )
                showsImageMasks(im, labeled, param.param["fileNameMD"], outputFileNameStarDist)

            showsImageMasks(im, output, param.param["fileNameMD"], outputFileName)

            # saves output 2d zProjection as matrix
            Im.saveImage2D(dataFolder.outputFolders["zProject"])
            saveImage2Dcmd(output, outputFileName + "_Masks")
        else:
            output = []
        del Im

        return output
    else:
        printLog(
            "# 2D aligned file does not exist:{}\n{}\n{}".format(
                fileName_2d_aligned,
                fileName in session1.data.keys(),
                os.path.exists(fileName_2d_aligned),
            ))
        return []


def segmentMasks(param, session1, fileName=None):
    sessionName = "segmentMasks"

    # processes folders and files
    printLog(
        "\n===================={}:{}====================\n".format(sessionName, param.param["acquisition"]["label"])
    )
    dataFolder = folders(param.param["rootFolder"])
    printLog("> folders read: {}".format(len(dataFolder.listFolders)))
    writeString2File(
        param.param["fileNameMD"], "## {}: {}\n".format(sessionName, param.param["acquisition"]["label"]), "a",
    )
    barcodesCoordinates = Table()

    for currentFolder in dataFolder.listFolders:
        # currentFolder=dataFolder.listFolders[0]
        filesFolder = glob.glob(currentFolder + os.sep + "*.tif")
        dataFolder.createsFolders(currentFolder, param)

        # generates lists of files to process
        param.files2Process(filesFolder)
        printLog("> Processing Folder: {}".format(currentFolder))
        printLog("> Files to Segment: {} \n".format(len(param.fileList2Process)))

        # fileName2ProcessList = [x for x in param.fileList2Process if fileName==None or (fileName!=None and os.path.basename(fileName)==os.path.basename(x))]
        label = param.param["acquisition"]["label"]
        outputFile = dataFolder.outputFiles["segmentedObjects"] + "_" + label + ".dat"

        if param.param["parallel"]:
            # running in parallel mode
            client = get_client()
            futures = list()

            for fileName2Process in param.fileList2Process:
                if fileName == None or (
                    fileName != None and os.path.basename(fileName) == os.path.basename(fileName2Process)
                ):
                    if label != "fiducial":
                        futures.append(
                            client.submit(makesSegmentations, fileName2Process, param, session1, dataFolder)
                        )
                        session1.add(fileName2Process, sessionName)

            printLog("Waiting for {} results to arrive".format(len(futures)))

            results = client.gather(futures)

            if label == "barcode":
                # gathers results from different barcodes and ROIs
                printLog("Retrieving {} results from cluster".format(len(results)))
                detectedSpots = []
                for result in results:
                    detectedSpots.append(len(result))
                    barcodesCoordinates = vstack([barcodesCoordinates, result])

                    # saves results together into a single Table
                    barcodesCoordinates.write(outputFile, format="ascii.ecsv", overwrite=True)
                printLog("$ File {} written to file.".format(outputFile))
                printLog("$ Detected spots: {}".format(",".join([str(x) for x in detectedSpots])))

        else:

            for fileName2Process in param.fileList2Process:
                if fileName == None or (
                    fileName != None and os.path.basename(fileName) == os.path.basename(fileName2Process)
                ):
                    if label != "fiducial":

                        # running in sequential mode
                        output = makesSegmentations(fileName2Process, param, session1, dataFolder)

                        # gathers results from different barcodes and ROIs
                        if label == "barcode":
                            barcodesCoordinates = vstack([barcodesCoordinates, output])
                            barcodesCoordinates.write(outputFile, format="ascii.ecsv", overwrite=True)
                            printLog("$ File {} written to file.".format(outputFile))

                        session1.add(fileName2Process, sessionName)

    return 0
