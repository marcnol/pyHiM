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
from __future__ import absolute_import, division, print_function, unicode_literals

import glob
import os
import time
import uuid

# to remove in a future version
import warnings
# ---- stardist
import matplotlib
import matplotlib.pylab as plt
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import SigmaClip, gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy.table import Column, Table, vstack
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from dask.distributed import get_client
from matplotlib.path import Path
from photutils import (
    Background2D,
    DAOStarFinder,
    MedianBackground,
    deblend_sources,
    detect_sources,
    detect_threshold,
)
from csbdeep.utils import normalize
from stardist.models import StarDist2D
from photutils.segmentation.core import SegmentationImage
from scipy.ndimage import gaussian_filter
from scipy.spatial import Voronoi
from skimage.measure import regionprops
from stardist import random_label_cmap
from fileProcessing.fileManagement import (
    Folders,
    print_log,
    write_string_to_file,
)
from imageProcessing.imageProcessing import Image, save_image_2d_cmd

matplotlib.rcParams["image.interpolation"] = None


warnings.filterwarnings("ignore")

# =============================================================================
# FUNCTIONS
# =============================================================================


def _show_image_sources(
    im, im1_bkg_substracted, x, y, flux, percent=99.5, vmin=0, vmax=2000
):

    fig, ax = plt.subplots()
    fig.set_size_inches((50, 50))

    norm = simple_norm(im, "sqrt", percent=percent)
    ax.imshow(im1_bkg_substracted, cmap="Greys", origin="lower", norm=norm)

    p_1 = ax.scatter(
        x,
        y,
        c=flux,
        s=50,
        facecolors="none",
        cmap="jet",
        marker="x",
        vmin=vmin,
        vmax=vmax,
    )
    fig.colorbar(p_1, ax=ax, fraction=0.046, pad=0.04)

    ax.set_xlim(0, im.shape[1] - 1)
    ax.set_ylim(0, im.shape[0] - 1)

    return fig


def show_image_sources(
    im, im1_bkg_substracted, sources, markdown_filename, output_filename
):

    percent = 99.5
    flux = sources["flux"]
    x = sources["xcentroid"] + 0.5
    y = sources["ycentroid"] + 0.5

    fig = _show_image_sources(im, im1_bkg_substracted, x, y, flux, percent=percent)
    fig.savefig(output_filename + "_segmentedSources.png")
    plt.close(fig)

    write_string_to_file(
        markdown_filename,
        "{}\n ![]({})\n".format(
            os.path.basename(output_filename), output_filename + "_segmentedSources.png"
        ),
        "a",
    )


def show_image_masks(im, segm_deblend, markdown_filename, output_filename):

    lbl_cmap = random_label_cmap()

    norm = ImageNormalize(stretch=SqrtStretch())
    cmap = lbl_cmap

    fig = plt.figure()
    fig.set_size_inches((30, 30))
    plt.imshow(im, cmap="Greys_r", origin="lower", norm=norm)
    plt.imshow(segm_deblend, origin="lower", cmap=cmap, alpha=0.5)
    plt.savefig(output_filename + "_segmentedMasks.png")
    plt.close()
    write_string_to_file(
        markdown_filename,
        "{}\n ![]({})\n".format(
            os.path.basename(output_filename), output_filename + "_segmentedMasks.png"
        ),
        "a",
    )


def _segment_source_inhomog_background(
    im, threshold_over_std, fwhm, brightest, sigma_clip
):
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
    bkg = Background2D(
        im,
        (64, 64),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )

    im1_bkg_substracted = im - bkg.background
    _, _, std = sigma_clipped_stats(im1_bkg_substracted, sigma=3.0)

    # estimates sources
    daofind = DAOStarFinder(
        fwhm=fwhm,
        threshold=threshold_over_std * std,
        brightest=brightest,
        exclude_border=True,
    )
    sources = daofind(im1_bkg_substracted)

    return sources, im1_bkg_substracted


def segment_source_inhomog_background(im, current_param):
    """
    Wrapper for function that segments barcodes by estimating inhomogeneous background
    Parameters
    ----------
    im : NPY 2D
        image to be segmented
    current_param : Parameters
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
      derivation matches that of `DAOFIND` if ``sky`` is 0.0.
    * ``mag``: the object instrumental magnitude calculated as
      ``-2.5 * log10(flux)``.  The derivation matches that of
      `DAOFIND` if ``sky`` is 0.0.

    `None` is returned if no stars are found.

    img_bkc_substracted: 2D NPY array with background substracted image
    """

    threshold_over_std = current_param.param_dict["segmentedObjects"][
        "threshold_over_std"
    ]
    fwhm = current_param.param_dict["segmentedObjects"]["fwhm"]
    brightest = current_param.param_dict["segmentedObjects"][
        "brightest"
    ]  # keeps brightest sources

    # sigma_clip = SigmaClip(sigma=3.0)
    sigma_clip = SigmaClip(
        sigma=current_param.param_dict["segmentedObjects"]["background_sigma"]
    )

    sources, im1_bkg_substracted = _segment_source_inhomog_background(
        im, threshold_over_std, fwhm, brightest, sigma_clip
    )
    return sources, im1_bkg_substracted


def segment_source_flat_background(im, current_param):
    """
    Segments barcodes using flat background
    Parameters
    ----------
    im : NPY 2D
        image to be segmented
    current_param : Parameters
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
      derivation matches that of `DAOFIND` if ``sky`` is 0.0.
    * ``mag``: the object instrumental magnitude calculated as
      ``-2.5 * log10(flux)``.  The derivation matches that of
      `DAOFIND` if ``sky`` is 0.0.

    `None` is returned if no stars are found.

    img_bkc_substracted: 2D NPY array with background substracted image
    """

    threshold_over_std = current_param.param_dict["segmentedObjects"][
        "threshold_over_std"
    ]
    fwhm = current_param.param_dict["segmentedObjects"]["fwhm"]

    # removes background
    _, median, std = sigma_clipped_stats(
        im, current_param.param_dict["segmentedObjects"]["background_sigma"]
    )
    im1_bkg_substracted = im - median

    # estimates sources
    daofind = DAOStarFinder(
        fwhm=fwhm, threshold=threshold_over_std * std, exclude_border=True
    )
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
    voronoi_data : TYPE
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

    num_masks = np.max(mask_labeled)
    centroid = np.zeros((num_masks + 1, 2))  # +1 as labels run from 0 to max

    for props in regions:
        y_0, x_0 = props.centroid
        label = props.label
        centroid[label, :] = x_0, y_0

    # tesselation
    # remove first centroid (this is the background label)
    xy = centroid[1:, :]
    voronoi_data = get_tessellation(xy, mask_labeled.shape)

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
    for label in range(0, num_masks):  # label is shifted by -1 now
        mask_id = label + 1

        idx_vor_region = voronoi_data.point_region[label]
        idx_vor_vertices = voronoi_data.regions[
            idx_vor_region
        ]  # list of indices of the Voronoi vertices

        vertices = np.full((len(idx_vor_vertices), 2), np.NaN)
        drop_vert = False
        # pylint: disable-next=consider-using-enumerate
        for i in range(len(idx_vor_vertices)):
            idx = idx_vor_vertices[i]
            if (
                idx == -1
            ):  # this means a "virtual point" at infinity as the vertex is not closed
                drop_vert = True
                print_log('$ Detected "virtual point" at infinity. Skipping this mask.')
                break
            vertices[i, :] = voronoi_data.vertices[idx]

        if drop_vert:  # region is not bounded
            continue

        poly_path = Path(vertices)
        mask = poly_path.contains_points(points)
        mask = mask.reshape((ny, nx))
        mask_voronoi[mask & mask_blurred] = mask_id

    # print_log("--- Took {:.2f}s seconds ---".format(time.time() - start_time))
    print_log("$ Tessellation took {:.2f}s seconds.".format(time.time() - start_time))

    return voronoi_data, mask_voronoi


def get_tessellation(xy, img_shape):
    """
    * runs the actual tesselation based on the xy position of the markers in an image of given shape

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
    #    ridge_vertices list of list of ints, shape (nridges, \*)
    #        Indices of the Voronoi vertices forming each Voronoi ridge.
    #
    #    regions list of list of ints, shape (nregions, \*)
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
    voronoi_data : TYPE
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
    voronoi_data = Voronoi(xy)

    return voronoi_data


def segment_mask_inhomog_background(im, current_param):
    """
    Function used for segmenting masks with the ASTROPY library that uses image processing

    Parameters
    ----------
    im : 2D np array
        image to be segmented.
    current_param : Parameters class
        parameters.

    Returns
    -------
    segm_deblend: 2D np array where each pixel contains the label of the mask segmented. Background: 0

    """
    # removes background
    threshold = detect_threshold(im, nsigma=2.0)
    sigma_clip = SigmaClip(
        sigma=current_param.param_dict["segmentedObjects"]["background_sigma"]
    )

    bkg_estimator = MedianBackground()
    bkg = Background2D(
        im,
        (64, 64),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )
    threshold = bkg.background + (
        current_param.param_dict["segmentedObjects"]["threshold_over_std"]
        * bkg.background_rms
    )  # background-only error image, typically 1.0

    sigma = (
        current_param.param_dict["segmentedObjects"]["fwhm"] * gaussian_fwhm_to_sigma
    )  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    # estimates masks and deblends
    segm = detect_sources(
        im,
        threshold,
        npixels=current_param.param_dict["segmentedObjects"]["area_min"],
        filter_kernel=kernel,
    )

    # removes masks too close to border
    segm.remove_border_labels(border_width=10)  # parameter to add to infoList

    segm_deblend = deblend_sources(
        im,
        segm,
        npixels=current_param.param_dict["segmentedObjects"][
            "area_min"
        ],  # typically 50 for masks
        filter_kernel=kernel,
        nlevels=32,
        contrast=0.001,  # try 0.2 or 0.3
        relabel=True,
    )

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        # print_log('label {}, with area {}'.format(label,area))
        if (
            area < current_param.param_dict["segmentedObjects"]["area_min"]
            or area > current_param.param_dict["segmentedObjects"]["area_max"]
        ):
            segm_deblend.remove_label(label=label)
            # print_log('label {} removed'.format(label))

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend


def segment_mask_stardist(im, current_param):
    """
    Function used for segmenting masks with the STARDIST package that uses Deep Convolutional Networks

    Parameters
    ----------
    im : 2D np array
        image to be segmented.
    current_param : Parameters class
        parameters.

    Returns
    -------
    segm_deblend: 2D np array where each pixel contains the label of the mask segmented. Background: 0

    """


    np.random.seed(6)

    # removes background
    # threshold = detect_threshold(im, nsigma=2.0)
    # sigma_clip = SigmaClip(sigma = current_param.param_dict["segmentedObjects"]["background_sigma"])

    # bkg_estimator = MedianBackground()
    # bkg = Background2D(im, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)
    # threshold = bkg.background + (
    #     current_param.param_dict["segmentedObjects"]["threshold_over_std"] * bkg.background_rms
    # )  # background-only error image, typically 1.0

    sigma = (
        current_param.param_dict["segmentedObjects"]["fwhm"] * gaussian_fwhm_to_sigma
    )  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    n_channel = 1 if im.ndim == 2 else im.shape[-1]
    axis_norm = (0, 1)  # normalize channels independently

    if n_channel > 1:
        print_log(
            "> Normalizing image channels %s."
            % ("jointly" if axis_norm is None or 2 in axis_norm else "independently")
        )

    model = StarDist2D(
        None,
        name=current_param.param_dict["segmentedObjects"]["stardist_network"],
        basedir=current_param.param_dict["segmentedObjects"]["stardist_basename"],
    )

    img = normalize(im, 1, 99.8, axis=axis_norm)
    labeled, _ = model.predict_instances(img)

    # estimates masks and deblends
    segm = SegmentationImage(labeled)

    # removes masks too close to border
    segm.remove_border_labels(border_width=10)  # parameter to add to infoList
    segm_deblend = segm

    # removes Masks too big or too small
    for label in segm_deblend.labels:
        # take regions with large enough areas
        area = segm_deblend.get_area(label)
        if (
            area < current_param.param_dict["segmentedObjects"]["area_min"]
            or area > current_param.param_dict["segmentedObjects"]["area_max"]
        ):
            segm_deblend.remove_label(label=label)

    # relabel so masks numbers are consecutive
    segm_deblend.relabel_consecutive()

    return segm_deblend, labeled


def make_segmentations(file_name, current_param, current_session, data_folder):

    root_filename = os.path.basename(file_name).split(".")[0]
    output_filename = (
        data_folder.output_folders["segmentedObjects"] + os.sep + root_filename
    )
    filename_2d_aligned = (
        data_folder.output_folders["alignImages"]
        + os.sep
        + root_filename
        + "_2d_registered.npy"
    )

    print_log("> searching for {}".format(filename_2d_aligned))
    if os.path.exists(filename_2d_aligned):  # file exists

        roi = os.path.basename(file_name).split("_")[
            current_param.param_dict["acquisition"]["positionROIinformation"]
        ]
        label = current_param.param_dict["acquisition"]["label"]

        # loading registered 2D projection
        im_obj = Image(current_param)
        im_obj.load_image_2d(
            file_name, data_folder.output_folders["alignImages"], tag="_2d_registered",
        )
        im = im_obj.data_2d
        print_log(
            "> [{}] Loaded 2D registered file: {}".format(
                label, os.path.basename(file_name)
            )
        )

        ##########################################
        #               Segments barcodes
        ##########################################

        if (
            label == "barcode"
            and len([i for i in root_filename.split("_") if "RT" in i]) > 0
        ):
            segmentation_method = current_param.param_dict["segmentedObjects"][
                "background_method"
            ]
            print_log(f"\n$ Segmenting barcodes using method: {segmentation_method }")
            if segmentation_method == "flat":
                output = segment_source_flat_background(im, current_param)
            elif segmentation_method == "inhomogeneous":
                output, im1_bkg_substracted = segment_source_inhomog_background(
                    im, current_param
                )
            else:
                print_log(
                    f"# Method <{segmentation_method}> not available for barcode segmentation!"
                )
                return Table()

            # show results
            show_image_sources(
                im,
                im1_bkg_substracted,
                output,
                current_param.param_dict["fileNameMD"],
                output_filename,
            )

            # [ formats results Table for output by adding buid, barcode_id, CellID and roi]

            # buid
            buid = []
            for i in range(len(output)):
                buid.append(str(uuid.uuid4()))
            col_buid = Column(buid, name="Buid", dtype=str)

            # barcode_id, cellID and roi
            barcode_id = os.path.basename(file_name).split("_")[2].split("RT")[1]
            col_roi = Column(int(roi) * np.ones(len(output)), name="ROI #", dtype=int)
            col_barcode = Column(
                int(barcode_id) * np.ones(len(output)), name="Barcode #", dtype=int
            )
            col_cell_id = Column(np.zeros(len(output)), name="CellID #", dtype=int)
            zcoord = Column(
                np.nan * np.zeros(len(output)), name="zcentroid", dtype=float
            )

            # adds to table
            output.add_column(col_barcode, index=0)
            output.add_column(col_roi, index=0)
            output.add_column(col_buid, index=0)
            output.add_column(col_cell_id, index=2)
            output.add_column(zcoord, index=5)

            # changes format of table
            # for col in output.colnames:
            #    output[col].info.format = '%.8g'  # for consistent table output

        #######################################
        #           Segments Masks
        #######################################
        elif (
            label in ("DAPI", "mask")
        ):  # and root_filename.split("_")[2] == "DAPI":
            if (
                current_param.param_dict["segmentedObjects"]["background_method"]
                == "flat"
            ):
                output = segment_mask_inhomog_background(im, current_param)
            elif (
                current_param.param_dict["segmentedObjects"]["background_method"]
                == "inhomogeneous"
            ):
                output = segment_mask_inhomog_background(im, current_param)
            elif (
                current_param.param_dict["segmentedObjects"]["background_method"]
                == "stardist"
            ):
                output, labeled = segment_mask_stardist(im, current_param)
            else:
                print_log(
                    "# segmentedObjects/background_method not specified in json file"
                )
                output = np.zeros(1)
                return output

            if "tesselation" in current_param.param_dict["segmentedObjects"].keys():
                if current_param.param_dict["segmentedObjects"]["tesselation"]:
                    _, output = tessellate_masks(output)

            # show results
            if "labeled" in locals():
                output_filename_stardist = (
                    data_folder.output_folders["segmentedObjects"]
                    + os.sep
                    + root_filename
                    + "_stardist"
                )
                show_image_masks(
                    im,
                    labeled,
                    current_param.param_dict["fileNameMD"],
                    output_filename_stardist,
                )

            show_image_masks(
                im,
                output,
                current_param.param_dict["fileNameMD"],
                output_filename,
            )

            # saves output 2d zProjection as matrix
            im_obj.save_image_2d(data_folder.output_folders["zProject"])
            save_image_2d_cmd(output, output_filename + "_Masks")
        else:
            output = []
        del im_obj

        return output
    else:
        print_log(
            "# 2D aligned file does not exist:{}\n{}\n{}".format(
                filename_2d_aligned,
                file_name in current_session.data.keys(),
                os.path.exists(filename_2d_aligned),
            )
        )
        return []


def segment_masks(current_param, current_session, file_name=None):
    session_name = "segmentMasks"

    # processes folders and files
    print_log(
        "\n===================={}:{}====================\n".format(
            session_name, current_param.param_dict["acquisition"]["label"]
        )
    )
    data_folder = Folders(current_param.param_dict["rootFolder"])
    print_log("> folders read: {}".format(len(data_folder.list_folders)))
    write_string_to_file(
        current_param.param_dict["fileNameMD"],
        "## {}: {}\n".format(
            session_name, current_param.param_dict["acquisition"]["label"]
        ),
        "a",
    )
    barcodes_coordinates = Table()

    for current_folder in data_folder.list_folders:
        # current_folder=data_folder.list_folders[0]
        files_folder = glob.glob(current_folder + os.sep + "*.tif")
        data_folder.create_folders(current_folder, current_param)

        # generates lists of files to process
        current_param.find_files_to_process(files_folder)
        print_log("> Processing Folder: {}".format(current_folder))
        print_log(
            "> Files to Segment: {} \n".format(len(current_param.files_to_process))
        )

        label = current_param.param_dict["acquisition"]["label"]
        output_file = (
            data_folder.output_files["segmentedObjects"] + "_" + label + ".dat"
        )

        if current_param.param_dict["parallel"]:
            # running in parallel mode
            client = get_client()
            futures = []

            for filename_to_process in current_param.files_to_process:
                if file_name is None or (
                    file_name is not None
                    and os.path.basename(file_name)
                    == os.path.basename(filename_to_process)
                ):
                    if label != "fiducial":
                        futures.append(
                            client.submit(
                                make_segmentations,
                                filename_to_process,
                                current_param,
                                current_session,
                                data_folder,
                            )
                        )
                        current_session.add(filename_to_process, session_name)

            print_log("Waiting for {} results to arrive".format(len(futures)))

            results = client.gather(futures)

            if label == "barcode":
                # gathers results from different barcodes and rois
                print_log("Retrieving {} results from cluster".format(len(results)))
                detected_spots = []
                for result in results:
                    detected_spots.append(len(result))
                    barcodes_coordinates = vstack([barcodes_coordinates, result])

                    # saves results together into a single Table
                    barcodes_coordinates.write(
                        output_file, format="ascii.ecsv", overwrite=True
                    )
                print_log("$ File {} written to file.".format(output_file))
                print_log(
                    "$ Detected spots: {}".format(
                        ",".join([str(x) for x in detected_spots])
                    )
                )

        else:

            for filename_to_process in current_param.files_to_process:
                if file_name is None or (
                    file_name is not None
                    and os.path.basename(file_name)
                    == os.path.basename(filename_to_process)
                ):
                    if label != "fiducial":

                        # running in sequential mode
                        output = make_segmentations(
                            filename_to_process,
                            current_param,
                            current_session,
                            data_folder,
                        )

                        # gathers results from different barcodes and rois
                        if label == "barcode":
                            barcodes_coordinates = vstack([barcodes_coordinates, output])
                            barcodes_coordinates.write(
                                output_file, format="ascii.ecsv", overwrite=True
                            )
                            print_log("$ File {} written to file.".format(output_file))

                        current_session.add(filename_to_process, session_name)

    return 0
