#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:00:52 2020

@author: marcnol

Classes and functions for common image processing
"""
# =============================================================================
# IMPORTS
# =============================================================================

import os
import sys
import warnings

import cv2
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo
from astropy.convolution import Gaussian2DKernel
from astropy.stats import SigmaClip
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from csbdeep.data import PadAndCropResizer
from csbdeep.utils import normalize
from csbdeep.utils.tf import limit_gpu_memory
from matplotlib import ticker
from numpy import linalg as LA
from photutils import Background2D, MedianBackground, deblend_sources, detect_sources
from scipy import ndimage as ndi
from scipy.ndimage import shift as shift_image
from scipy.stats import sigmaclip
from skimage import color, exposure, io, measure
from skimage.feature import peak_local_max
from skimage.metrics import mean_squared_error, normalized_root_mse
from skimage.metrics import structural_similarity as ssim
from skimage.registration import phase_cross_correlation
from skimage.segmentation import watershed
from skimage.util.apply_parallel import apply_parallel
from skimage.util.shape import view_as_blocks
from stardist.models import StarDist3D
from tifffile import imsave
from tqdm import tqdm, trange
from stardist import random_label_cmap

from fileProcessing.fileManagement import print_log, try_get_client

warnings.filterwarnings("ignore")

np.seterr(divide="ignore", invalid="ignore")

# =============================================================================
# CLASSES
# =============================================================================


class Image:
    def __init__(self, param=None, current_log=None):
        self.current_param = param
        self.log = current_log
        self.data = []
        self.file_name = ""
        self.data_2d = np.zeros((1, 1))
        self.stage_coordinates = [0.0, 0.0]
        self.image_size = -1
        self.focus_plane = -1
        self.extension = ""
        self.z_range = None
        self.focal_plane_matrix = None

    # read an image as a numpy array
    def load_image(self, file_name):
        self.data = io.imread(file_name).squeeze()
        self.file_name = file_name
        self.image_size = self.data.shape
        self.extension = file_name.split(".")[-1]

    # save 2D projection as numpy array
    def save_image_2d(self, root_folder, tag="_2d"):
        file_name = self.get_image_filename(root_folder, tag)
        save_image_2d_cmd(self.data_2d, file_name)

    def get_image_filename(self, root_folder, tag):
        file_name = (
            root_folder + os.sep + os.path.basename(self.file_name).split(".")[0] + tag
        )
        return file_name

    # read an image as a numpy array
    def load_image_2d(self, file_name, master_folder, tag="_2d"):
        self.file_name = file_name
        file_name = self.get_image_filename(master_folder, tag) + ".npy"

        self.data_2d = np.load(file_name)
        print_log("$ Loading from disk:{}".format(os.path.basename(file_name)))

    # max intensity projection using all z planes
    def max_intensity_projection(self):
        self.data_2d = np.max(self.data, axis=0)

    # Normalize a 3d image <im> by subtracting local gaussian blur of std <sz>
    def normalize_image(self):
        background = self.data_2d.min()
        _im = self.data_2d - background
        _max = _im.max()
        _im = _im / _max
        return _im

    # returns the picture x,y location, if available
    def get_image_location(self):
        if hasattr(self, "stage_coordinates"):
            return self.stage_coordinates
        else:
            return [0.0, 0.0]

    # returns the film size
    def get_image_size(self):
        return self.image_size

    # returns the film focus
    def get_focus_plane(self):
        if hasattr(self, "focus_plane"):
            return self.focus_plane
        else:
            return 0.0

    # Outputs image properties to command line
    def print_image_properties(self):
        # print_log("Image Name={}".format(self.file_name))
        print_log("$ Image Size={}".format(self.image_size))
        # self.log.report("Stage position={}".format(self.stage_coordinates))
        print_log("$ Focal plane={}".format(self.focus_plane))

    # processes sum image in axial direction given range
    # @jit(nopython=True)
    def z_projection_range(self):

        # find the correct range for the projection
        if self.current_param.param_dict["zProject"]["zmax"] > self.image_size[0]:
            print_log("$ Setting z max to the last plane")
            self.current_param.param_dict["zProject"]["zmax"] = self.image_size[0]

        if self.current_param.param_dict["zProject"]["mode"] == "automatic":
            print_log("> Calculating planes...")
            z_range = calculate_zrange(self.data, self.current_param)

        elif self.current_param.param_dict["zProject"]["mode"] == "full":
            (zmin, zmax) = (0, self.image_size[0])
            z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

        if self.current_param.param_dict["zProject"]["mode"] == "laplacian":
            print_log("Stacking using Laplacian variance...")
            (
                self.data_2d,
                self.focal_plane_matrix,
                self.z_range,
            ) = reinterpolate_focal_plane(self.data, self.current_param.param_dict)
            self.focus_plane = self.z_range[0]

        else:
            # Manual: reads from parameters file
            (zmin, zmax) = (
                self.current_param.param_dict["zProject"]["zmin"],
                self.current_param.param_dict["zProject"]["zmax"],
            )
            if zmin >= zmax:
                raise SystemExit(
                    "zmin is equal or larger than zmax in configuration file. Cannot proceed."
                )
            z_range = (round((zmin + zmax) / 2), range(zmin, zmax))

        if "laplacian" not in self.current_param.param_dict["zProject"]["mode"]:
            self.data_2d = project_image_2d(
                self.data,
                z_range,
                self.current_param.param_dict["zProject"]["zProjectOption"],
            )
            self.focus_plane = z_range[0]
            self.z_range = z_range[1]

        print_log("> Processing z_range:{}".format(self.z_range))

    # displays image and shows it
    def show_image(
        self,
        show=False,
        size=(10, 10),
        output_name="tmp.png",
        save=True,
        normalization="stretch",
    ):
        fig = plt.figure()
        fig.set_size_inches(size)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()

        if normalization == "simple":
            norm = simple_norm(self.data_2d, "sqrt", percent=99.9)
        else:
            norm = ImageNormalize(stretch=SqrtStretch())

        ax.set_title("2D Data")

        if show:
            fig.add_axes(ax)
            ax.imshow(self.data_2d, origin="lower", cmap="Greys_r", norm=norm)
            return ax

        if save and not show:
            fig.add_axes(ax)
            ax.imshow(self.data_2d, origin="lower", cmap="Greys_r", norm=norm)
            fig.savefig(output_name)
            plt.close(fig)

    def remove_background_2d(self, normalize=False):
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(
            self.data_2d,
            (64, 64),
            filter_size=(3, 3),
            sigma_clip=sigma_clip,
            bkg_estimator=bkg_estimator,
        )

        im_bkg_substracted = self.data_2d - bkg.background

        if normalize:
            im_bkg_substracted = (im_bkg_substracted - im_bkg_substracted.min()) / (
                im_bkg_substracted.max()
            )

        return im_bkg_substracted

    def image_show_with_values(self, output_name):
        image_show_with_values(
            [self.focal_plane_matrix],
            title="focal plane = " + "{:.2f}".format(self.focus_plane),
            output_name=output_name,
        )


# =============================================================================
# GENERAL FUNCTIONS
# =============================================================================


def fit_1d_gaussian_scipy(x, y, title="", verbose=False):
    """
    Fits a function using a 1D Gaussian and returns parameters if successfull.
    Otherwise will return an empty dict
    Uses scipy spo package

    Parameters
    ----------
    x : numpy 1D array
        x data.
    y : numpy 1D array
        y data.
    title : str, optional
        figure title. The default is ''.
    verbose : Boolean, optional
        whether fig and results will be shown. The default is True.

    Returns
    -------
    {}
        dictionary with fitting parameters.
    fig
        matplotlib figure for saving

    """
    fit_result = {}

    try:
        fitgauss = spo.curve_fit(gaussian, x, y)
        fit_result["gauss1d.pos"] = fitgauss[0][1]
        fit_result["gauss1d.ampl"] = fitgauss[0][0]
        fit_result["gauss1d.fwhm"] = 2.355 * fitgauss[0][2]
    except RuntimeError:
        return {}, []
    except ValueError:
        fit_result["gauss1d.pos"] = np.mean(x)
        fit_result["gauss1d.ampl"] = 0.0
        fit_result["gauss1d.fwhm"] = 0.0
        print_log("Warning: returned middle plane!")
        return fit_result, []

    if verbose:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        print_log("<<Fitting successful>>")

        ax.plot(x, y, "ko", label="data")
        ax.plot(
            x,
            gaussian(x, fitgauss[0][0], fitgauss[0][1], fitgauss[0][2]),
            linewidth=2,
            label="gaussian fit",
        )
        ax.legend(loc=2)
        ax.set_title(title)
        return fit_result, fig

    return fit_result, []


def make_shift_matrix_hi_res(shift_matrices, block_ref_shape):
    """
    Reinterpolates a block matrix to the full size of a larger image

    Parameters
    ----------
    shift_matrices : list of numpy arrays
        list containing block matrices.
    block_ref_shape : tuple
        shape of block matrix.

    Returns
    -------
    shift_matrix : numpy array
        Reinterpolated (larger) image.
        Size will be N x n x n
        where n is block_ref_shape[3]
        and N is len(shift_matrices)

    """
    number_blocks = block_ref_shape[0]
    block_size_xy = block_ref_shape[3]

    shift_matrix = np.zeros(
        (
            len(shift_matrices),
            block_size_xy * shift_matrices[0].shape[0],
            block_size_xy * shift_matrices[0].shape[1],
        )
    )
    for _ax, m in enumerate(shift_matrices):
        # print_log("size={}".format(m.shape))
        for i in range(number_blocks):
            for j in range(number_blocks):
                shift_matrix[
                    _ax,
                    i * block_size_xy : (i + 1) * block_size_xy,
                    j * block_size_xy : (j + 1) * block_size_xy,
                ] = m[i, j]
    return shift_matrix


def project_image_2d(img, z_range, mode):

    # sums images
    image_size = img.shape
    i_collapsed = np.zeros((image_size[1], image_size[2]))

    if "MIP" in mode:
        # Max projection of selected planes
        i_collapsed = np.max(img[z_range[1][0] : z_range[1][-1]], axis=0)
    elif "sum" in mode:
        # Sums selected planes
        for i in z_range[1]:
            i_collapsed += img[i]
    else:
        print_log(
            "ERROR: mode not recognized. Expected: MIP or sum. Read: {}".format(mode)
        )

    return i_collapsed


# Gaussian function
# @jit(nopython=True)
def gaussian(x, a=1, mean=0, std=0.5):
    return (
        a
        * (1 / (std * (np.sqrt(2 * np.pi))))
        * (np.exp(-((x - mean) ** 2) / ((2 * std) ** 2)))
    )


# Finds best focal plane by determining the max of the std deviation vs z curve

# @jit(nopython=True)
def calculate_zrange(idata, parameters):
    """
    Calculates the focal planes based max standard deviation
    """
    num_planes = (
        parameters.param_dict["zProject"]["zmax"]
        - parameters.param_dict["zProject"]["zmin"]
    )
    std_matrix = np.zeros(num_planes)
    mean_matrix = np.zeros(num_planes)

    # calculate STD in each plane
    for i in range(0, num_planes):
        std_matrix[i] = np.std(idata[i])
        mean_matrix[i] = np.mean(idata[i])

    max_std = np.max(std_matrix)
    i_focus_plane = np.where(std_matrix == max_std)[0][0]
    # Select a window to avoid being on the edges of the stack

    if i_focus_plane < parameters.param_dict["zProject"]["windowSecurity"] or (
        i_focus_plane > num_planes - parameters.param_dict["zProject"]["windowSecurity"]
    ):
        focus_plane = i_focus_plane
    else:
        # interpolate zfocus
        axis_z = range(
            max(
                parameters.param_dict["zProject"]["zmin"],
                i_focus_plane - parameters.param_dict["zProject"]["windowSecurity"],
                min(
                    parameters.param_dict["zProject"]["zmax"],
                    i_focus_plane + parameters.param_dict["zProject"]["windowSecurity"],
                ),
            )
        )

        std_matrix -= np.min(std_matrix)
        std_matrix /= np.max(std_matrix)

        try:
            fitgauss = spo.curve_fit(
                gaussian, axis_z, std_matrix[axis_z[0] : axis_z[-1] + 1]
            )
            # print_log("Estimation of focal plane (px): ", int(fitgauss[0][1]))
            focus_plane = int(fitgauss[0][1])
        except RuntimeError:
            print_log("Warning, too many iterations")
            focus_plane = i_focus_plane

    zmin = max(
        parameters.param_dict["zProject"]["windowSecurity"],
        focus_plane - parameters.param_dict["zProject"]["zwindows"],
    )
    zmax = min(
        num_planes,
        parameters.param_dict["zProject"]["windowSecurity"] + num_planes,
        focus_plane + parameters.param_dict["zProject"]["zwindows"],
    )
    zrange = range(zmin, zmax + 1)

    return focus_plane, zrange


def find_transform(im_src, im_dst):
    warp = np.eye(3, dtype=np.float32)
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 0.001)
    try:
        _, warp = cv2.findTransformECC(
            im_src, im_dst, warp, cv2.MOTION_HOMOGRAPHY, criteria
        )
    except:
        print_log("Warning: find transform failed. Set warp as identity")
    return warp


def variance_of_laplacian(image):
    # compute the Laplacian of the image and then return the focus
    # measure, which is simply the variance of the Laplacian
    return cv2.Laplacian(image, cv2.CV_64F).var()


def scatter_3d_image(image):
    """
    splits 3D image plane by plane and scatteres them to a cluster

    Parameters
    ----------
    image : numpy array
        3D image.

    Returns
    -------
    image_list_scattered :  List, dict, iterator, or queue of futures matching the type of input.
        scattered image.

    """
    number_planes = image.shape[0]
    image_list_scattered = [image[z, :, :] for z in range(number_planes)]
    return image_list_scattered


def reassemble_3d_image(client, futures, output_shape):
    """
    waits for futures to arrive
    collects them into a results list
    reassembles 3D image plane by plane

    Parameters
    ----------
    client : dask CLient()
        result of get_client()
    futures : []
        list of futures
    output_shape : tuple
        result of image.shape

    Returns
    -------
    output : numpy array
        contains reassembled 3D image.

    """
    results = client.gather(futures)
    print_log(" > Retrieving {} results from cluster".format(len(results)))

    output = np.zeros(output_shape)
    for z, result in enumerate(results):
        output[z, :, :] = result

    del results
    return output


# =============================================================================
# CONTRAST and PIXEL INTENSITY NORMALIZATION, INHOMOGENEOUS BACKGROUND
# =============================================================================


def preprocess_3d_image(x, lower_threshold, higher_threshold, parallel_execution=True):
    """
    3D stack pre-procesing:
        - rescales intensities to 0->1
        - removes inhomogeneous background plane by plane
        - adjusts image levels by thresholding

    Parameters
    ----------
    x : numpy array
        3D image.
    lower_threshold : float
        lower threshold for adjusting image levels.
    higher_threshold : float
        higher threshold for adjusting image levels..

    Returns
    -------
    image : numpy array
        pre-processed 3D image.

    """
    image = exposure.rescale_intensity(x, out_range=(0, 1))

    # print_log("Removing inhomogeneous background...")
    image = _remove_inhomogeneous_background(
        image, parallel_execution=parallel_execution
    )

    # print_log("Rescaling grey levels...")
    image = image_adjust(
        image, lower_threshold=lower_threshold, higher_threshold=higher_threshold
    )[0]

    return image


def image_adjust(image, lower_threshold=0.3, higher_threshold=0.9999):
    """
    Adjust intensity levels:
        - rescales exposures
        - gets histogram of pixel intensities to define cutoffs
        - applies thresholds

    Parameters
    ----------
    image : numpy array
        input 3D image.
    lower_threshold : float, optional
        lower threshold for adjusting image levels. The default is 0.3.
    higher_threshold : float, optional
        higher threshold for adjusting image levels.. The default is 0.9999.

    Returns
    -------
    image1 : numpy array
        adjusted 3D image.
    hist1_before : numpy array
        histogram of pixel intensities before adjusting levels.
    hist1 : numpy array
        histogram of pixel intensities after adjusting levels.
    lower_cutoff : float
        lower cutoff used for thresholding.
    higher_cutoff : float
        higher cutoff used for thresholding.

    """
    # print_log("> Rescaling grey levels...")

    # rescales image to [0,1]
    image1 = exposure.rescale_intensity(image, out_range=(0, 1))

    # calculates histogram of intensities
    hist1_before = exposure.histogram(image1)

    sum = np.zeros(len(hist1_before[0]))
    for i in range(len(hist1_before[0]) - 1):
        sum[i + 1] = sum[i] + hist1_before[0][i]

    sum_normalized = sum / sum.max()
    lower_cutoff = np.where(sum_normalized > lower_threshold)[0][0] / 255
    higher_cutoff = np.where(sum_normalized > higher_threshold)[0][0] / 255

    # adjusts image intensities from (lower_threshold,higher_threshold) --> [0,1]
    image1 = exposure.rescale_intensity(
        image1, in_range=(lower_cutoff, higher_cutoff), out_range=(0, 1)
    )

    # calculates histogram of intensities of adjusted image
    hist1 = exposure.histogram(image1)

    return image1, hist1_before, hist1, lower_cutoff, higher_cutoff


def save_image_differences(img_1, img_2, img_3, img_4, output_filename):
    """
    Overlays two images as R and B and saves them to output file
    """

    img_1, img_2 = img_1 / img_1.max(), img_2 / img_2.max()
    img_3, img_4 = img_3 / img_3.max(), img_4 / img_4.max()

    img_1, _, _, _, _ = image_adjust(img_1, lower_threshold=0.5, higher_threshold=0.9999)
    img_2, _, _, _, _ = image_adjust(img_2, lower_threshold=0.5, higher_threshold=0.9999)
    img_3, _, _, _, _ = image_adjust(img_3, lower_threshold=0.5, higher_threshold=0.9999)
    img_4, _, _, _, _ = image_adjust(img_4, lower_threshold=0.5, higher_threshold=0.9999)

    cmap = "seismic"

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((60, 30))

    ax1.imshow(img_1 - img_2, cmap=cmap)
    ax1.axis("off")
    ax1.set_title("uncorrected")

    ax2.imshow(img_3 - img_4, cmap=cmap)
    ax2.axis("off")
    ax2.set_title("corrected")

    fig.savefig(output_filename)

    plt.close(fig)


def _remove_inhomogeneous_background(
    im,
    box_size=(32, 32),
    filter_size=(3, 3),
    verbose=True,
    parallel_execution=True,
    background=False,
):
    """
    wrapper to remove inhomogeneous backgrounds for 2D and 3D images

    Parameters
    ----------
    im : numpy array
        input image.
    box_size : tuple of ints, optional
        size of box_size used for block decomposition. The default is (32, 32).
    filter_size : tuple of ints, optional
        Size of gaussian filter used for smoothing results. The default is (3, 3).
    verbose : boolean, optional
        The default is True.

    Returns
    -------
    output : TYPE
        DESCRIPTION.

    """
    if len(im.shape) == 2:
        output = _remove_inhomogeneous_background_2d(
            im,
            filter_size=filter_size,
            background=background,
        )
    elif len(im.shape) == 3:
        output = _remove_inhomogeneous_background_3d(
            im,
            box_size=box_size,
            filter_size=filter_size,
            parallel_execution=parallel_execution,
            background=background,
        )
    else:
        return None
    return output


def _remove_inhomogeneous_background_2d(im, filter_size=(3, 3), background=False):
    """
    Calls Background2D() from ASTROPY to perform background substraction in a 2D image

    Parameters
    ----------
    im : numpy array
        2D input image.
    filter_size : tuple of ints, optional
        Size of gaussian filter used for smoothing results. The default is (3, 3).
    background : boolean, optional
        if True returs the substracted image and the background
        otherwise only the former. The default is False.

    Returns
    -------
    numpy array
        background substracted 2D image.

    """
    print_log("Removing inhomogeneous background from 2D image...")

    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        im,
        (64, 64),
        filter_size=filter_size,
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )

    im1_bkg_substracted = im - bkg.background

    if background:
        return im1_bkg_substracted, bkg
    else:
        return im1_bkg_substracted


def _remove_inhomogeneous_background_3d(
    image_3d,
    box_size=(64, 64),
    filter_size=(3, 3),
    parallel_execution=True,
    background=False,
):
    """
    Wrapper to remove inhomogeneous background in a 3D image by recursively calling _remove_inhomogeneous_background_2d():
        - addresses output
        - iterates over planes and calls _remove_inhomogeneous_background_2d in each plane
        - reassembles results into a 3D image

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    box_size : tuple of ints, optional
        size of box_size used for block decomposition. The default is (32, 32).
    filter_size : tuple of ints, optional
        Size of gaussian filter used for smoothing results. The default is (3, 3).

    Returns
    -------
    output : numpy array
        processed 3D image.

    """
    if parallel_execution:
        client = try_get_client()
    else:
        client = None

    number_planes = image_3d.shape[0]
    output = np.zeros(image_3d.shape)

    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    if client is not None:
        print_log(
            "> Removing inhomogeneous background from {} planes using {} workers...".format(
                number_planes, len(client.scheduler_info()["workers"])
            )
        )
        image_list = [image_3d[z, :, :] for z in range(number_planes)]
        # image_list_scattered = client.scatter(image_list)

        futures = [
            client.submit(
                Background2D,
                img,
                box_size,
                filter_size=filter_size,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
            )
            for img in image_list
        ]

        results = client.gather(futures)
        print_log(" > Retrieving {} results from cluster".format(len(results)))

        for z, img, bkg in zip(range(number_planes), image_list, results):
            output[z, :, :] = img - bkg.background
        del results, futures, image_list
        # del image_list_scattered

    else:
        print_log(
            "> Removing inhomogeneous background from {} planes using 1 worker...".format(
                number_planes
            )
        )
        z_range = trange(number_planes)
        for z in z_range:
            image_2d = image_3d[z, :, :]
            bkg = Background2D(
                image_2d,
                box_size,
                filter_size=filter_size,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
            )
            output[z, :, :] = image_2d - bkg.background

    if background:
        return output, bkg.background
    else:
        return output


# =============================================================================
# IMAGE ALIGNMENT
# =============================================================================


def apply_xy_shift_3d_images(image, shift, parallel_execution=True):
    """
    Applies XY shift to a 3D stack

    Parameters
    ----------
    images : 3D numpy array
        image to process.

    Returns
    -------
    shifted 3D image.

    """
    if parallel_execution:
        client = try_get_client()
    else:
        client = None

    number_planes = image.shape[0]

    if client is None:
        print_log("> Shifting {} planes with 1 thread...".format(number_planes))
        shift_3d = np.zeros((3))
        shift_3d[0], shift_3d[1], shift_3d[2] = 0, shift[0], shift[1]
        output = shift_image(image, shift_3d)
    else:
        print_log(
            "> Shifting {} planes using {} workers...".format(
                number_planes, len(client.scheduler_info()["workers"])
            )
        )

        image_list_scattered = scatter_3d_image(image)

        futures = [
            client.submit(shift_image, img, shift) for img in image_list_scattered
        ]

        output = reassemble_3d_image(client, futures, image.shape)

        del futures
        del image_list_scattered

    print_log("$ Done shifting 3D image.")

    return output


def image_block_alignment_3d(images, block_size_xy=256, upsample_factor=100):

    # sanity checks
    if len(images) < 2:
        sys.exit("# Error, number of images must be 2, not {}".format(len(images)))

    # - break in blocks
    num_planes = images[0].shape[0]
    block_size = (num_planes, block_size_xy, block_size_xy)

    print_log("$ Breaking images into blocks")
    blocks = [view_as_blocks(x, block_shape=block_size).squeeze() for x in images]

    block_ref = blocks[0]
    block_target = blocks[1]

    # - loop thru blocks and calculates block shift in xyz:
    shift_matrices = [np.zeros(block_ref.shape[0:2]) for x in range(3)]

    # print_log("$ Aligning {} blocks".format(len(block_ref.shape[0])))
    for i in trange(block_ref.shape[0]):
        for j in range(block_ref.shape[1]):
            # - cross correlate in 3D to find 3D shift
            shifts_xyz, _, _ = phase_cross_correlation(
                block_ref[i, j], block_target[i, j], upsample_factor=upsample_factor
            )
            for matrix, _shift in zip(shift_matrices, shifts_xyz):
                matrix[i, j] = _shift

    return shift_matrices, block_ref, block_target


def combine_blocks_image_by_reprojection(
    block_ref, block_target, shift_matrices=None, axis1=0
):
    """
    This routine will overlap block_ref and block_target images block by block.
    block_ref will be used as a template.
    - block_target will be first translated in ZXY using the corresponding values in shift_matrices
    to realign each block
    - then an rgb image will be created with block_ref in the red channel, and the reinterpolated
    block_target block in the green channel.
    - the Blue channel is used for the grid to improve visualization of blocks.


    Parameters
    ----------
    block_ref : npy array
        return of view_as_blocks()
    block_target : npy array
        return of view_as_blocks()
    shift_matrices : list of npy arrays
        index 0 contains Z, index 1 X and index 2 Y
    axis1 : int
        axis used for the reprojection: The default is 0.
        - 0 means an XY projection
        - 1 an ZX projection
        - 2 an ZY projection

    Returns
    -------
    output : NPY array of size im_size x im_size x 3
        rgb image.
    ssim_as_blocks = NPY array of size number_blocks x number_blocks
        Structural similarity index between ref and target blocks
    """
    number_blocks = block_ref.shape[0]
    block_sizes = list(block_ref.shape[2:])
    block_sizes.pop(axis1)
    img_sizes = [x * number_blocks for x in block_sizes]

    # gets ranges for slicing
    slice_coordinates = []
    for block_size in block_sizes:
        slice_coordinates.append(
            [range(x * block_size, (x + 1) * block_size) for x in range(number_blocks)]
        )

    # creates output images
    output = np.zeros((img_sizes[0], img_sizes[1], 3))
    ssim_as_blocks = np.zeros((number_blocks, number_blocks))
    mse_as_blocks = np.zeros((number_blocks, number_blocks))
    nrmse_as_blocks = np.zeros((number_blocks, number_blocks))

    # blank image for blue channel to show borders between blocks
    blue = np.zeros(block_sizes)
    blue[0, :], blue[:, 0], blue[:, -1], blue[-1, :] = [0.5] * 4

    # reassembles image
    # takes one plane block
    for i, i_slice in enumerate(tqdm(slice_coordinates[0])):
        for j, j_slice in enumerate(slice_coordinates[1]):
            imgs = []
            imgs.append(block_ref[i, j])  # appends reference image to image list

            if shift_matrices is not None:
                shift_3d = np.array(
                    [x[i, j] for x in shift_matrices]
                )  # gets 3D shift from block decomposition
                imgs.append(
                    shift_image(block_target[i, j], shift_3d)
                )  # realigns and appends to image list
            else:
                imgs.append(
                    block_target[i, j]
                )  # appends original target with no re-alignment

            imgs = [np.sum(x, axis=axis1) for x in imgs]  # projects along axis1
            imgs = [
                exposure.rescale_intensity(x, out_range=(0, 1)) for x in imgs
            ]  # rescales intensity values
            imgs = [
                image_adjust(x, lower_threshold=0.5, higher_threshold=0.9999)[0]
                for x in imgs
            ]  # adjusts pixel intensities

            nrmse_as_blocks[i, j] = normalized_root_mse(
                imgs[0], imgs[1], normalization="euclidean"
            )
            mse_as_blocks[i, j] = mean_squared_error(imgs[0], imgs[1])
            ssim_as_blocks[i, j] = ssim(
                imgs[0], imgs[1], data_range=imgs[1].max() - imgs[1].min()
            )

            imgs.append(blue)  # appends last channel with grid

            rgb = np.dstack(imgs)  # makes block rgb image

            output[
                i_slice[0] : i_slice[-1] + 1, j_slice[0] : j_slice[-1] + 1, :
            ] = rgb  # inserts block into final rgb stack

    return output, ssim_as_blocks, mse_as_blocks, nrmse_as_blocks


def align_2_images_cross_correlation(
    image1_uncorrected,
    image2_uncorrected,
    lower_threshold=0.999,
    higher_threshold=0.9999999,
    upsample_factor=100,
):
    """
    Aligns 2 images by contrast adjust and cross correlation
    Parameters
    ----------
    img_reference : TYPE
        DESCRIPTION.
    img_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    shift : TYPE
        DESCRIPTION.
    error : TYPE
        DESCRIPTION.
    diffphase : TYPE
        DESCRIPTION.
    lower_threshold : TYPE
        DESCRIPTION.
    i_histogram : TYPE
        DESCRIPTION.
    image2_corrected : TYPE
        DESCRIPTION.
    image1_adjusted : TYPE
        DESCRIPTION.
    image2_adjusted : TYPE
        DESCRIPTION.
    image2_corrected_raw : TYPE
        DESCRIPTION.

    """

    (
        image1_adjusted,
        hist1_before,
        hist1_after,
        lower_cutoff1,
        _,
    ) = image_adjust(
        image1_uncorrected,
        lower_threshold=lower_threshold,
        higher_threshold=higher_threshold,
    )
    (
        image2_adjusted,
        hist2_before,
        hist2_after,
        lower_cutoff2,
        _,
    ) = image_adjust(
        image2_uncorrected,
        lower_threshold=lower_threshold,
        higher_threshold=higher_threshold,
    )

    # zips histograms
    lower_threshold = {"Im1": lower_cutoff1, "Im2": lower_cutoff2}
    i_histogram = {
        "Im1": (hist1_before, hist1_after),
        "Im2": (hist2_before, hist2_after),
    }

    # calculates shift
    shift, error, diffphase = phase_cross_correlation(
        image1_adjusted, image2_adjusted, upsample_factor=upsample_factor
    )

    # corrects image
    # The shift corresponds to the pixel offset relative to the reference image
    image2_corrected = shift_image(image2_adjusted, shift)
    image2_corrected = exposure.rescale_intensity(image2_corrected, out_range=(0, 1))

    results = (
        shift,
        error,
        diffphase,
        lower_threshold,
        i_histogram,
        image2_corrected,
        image1_adjusted,
        image2_adjusted,
    )

    return results


def align_cv2(im1, im2, warp_mode):

    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY:
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else:
        warp_matrix = np.eye(2, 3, dtype=np.float32)

    # Specify the number of iterations.
    number_of_iterations = 1000  # 5000

    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-10

    # Define termination criteria
    criteria = (
        cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT,
        number_of_iterations,
        termination_eps,
    )

    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1
        )
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria
        )
    except cv2.error:
        cc = 0
        # print_log('Warning: find transform failed. Set warp as identity')

    return cc, warp_matrix


def apply_correction(im2, warp_matrix):

    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(
        im2, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP
    )

    return im2_aligned


def align_images_by_blocks(
    img_1,
    img_2,
    block_size,
    upsample_factor=100,
    min_number_pollsters=4,
    tolerance=0.1,
    use_cv2=False,
    shift_error_tolerance=5,
):

    block_1 = view_as_blocks(img_1, block_size)
    block_2 = view_as_blocks(img_2, block_size)

    if use_cv2:
        warp_matrix = np.eye(2, 3, dtype=np.float32)
        warp_mode = cv2.MOTION_TRANSLATION

    shift_image_norm = np.zeros((block_1.shape[0], block_1.shape[1]))
    shifted_image = np.zeros((block_1.shape[0], block_1.shape[1], 2))
    rms_image = np.zeros((block_1.shape[0], block_1.shape[1]))

    for i in trange(block_1.shape[0]):
        for j in range(block_1.shape[1]):
            if not use_cv2:
                # using Scimage registration functions
                shift, _, _ = phase_cross_correlation(
                    block_1[i, j], block_2[i, j], upsample_factor=upsample_factor
                )
                shift_image_norm[i, j] = LA.norm(shift)
                shifted_image[i, j, 0], shifted_image[i, j, 1] = shift[0], shift[1]
                img_2_aligned = shift_image(img_2, shift)
            else:
                # uses CV2 cause it is 20 times faster than Scimage
                _, warp_matrix = align_cv2(block_1[i, j], block_2[i, j], warp_mode)
                shift_image_norm[i, j] = LA.norm(warp_matrix[:, 2])
                shifted_image[i, j, 0], shifted_image[i, j, 1] = (
                    warp_matrix[:, 2][0],
                    warp_matrix[:, 2][1],
                )
                img_2_aligned = apply_correction(img_2, warp_matrix)

            rms_image[i, j] = np.sum(np.sum(np.abs(img_1 - img_2_aligned), axis=1))

    # [calculates optimal shifts by polling blocks showing the best RMS]

    # threshold = filters.threshold_otsu(rms_image)
    threshold = (1 + tolerance) * np.min(rms_image)
    mask = rms_image < threshold

    contours = measure.find_contours(rms_image, threshold)

    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour = np.array([0, 0])

    # [Averages shifts and errors from regions within the tolerated blocks]
    mean_shifts = [np.mean(shifted_image[mask, 0]), np.mean(shifted_image[mask, 1])]
    std_shifts = [np.std(shifted_image[mask, 0]), np.std(shifted_image[mask, 1])]
    mean_shift_norm = np.mean(shift_image_norm[mask])
    mean_error = np.mean(rms_image[mask])
    relative_shifts = np.abs(shift_image_norm - mean_shift_norm)

    # [calculates global shift, if it is better than the polled shift, or
    # if we do not have enough pollsters to fall back to then it does a global cross correlation!]
    mean_shifts_global, _, _ = phase_cross_correlation(img_1, img_2, upsample_factor=100)
    img_2_aligned_global = shift_image(img_2, shift)
    mean_error_global = np.sum(np.sum(np.abs(img_1 - img_2_aligned_global), axis=1))

    print_log(
        "Block alignment error: {}, global alignment error: {}".format(
            mean_error, mean_error_global
        )
    )

    if (
        np.sum(mask) < min_number_pollsters
        or mean_error_global < mean_error
        or np.max(std_shifts) > shift_error_tolerance
    ):
        mean_shifts = mean_shifts_global
        mean_error = mean_error_global
        print_log("Falling back to global registration")

    print_log(
        "*** Global XY shifts: {:.2f} px | {:.2f} px".format(
            mean_shifts_global[0], mean_shifts_global[1]
        )
    )
    print_log(
        "*** Mean polled XY shifts: {:.2f}({:.2f}) px | {:.2f}({:.2f}) px".format(
            mean_shifts[0], std_shifts[0], mean_shifts[1], std_shifts[1]
        )
    )

    return np.array(mean_shifts), mean_error, relative_shifts, rms_image, contour


# =============================================================================
# FOCAL PLANE INTERPOLATION
# =============================================================================


def focal_plane(data, threshold_fwhm=20, verbose=False):
    """
    This function will find the focal plane of a 3D image
    - calculates the laplacian variance of the image for each z plane
    - fits 1D gaussian profile on the laplacian variance
    - to get the maximum (focal plane) and the full width at half maximum
    - it returns nan if the fwhm > threshold_fwhm (means fit did not converge)
    - threshold_fwhm should represend the width of the laplacian variance curve
    - which is often 5-10 planes depending on sample.

    Parameters
    ----------
    data : numpy array
        input 3D image ZYX.
    threshold_fwhm : float, optional
        threshold fwhm used to remove outliers. The default is 20.

    Returns
    -------
    focal_plane : float
        focal plane: max of fitted z-profile.
    fwhm : float
        full width hald maximum of fitted z-profile.

    """
    # finds focal plane
    raw_images = [data[i, :, :] for i in range(data.shape[0])]
    laplacian_variance = [cv2.Laplacian(img, cv2.CV_64F).var() for img in raw_images]
    laplacian_variance = laplacian_variance / max(laplacian_variance)
    x_coord = range(len(laplacian_variance))
    fit_result, _ = fit_1d_gaussian_scipy(
        x_coord,
        laplacian_variance,
        title="laplacian variance z-profile",
        verbose=verbose,
    )

    if len(fit_result) > 0:
        focal_plane = fit_result["gauss1d.pos"]
        fwhm = fit_result["gauss1d.fwhm"]

        if fwhm > threshold_fwhm:
            fwhm, focal_plane = np.nan, np.nan
    else:
        fwhm, focal_plane = np.nan, np.nan

    if verbose:
        print_log("Focal plane found: {} with fwhm: {}".format(focal_plane, fwhm))

    return focal_plane, fwhm


def calculate_focus_per_block(data, block_size_xy=128):
    """
    Calculates the most likely focal plane of an image by breaking into blocks and calculating
    the focal plane in each block

    - breaks image into blocks
    - returns the focal plane + fwhm for each block

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    block_size_xy : TYPE, optional
        DESCRIPTION. The default is 512.

    Returns
    -------
    focal_plane_matrix: np array
        matrix containing the maximum of the laplacian variance per block
    fwhm: np array
        matrix containing the fwhm of the laplacian variance per block
    block: np array
        3D block reconstruction of matrix

    """
    n_planes = data.shape[0]

    block_size = (n_planes, block_size_xy, block_size_xy)

    block = view_as_blocks(data, block_size).squeeze()
    focal_plane_matrix = np.zeros(block.shape[0:2])
    fwhm = np.zeros(block.shape[0:2])

    fwhm = {}

    for i in trange(block.shape[0]):
        # fwhm[str(i)] = {}
        for j in range(block.shape[1]):
            focal_plane_matrix[i, j], fwhm[i, j] = focal_plane(block[i, j])

    return focal_plane_matrix, fwhm, block


# Program to find most frequent
# element in a list
def most_frequent(elt_list):
    return max(set(elt_list), key=elt_list.count)


def reassemble_images(focal_plane_matrix, block, window=0):
    """
    Makes 2D image from 3D stack by reassembling sub-blocks
    For each sub-block we know the optimal focal plane, which is
    selected for the assembly of the while image

    Parameters
    ----------
    focal_plane_matrix : numpy 2D array
        matrix containing the focal plane selected for each block.
    block : numpy matrix
        original 3D image sorted by blocks.

    Returns
    -------
    output : numpy 2D array
        output 2D projection

    """
    # gets image size from block image
    number_blocks = block.shape[0]
    block_size_xy = block.shape[3]
    im_size = number_blocks * block_size_xy

    # gets ranges for slicing
    slice_coordinates = [
        range(x * block_size_xy, (x + 1) * block_size_xy) for x in range(number_blocks)
    ]

    # creates output image
    output = np.zeros((im_size, im_size))

    # gets more common plane
    focal_planes = []
    for i, i_slice in enumerate(slice_coordinates):
        for j, j_slice in enumerate(slice_coordinates):
            focal_planes.append(int(focal_plane_matrix[i, j]))
    most_common_focal_plane = most_frequent(focal_planes)

    # reassembles image
    if window == 0:
        # takes one plane block
        for i, i_slice in enumerate(slice_coordinates):
            for j, j_slice in enumerate(slice_coordinates):
                focus = int(focal_plane_matrix[i, j])
                if np.abs(focus - most_common_focal_plane) > 1:
                    focus = int(most_common_focal_plane)
                output[
                    i_slice[0] : i_slice[-1] + 1, j_slice[0] : j_slice[-1] + 1
                ] = block[i, j][focus, :, :]
    else:
        # takes neighboring planes by projecting
        for i, i_slice in enumerate(slice_coordinates):
            for j, j_slice in enumerate(slice_coordinates):
                focus = int(focal_plane_matrix[i, j])
                if np.abs(focus - most_common_focal_plane) > 1:
                    focus = int(most_common_focal_plane)
                zmin = np.max((0, focus - round(window / 2)))
                zmax = np.min((block[i, j].shape[0], focus + round(window / 2)))
                z_range = (focus, range(zmin, zmax))
                output[
                    i_slice[0] : i_slice[-1] + 1, j_slice[0] : j_slice[-1] + 1
                ] = project_image_2d(block[i, j][:, :, :], z_range, "MIP")

    return output


def reinterpolate_focal_plane(data, param_dict):

    if "blockSize" in param_dict["zProject"]:
        block_size_xy = param_dict["zProject"]["blockSize"]
    else:
        block_size_xy = 128

    if "zwindows" in param_dict["zProject"]:
        window = param_dict["zProject"]["zwindows"]
    else:
        window = 0

    focal_plane_matrix, z_range, block = _reinterpolate_focal_plane(
        data, block_size_xy=block_size_xy, window=window
    )

    # reassembles image
    output = reassemble_images(focal_plane_matrix, block, window=window)

    return output, focal_plane_matrix, z_range


def _reinterpolate_focal_plane(data, block_size_xy=256, window=10):
    """
    Reinterpolates the focal plane of a 3D image by breking it into blocks
    - Calculates the focal_plane and fwhm matrices by block
    - removes outliers
    - calculates the focal plane for each block using sigmaClip statistics
    - returns a tuple with focal plane and the range to use

    Parameters
    ----------
    data : numpy array
        input 3D image.
    block_size_xy : int
        size of blocks in XY, typically 256.
    window : int, optional
        number of planes before and after the focal plane to construct the z_range. The default is 0.

    Returns
    -------
    focal_plane_matrix : numpy array
        focal plane matrix.
    z_range : tuple
        focus_plane, z_range.
    block : numpy array
        block representation of 3D image.

    """
    print_log("> Reducing 3D image size by slicing around focal plane")
    # breaks into subplanes, iterates over them and calculates the focal_plane in each subplane.

    focal_plane_matrix, _, block = calculate_focus_per_block(
        data, block_size_xy=block_size_xy
    )

    focal_planes_to_process = focal_plane_matrix[~np.isnan(focal_plane_matrix)]

    focal_plane, _, _ = sigmaclip(focal_planes_to_process, high=3, low=3)
    focus_plane = np.mean(focal_plane)
    if np.isnan(focus_plane):
        print_log("# focus_plane detection failed. Using full stack.")
        focus_plane = data.shape[0] // 2
        z_range = focus_plane, range(0, data.shape[0])
    else:
        focus_plane = np.mean(focal_plane).astype("int64")
        zmin = np.max([focus_plane - window, 0])
        zmax = np.min([focus_plane + window, data.shape[0]])
        z_range = focus_plane, range(zmin, zmax)

    return focal_plane_matrix, z_range, block


def _remove_z_planes(image_3d, z_range):
    """
    Removes planes in input image.
    For instance, if you provide a z_range = range(0,image_3d.shape[0],2)

    then the routine will remove any other plane. Number of planes skipped
    can be programmed by tuning z_range.

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))
    for i, index in enumerate(z_range):
        output[i, :, :] = image_3d[index, :, :]

    return output


def _average_z_planes(image_3d, z_range):
    """
    Removes z-planes by calculating the average between successive planes

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """

    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))
    for i, index in enumerate(z_range):
        average = (
            image_3d[index, :, :].astype(np.float)
            + image_3d[index + 1, :, :].astype(np.float)
        ) / 2
        output[i, :, :] = average.astype(np.uint16)

    return output


def _interpolate_z_planes(image_3d, z_range):
    """
    Removes z planes by reinterpolation

    TODO

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """

    output = np.zeros((len(z_range), image_3d.shape[1], image_3d.shape[2]))

    # need to code using interpn
    output = image_3d

    return output


def reinterpolate_z(image_3d, z_range, mode="average"):
    """
    wrapper function for any kind of z-interpolation
    to reduce the number of planes in an image

    Parameters
    ----------
    image_3d : numpy array
        input 3D image.
    z_range : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    if "interpolate" in mode:
        output = _interpolate_z_planes(image_3d, z_range)
    elif "remove" in mode:
        output = _remove_z_planes(image_3d, z_range)
    elif "average" in mode:
        output = _average_z_planes(image_3d, z_range)

    print_log(
        "$ Reduced Z-planes from {} to {}".format(image_3d.shape[0], output.shape[0])
    )

    return output


# =============================================================================
# SEGMENTATION FUNCTIONS
# =============================================================================


def _segment_2d_image_by_thresholding(
    image_2d,
    threshold_over_std=10,
    area_min=3,
    area_max=1000,
    nlevels=64,
    contrast=0.001,
    kernel=None,
):

    # makes threshold matrix
    threshold = np.zeros(image_2d.shape)
    threshold[:] = threshold_over_std * image_2d.max() / 100

    # segments objects
    segm = detect_sources(image_2d, threshold, npixels=area_min, filter_kernel=kernel,)

    if segm.nlabels > 0:
        # removes masks too close to border
        segm.remove_border_labels(border_width=10)  # parameter to add to infoList

        if segm.nlabels > 0:
            segm_deblend = deblend_sources(
                image_2d,
                segm,
                npixels=area_min,  # watch out, this is per plane!
                filter_kernel=kernel,
                nlevels=nlevels,
                contrast=contrast,
                relabel=True,
                mode="exponential",
            )
            if segm_deblend.nlabels > 0:
                # removes Masks too big or too small
                for label in segm_deblend.labels:
                    # take regions with large enough areas
                    area = segm_deblend.get_area(label)
                    if area < area_min or area > area_max:
                        segm_deblend.remove_label(label=label)

                # relabel so masks numbers are consecutive
                # segm_deblend.relabel_consecutive()

            # image_2d_segmented = segm.data % changed during recoding function
            image_2d_segmented = segm_deblend.data

            image_2d_segmented[image_2d_segmented > 0] = 1
            return image_2d_segmented

        else:
            # returns empty image as no objects were detected
            return segm.data

    else:
        # returns empty image as no objects were detected
        return segm.data


def _segment_3d_volumes_stardist(
    image_3d,
    deblend_3d=False,
    axis_norm=(0, 1, 2),
    model_dir="/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/",
    model_name="stardist_18032021_single_loci",
):

    number_planes = image_3d.shape[0]

    print_log("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print_log("> Segmenting {} planes using 1 worker...".format(number_planes))
    print_log("> Loading model {} from {}...".format(model_name, model_dir))
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    model = StarDist3D(None, name=model_name, basedir=model_dir)
    limit_gpu_memory(None, allow_growth=True)

    im = normalize(image_3d, 1, 99.8, axis=axis_norm)
    l_x = im.shape[1]

    if l_x < 1000:
        labels, _ = model.predict_instances(im)

    else:
        resizer = PadAndCropResizer()
        axes = "ZYX"

        im = resizer.before(im, axes, model._axes_div_by(axes))
        labels, _ = model.predict_instances(im, n_tiles=(1, 8, 8))
        labels = resizer.after(labels, axes)

    mask = np.array(labels > 0, dtype=int)

    # Now we want to separate objects in 3D using watersheding
    if deblend_3d:
        labeled_image = _deblend_3d_segmentation(mask)
    else:
        labeled_image = labels

    return mask, labeled_image


def _segment_3d_volumes_by_thresholding(
    image_3d,
    threshold_over_std=10,
    sigma=3,
    box_size=(32, 32),
    filter_size=(3, 3),
    area_min=3,
    area_max=1000,
    nlevels=64,
    contrast=0.001,
    deblend_3d=False,
    parallel_execution=True,
):
    if parallel_execution:
        client = try_get_client()
    else:
        client = None

    number_planes = image_3d.shape[0]

    kernel = Gaussian2DKernel(sigma, x_size=sigma, y_size=sigma)
    kernel.normalize()

    print_log("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

    parallel = True
    if client is None:
        parallel = False
    else:
        if len(client.scheduler_info()["workers"]) < 1:
            parallel = False
            print_log("# Failed getting workers. Report of scheduler:")
            for key in client.scheduler_info().keys():
                print_log("{}:{}".format(key, client.scheduler_info()[key]))

    if not parallel:
        print_log("> Segmenting {} planes using 1 worker...".format(number_planes))

        output = np.zeros(image_3d.shape)

        for z in trange(number_planes):
            image_2d = image_3d[z, :, :]
            image_2d_segmented = _segment_2d_image_by_thresholding(
                image_2d,
                threshold_over_std=threshold_over_std,
                area_min=area_min,
                area_max=area_max,
                nlevels=nlevels,
                contrast=contrast,
                kernel=kernel,
            )
            output[z, :, :] = image_2d_segmented

    else:

        print_log(
            "> Segmenting {} planes using {} workers...".format(
                number_planes, len(client.scheduler_info()["workers"])
            )
        )

        image_list_scattered = scatter_3d_image(image_3d)

        futures = [
            client.submit(
                _segment_2d_image_by_thresholding,
                img,
                threshold_over_std=threshold_over_std,
                sigma=sigma,
                box_size=box_size,
                filter_size=filter_size,
                area_min=area_min,
                area_max=area_max,
                nlevels=nlevels,
                contrast=contrast,
                deblend_3d=deblend_3d,
                kernel=kernel,
            )
            for img in image_list_scattered
        ]

        output = reassemble_3d_image(client, futures, image_3d.shape)

        del futures, image_list_scattered

    labels = measure.label(output)

    # Now we want to separate objects in 3D using watersheding
    if deblend_3d:
        labels = _deblend_3d_segmentation(output)

    return output, labels


def _deblend_3d_segmentation(binary):
    """
    Deblends objects in 3D using watersheding

    Parameters
    ----------
    binary : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    binary = binary > 0
    print_log(" > Constructing distance matrix from 3D binary mask...")

    distance = apply_parallel(ndi.distance_transform_edt, binary)

    print_log(" > Deblending sources in 3D by watersheding...")
    coords = peak_local_max(distance, footprint=np.ones((10, 10, 25)), labels=binary)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)

    labels = watershed(-distance, markers, mask=binary)
    return labels


def _segment_3d_masks(
    image_3d,
    axis_norm=(0, 1, 2),
    pmin=1,
    pmax=99.8,
    model_dir="/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
    model_name="stardist_20210625_deconvolved",
):

    """
    Parameters
    ----------
    image_3d : numpy ndarray (N-dimensional array)
        3D raw image to be segmented

    model_dir : List of strings, optional
        paths of all models directory, the default is ["/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks"]

    model_name : List of strings, optional
        names of all models, the default is ['stardist_20210625_deconvolved']

    """

    np.random.seed(6)

    number_planes = image_3d.shape[0]

    print_log("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print_log("> Segmenting {} planes using 1 worker...".format(number_planes))
    print_log("> Loading model {} from {}...".format(model_name, model_dir))
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"  # why do we need this?

    # Load the model
    # --------------

    model = StarDist3D(None, name=model_name, basedir=model_dir)
    limit_gpu_memory(None, allow_growth=True)

    im = normalize(image_3d, pmin=pmin, pmax=pmax, axis=axis_norm)
    l_x = im.shape[1]

    if l_x < 351:  # what is this value? should it be a k-arg?
        labels, _ = model.predict_instances(im)

    else:
        resizer = PadAndCropResizer()
        axes = "ZYX"

        im = resizer.before(im, axes, model._axes_div_by(axes))
        labels, _ = model.predict_instances(im, n_tiles=(1, 8, 8))
        labels = resizer.after(labels, axes)

    mask = np.array(labels > 0, dtype=int)
    mask[mask > 0] = 1

    return mask, labels


def plot_raw_images_and_labels(image, label):

    """
    Parameters
    ----------
    image : List of numpy ndarray (N-dimensional array)
        3D raw image of format .tif

    label : List of numpy ndarray (N-dimensional array)
        3D labeled image of format .tif
    """

    cmap = random_label_cmap()

    moy = np.mean(image, axis=0)
    lbl_moy = np.max(label, axis=0)

    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches((50, 50))
    ax = axes.ravel()
    titles = ["raw image", "projected labeled image"]

    ax[0].imshow(moy, cmap="Greys_r", origin="lower")

    ax[1].imshow(lbl_moy, cmap=cmap, origin="lower")

    for axis, title in zip(ax, titles):
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(title)

    return fig


########################################################
# SAVING ROUTINES
########################################################


def save_2_images_rgb(img_1, img_2, output_filename):
    """
    Overlays two images as R and B and saves them to output file
    """

    sz = img_1.shape
    img_1, img_2 = img_1 / img_1.max(), img_2 / img_2.max()

    img_1, _, _, _, _ = image_adjust(img_1, lower_threshold=0.5, higher_threshold=0.9999)
    img_2, _, _, _, _ = image_adjust(img_2, lower_threshold=0.5, higher_threshold=0.9999)

    fig, ax1 = plt.subplots()
    fig.set_size_inches((30, 30))

    null_image = np.zeros(sz)

    rgb = np.dstack([img_1, img_2, null_image])
    ax1.imshow(rgb)
    ax1.axis("off")

    fig.savefig(output_filename)

    plt.close(fig)


def save_image_2d_cmd(image, file_name):
    if image.shape > (1, 1):
        np.save(file_name, image)
        # log.report("Saving 2d projection to disk:{}\n".format(os.path.basename(file_name)),'info')
        print_log("$ Image saved to disk: {}".format(file_name + ".npy"), "info")
    else:
        print_log("# Warning, image is empty", "Warning")


def save_image_as_blocks(img, full_filename, block_size_xy=256, label="raw_image"):
    num_planes = img.shape[0]
    block_size = (num_planes, block_size_xy, block_size_xy)
    blocks = view_as_blocks(img, block_shape=block_size).squeeze()
    print_log(
        "\nDecomposing image into {} blocks".format(blocks.shape[0] * blocks.shape[1])
    )

    folder = full_filename.split(".")[0]
    file_name = os.path.basename(full_filename).split(".")[0]

    if not os.path.exists(folder):
        os.mkdir(folder)
        print_log("Folder created: {}".format(folder))

    for i in trange(blocks.shape[0]):
        for j in range(blocks.shape[1]):
            outfile = (
                folder
                + os.sep
                + file_name
                + "_"
                + label
                + "_block_"
                + str(i)
                + "_"
                + str(j)
                + ".tif"
            )
            imsave(outfile, blocks[i, j])

    # cmap = "seismic"


def image_show_with_values_single(
    ax, matrix, cbarlabel, fontsize, cbar_kw, valfmt="{x:.0f}", cmap="YlGn"
):
    row = ["".format(x) for x in range(matrix.shape[0])]
    im, _ = heatmap(
        matrix,
        row,
        row,
        ax=ax,
        cmap=cmap,
        cbarlabel=cbarlabel,
        fontsize=fontsize,
        cbar_kw=cbar_kw,
    )
    _ = annotate_heatmap(
        im, valfmt=valfmt, size=fontsize, threshold=None, textcolors=("black", "white")
    )  # , fontsize=fontsize


def image_show_with_values(
    matrices,
    output_name="tmp.png",
    cbarlabels=["focalPlane"],
    fontsize=6,
    verbose=False,
    title="",
):
    """
    Plots a list of matrices with their values in each pixel.

    Parameters
    ----------
    matrices : list
        matrices to plot. Should be 2D numpy arrays
    output_name : TYPE, optional
        DESCRIPTION. The default is "tmp.png".
    cbarlabels : list, optional
        titles of subplots. The default is ["focalPlane"].
    fontsize : float, optional
        fontsize. The default is 6.
    verbose : Boolean, optional
        self explanatory. The default is False.
    title : str, optional
        figure title. The default is "".

    Returns
    -------
    None.

    """
    number_images = len(matrices)
    fig, axes = plt.subplots(1, number_images)
    fig.set_size_inches((number_images * 2, 5))
    ax = axes.ravel()
    fig.suptitle(title)
    cbar_kw = {}
    cbar_kw["fraction"] = 0.046
    cbar_kw["pad"] = 0.04

    if len(cbarlabels) != number_images:
        cbarlabels = cbarlabels[0] * number_images

    for matrix, axis, cbarlabel in zip(matrices, ax, cbarlabels):
        image_show_with_values_single(axis, matrix, cbarlabel, fontsize, cbar_kw)

    fig.tight_layout()
    plt.savefig(output_name)

    if not verbose:
        plt.close(fig)


def display_3d(
    image_3d=None,
    labels=None,
    localizations_list=None,
    z=40,
    range_xy=1000,
    norm=True,
    cmap="Greys",
):

    if image_3d is not None:
        images = []
        images.append(image_3d[z, :, :])
        images.append(image_3d[:, range_xy, :])
        images.append(image_3d[:, :, range_xy])
    else:
        images = [1, 1, 1]

    if labels is not None:
        segmented = []
        segmented.append(labels[z, :, :])
        segmented.append(labels[:, range_xy, :])
        segmented.append(labels[:, :, range_xy])
    else:
        segmented = [1, 1, 1]

    if localizations_list is not None:
        localized_list = []

        for localizations in localizations_list:
            localized = []
            localized.append(localizations[:, [2, 1]])
            localized.append(localizations[:, [2, 0]])
            localized.append(localizations[:, [1, 0]])
            localized_list.append(localized)

    else:
        localized_list = [1, 1, 1]

    percent = 99.5
    symbols = ["+", "o", "*", "^"]
    colors = ["r", "b", "g", "y"]

    fig, axes = plt.subplots(1, len(images))
    fig.set_size_inches(len(images) * 50, 50)
    ax = axes.ravel()

    for image, segm, axis, i_plane in zip(images, segmented, ax, range(len(ax))):
        if image_3d is not None:
            if norm:
                norm = simple_norm(image, "sqrt", percent=percent)
                axis.imshow(image, cmap=cmap, origin="lower", norm=norm)
            else:
                axis.imshow(image, cmap=cmap, origin="lower")
        if labels is not None:
            axis.imshow(color.label2rgb(segm, bg_label=0), alpha=0.3)
        if localizations_list is not None:
            for i_loc_list, symbol, color in zip(
                range(len(localized_list)), symbols, colors
            ):
                locs = localized_list[i_loc_list][i_plane]
                axis.plot(locs[:, 0], locs[:, 1], symbol, color=color, alpha=0.7)

    return fig


def display_3d_assembled(
    images, localizations=None, plotting_range=None, normalize=True, masks=None
):
    wspace = 25

    rows_xy, cols_xy = images[0].shape[1], images[0].shape[0]
    rows_yz, cols_zx = images[2].shape[0], images[1].shape[0]
    rows = rows_xy + rows_yz + wspace
    cols = cols_xy + cols_zx + wspace

    fig = plt.figure(figsize=(50, 50), tight_layout=True)
    axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    display_image = np.zeros((rows, cols))
    display_image[0:rows_xy, 0:cols_xy] = images[0]
    display_image[rows_xy + wspace : rows_xy + rows_yz + wspace, 0:cols_xy] = images[2]
    display_image[0:rows_xy, cols_xy + wspace : cols_xy + cols_zx + wspace] = images[
        1
    ].transpose()

    if normalize:
        norm = simple_norm(display_image[:, :], "sqrt", percent=99)
        axis.imshow(display_image[:, :], cmap="Greys", alpha=1, norm=norm)
    else:
        axis.imshow(display_image[:, :], cmap="Greys", alpha=1)

    if masks is not None:
        display_mask = np.zeros((rows, cols))
        display_mask[0:rows_xy, 0:cols_xy] = masks[0]
        display_mask[rows_xy + wspace : rows_xy + rows_yz + wspace, 0:cols_xy] = masks[
            2
        ]
        display_mask[0:rows_xy, cols_xy + wspace : cols_xy + cols_zx + wspace] = masks[
            1
        ].transpose()

        axis.imshow(color.label2rgb(display_mask[:, :], bg_label=0), alpha=0.3)

    colors = ["r", "g", "b", "y"]
    markersizes = [2, 1, 1, 1, 1]
    if localizations is not None:
        for i, loc in enumerate(localizations):
            axis.plot(loc[:, 2], loc[:, 1], "+", color=colors[i], markersize=1)
            if plotting_range is not None:
                selections = [
                    np.abs(loc[:, a] - plotting_range[0]) < plotting_range[1]
                    for a in [2, 1]
                ]
                axis.plot(
                    loc[selections[0], 0] + cols_xy + wspace,
                    loc[selections[0], 1],
                    "+",
                    color=colors[i],
                    markersize=markersizes[i],
                    alpha=0.9,
                )
                axis.plot(
                    loc[selections[1], 2],
                    loc[selections[1], 0] + rows_xy + wspace,
                    "+",
                    color=colors[i],
                    markersize=markersizes[i],
                    alpha=0.9,
                )
            else:
                axis.plot(
                    loc[:, 0] + cols_xy, loc[:, 1], "+", color=colors[i], markersize=1
                )

    axis.axes.yaxis.set_visible(False)
    axis.axes.xaxis.set_visible(False)
    axis.axis("off")

    return fig


def _plot_image_3d(image_3d, localizations=None, masks=None, normalize=True, window=3):
    """
    makes list with XY, XZ and ZY projections and sends for plotting

    Parameters
    ----------
    image_3d : numpy array
        image in 3D.
    localizations : list
        list of localizations to overlay onto 3D image. The default is None.

    Returns
    -------
    figure handle

    """
    img = image_3d
    center = int(img.shape[1] / 2)

    images = []
    images.append(np.sum(img, axis=0))
    images.append(np.sum(img[:, :, center - window : center + window], axis=2))
    images.append(np.sum(img[:, center - window : center + window, :], axis=1))

    if masks is not None:
        labels = []
        labels.append(np.max(masks, axis=0))
        labels.append(np.max(masks[:, :, center - window : center + window], axis=2))
        labels.append(np.max(masks[:, center - window : center + window, :], axis=1))
    else:
        labels = None

    fig1 = display_3d_assembled(
        images,
        localizations=localizations,
        masks=labels,
        plotting_range=[center, window],
        normalize=normalize,
    )

    return fig1


def plot_3d_shift_matrices(shift_matrices, fontsize=8, log=False, valfmt="{x:.1f}"):

    cbar_kw = {}
    cbar_kw["fraction"] = 0.046
    cbar_kw["pad"] = 0.04

    fig, axes = plt.subplots(1, len(shift_matrices))
    fig.set_size_inches((len(shift_matrices) * 10, 10))
    ax = axes.ravel()
    titles = ["z shift matrix", "x shift matrix", "y shift matrix"]

    for axis, title, x in zip(ax, titles, shift_matrices):
        if log:
            x = np.log10(x)
        image_show_with_values_single(
            axis, x, title, fontsize, cbar_kw, valfmt=valfmt, cmap="YlGn"
        )  # YlGnBu
        axis.set_title(title)

    return fig


def plot_4_images(
    allimages,
    titles=["reference", "cycle <i>", "processed reference", "processed cycle <i>"],
):

    fig, axes = plt.subplots(2, 2)
    fig.set_size_inches((10, 10))
    ax = axes.ravel()

    for axis, img, title in zip(ax, allimages, titles):
        # im = np.sum(img, axis=0)
        axis.imshow(img, cmap="Greys")
        axis.set_title(title)
    fig.tight_layout()

    return fig


def plotting_block_alignment_results(
    relative_shifts, rms_image, contour, file_name="BlockALignmentResults.png"
):

    # plotting
    fig, axes = plt.subplots(1, 2)
    ax = axes.ravel()
    fig.set_size_inches((10, 5))

    cbwindow = 3
    p_1 = ax[0].imshow(relative_shifts, cmap="terrain", vmin=0, vmax=cbwindow)
    ax[0].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[0].set_title("abs(global-block) shifts, px")
    fig.colorbar(p_1, ax=ax[0], fraction=0.046, pad=0.04)

    p_2 = ax[1].imshow(
        rms_image, cmap="terrain", vmin=np.min(rms_image), vmax=np.max(rms_image)
    )
    ax[1].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[1].set_title("RMS")
    fig.colorbar(p_2, ax=ax[1], fraction=0.046, pad=0.04)

    for axe in ax:
        axe.axis("off")

    fig.savefig(file_name)

    plt.close(fig)


def heatmap(
    data,
    row_labels,
    col_labels,
    ax=None,
    cbar_kw=None,
    cbarlabel="",
    fontsize=12,
    **kwargs
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)
    # im = ax.imshow(data)
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", size=fontsize)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, size=fontsize)
    ax.set_yticklabels(row_labels, size=fontsize)

    # Let the horizontal axes labeling appear on top.
    # ax.tick_params(top=True, bottom=False,
    #                 labeltop=True, labelbottom=False)

    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for _, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.1f}",
    textcolors=("black", "white"),
    threshold=None,
    **textkw
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    new_kwargs = dict(horizontalalignment="center", verticalalignment="center")
    new_kwargs.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            new_kwargs.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **new_kwargs)
            texts.append(text)

    return texts
