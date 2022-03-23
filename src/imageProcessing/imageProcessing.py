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

import os,sys
from tqdm import trange, tqdm
import numpy as np
from numpy import linalg as LA
from matplotlib import ticker
import matplotlib.pyplot as plt

import scipy.optimize as spo
from scipy.ndimage import shift as shiftImage
from scipy import ndimage as ndi
from scipy.stats import sigmaclip

import cv2
from tifffile import imsave

from skimage import io
from skimage import measure
from skimage import exposure,color
from skimage.util.shape import view_as_blocks
from skimage.registration import phase_cross_correlation
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error, normalized_root_mse
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.util.apply_parallel import apply_parallel

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch, simple_norm
from astropy.stats import SigmaClip
from astropy.convolution import Gaussian2DKernel

from photutils import detect_sources
from photutils import deblend_sources
from photutils import Background2D, MedianBackground

from csbdeep.utils import normalize
from stardist.models import StarDist3D
from csbdeep.utils.tf import limit_gpu_memory
from csbdeep.data import PadAndCropResizer

from fileProcessing.fileManagement import try_get_client, printLog
import warnings

warnings.filterwarnings("ignore")

np.seterr(divide='ignore', invalid='ignore')

# =============================================================================
# CLASSES
# =============================================================================


class Image:
    def __init__(self, param=dict(), log1=[]):
        self.param = param
        self.log = log1
        self.data = []
        self.fileName = ""
        self.data_2D = np.zeros((1, 1))
        self.stageCoordinates = [0.0, 0.0]
        self.imageSize = -1
        self.focusPlane = -1
        self.extension = ""

    # read an image as a numpy array
    def loadImage(self, fileName):
        self.data = io.imread(fileName).squeeze()
        self.fileName = fileName
        self.imageSize = self.data.shape
        self.extension = fileName.split(".")[-1]

    # save 2D projection as numpy array
    def saveImage2D(self, rootFolder, tag="_2d"):
        fileName = self.getImageFileName(rootFolder, tag)
        saveImage2Dcmd(self.data_2D, fileName)

    def getImageFileName(self, rootFolder, tag):
        fileName = rootFolder + os.sep + os.path.basename(self.fileName).split(".")[0] + tag
        return fileName

    # read an image as a numpy array
    def loadImage2D(self, fileName, masterFolder, tag="_2d"):
        self.fileName = fileName
        fileName = self.getImageFileName(masterFolder, tag) + ".npy"

        self.data_2D = np.load(fileName)
        printLog("$ Loading from disk:{}".format(os.path.basename(fileName)))

    # max intensity projection using all z planes
    def maxIntensityProjection(self):
        self.data_2D = np.max(self.data, axis=0)

    # Normalize a 3d image <im> by subtracting local gaussian blur of std <sz>
    def normalizeImage(self, sz=30, ratio=False):

        # if not ratio:
        background = self.data_2D.min()
        _im = self.data_2D - background
        _max = _im.max()
        _im = _im / _max
        return _im

    # returns the imagename
    def getImageFilename(self):
        return self.fileName

    # returns the picture x,y location, if available
    def getImageLocation(self):
        if hasattr(self, "stageCoordinates"):
            return self.stageCoordinates
        else:
            return [0.0, 0.0]

    # returns the film size
    def getImageSize(self):
        return self.imageSize

    # returns the film focus
    def getFocusPlane(self):
        if hasattr(self, "focusPlane"):
            return self.focusPlane
        else:
            return 0.0

    # Outputs image properties to command line
    def printImageProperties(self):
        # printLog("Image Name={}".format(self.fileName))
        printLog("$ Image Size={}".format(self.imageSize))
        # self.log.report("Stage position={}".format(self.stageCoordinates))
        printLog("$ Focal plane={}".format(self.focusPlane))

    # processes sum image in axial direction given range
    # @jit(nopython=True)
    def zProjectionRange(self):

        # find the correct range for the projection
        if self.param.param["zProject"]["zmax"] > self.imageSize[0]:
            printLog("$ Setting z max to the last plane")
            self.param.param["zProject"]["zmax"] = self.imageSize[0]

        if self.param.param["zProject"]["mode"] == "automatic":
            printLog("> Calculating planes...")
            zRange = calculate_zrange(self.data, self.param)

        elif self.param.param["zProject"]["mode"] == "full":
            (zmin, zmax) = (0, self.imageSize[0])
            zRange = (round((zmin + zmax) / 2), range(zmin, zmax))

        if self.param.param["zProject"]["mode"] == "laplacian":
            printLog("Stacking using Laplacian variance...")
            (
                self.data_2D,
                self.focalPlaneMatrix,
                self.zRange,
            ) = reinterpolatesFocalPlane(self.data, self.param.param)
            self.focusPlane = self.zRange[0]

        else:
            # Manual: reads from parameters file
            (zmin, zmax) = (
                self.param.param["zProject"]["zmin"],
                self.param.param["zProject"]["zmax"],
            )
            if zmin >= zmax:
                raise SystemExit("zmin is equal or larger than zmax in configuration file. Cannot proceed.")
            zRange = (round((zmin + zmax) / 2), range(zmin, zmax))

        if "laplacian" not in self.param.param["zProject"]["mode"]:
            self.data_2D = projectsImage2D(self.data, zRange, self.param.param["zProject"]["zProjectOption"])
            self.focusPlane = zRange[0]
            self.zRange = zRange[1]

        printLog("> Processing zRange:{}".format(self.zRange))

    # displays image and shows it
    def imageShow(
        self,
        show=False,
        cmap="plasma",
        size=(10, 10),
        dpi=300,
        outputName="tmp.png",
        save=True,
        normalization="stretch",
    ):
        fig = plt.figure()
        fig.set_size_inches(size)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()

        if normalization == "simple":
            norm = simple_norm(self.data_2D, "sqrt", percent=99.9)
        else:
            norm = ImageNormalize(stretch=SqrtStretch())

        ax.set_title("2D Data")

        if show:
            fig.add_axes(ax)
            ax.imshow(self.data_2D, origin="lower", cmap="Greys_r", norm=norm)
            return ax

        if save and not show:
            fig.add_axes(ax)
            ax.imshow(self.data_2D, origin="lower", cmap="Greys_r", norm=norm)
            fig.savefig(outputName)
            plt.close(fig)

    def removesBackground2D(self, normalize=False):
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(
            self.data_2D, (64, 64), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
        )

        im_bkg_substracted = self.data_2D - bkg.background

        if normalize:
            im_bkg_substracted = (im_bkg_substracted - im_bkg_substracted.min()) / (im_bkg_substracted.max())

        return im_bkg_substracted

    def imageShowWithValues(self, outputName):
        # _, mask = findsFocusFromBlocks(self.focalPlaneMatrix, self.LaplacianMatrix)

        imageShowWithValues(
            [self.focalPlaneMatrix],
            title="focal plane = " + "{:.2f}".format(self.focusPlane),
            outputName=outputName,
        )

# =============================================================================
# GENERAL FUNCTIONS
# =============================================================================

def fit1DGaussian_scipy(x,y,title='',verbose=False):
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
    dict()
        dictionary with fitting parameters.
    fig
        matplotlib figure for saving

    """
    fitResult=dict()

    try:
        fitgauss = spo.curve_fit(gaussian, x, y)
        fitResult["gauss1d.pos"] = fitgauss[0][1]
        fitResult["gauss1d.ampl"] = fitgauss[0][0]
        fitResult["gauss1d.fwhm"] = 2.355*fitgauss[0][2]
    except RuntimeError:
        return dict(), []
    except ValueError:    
        fitResult["gauss1d.pos"] = np.mean(x)
        fitResult["gauss1d.ampl"] = 0.0
        fitResult["gauss1d.fwhm"] = 0.0
        printLog("Warning: returned middle plane!")
        return fitResult, []
    
    if verbose:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)

        printLog("<<Fitting successful>>")

        ax.plot(x,y,'ko',label='data')
        ax.plot(x,gaussian(x,fitgauss[0][0],fitgauss[0][1],fitgauss[0][2]),linewidth=2, label='gaussian fit')
        ax.legend(loc=2)
        ax.set_title(title)
        return fitResult, fig

    return fitResult, []

def makesShiftMatrixHiRes(shiftMatrices, block_ref_shape):
    """
    Reinterpolates a block matrix to the full size of a larger image

    Parameters
    ----------
    shiftMatrices : list of numpy arrays
        list containing block matrices.
    block_ref_shape : tuple
        shape of block matrix.

    Returns
    -------
    shiftMatrix : numpy array
        Reinterpolated (larger) image.
        Size will be N x n x n
        where n is block_ref_shape[3]
        and N is len(shiftMatrices)

    """
    numberBlocks = block_ref_shape[0]
    blockSizeXY = block_ref_shape[3]

    shiftMatrix=np.zeros((len(shiftMatrices),blockSizeXY*shiftMatrices[0].shape[0],blockSizeXY*shiftMatrices[0].shape[1]))
    for _ax,m in enumerate(shiftMatrices):
        # printLog("size={}".format(m.shape))
        for i in range(numberBlocks):
            for j in range(numberBlocks):
                shiftMatrix[_ax,i * blockSizeXY: (i + 1) * blockSizeXY,j * blockSizeXY: (j + 1) * blockSizeXY] = m[i,j]
    return shiftMatrix


def projectsImage2D(img, zRange, mode):

    # sums images
    imageSize = img.shape
    I_collapsed = np.zeros((imageSize[1], imageSize[2]))

    if "MIP" in mode:
        # Max projection of selected planes
        I_collapsed = np.max(img[zRange[1][0] : zRange[1][-1]], axis=0)
    elif "sum" in mode:
        # Sums selected planes
        for i in zRange[1]:
            I_collapsed += img[i]
    else:
        printLog("ERROR: mode not recognized. Expected: MIP or sum. Read: {}".format(mode))

    return I_collapsed


# Gaussian function
# @jit(nopython=True)
def gaussian(x, a=1, mean=0, std=0.5):
    return a * (1 / (std * (np.sqrt(2 * np.pi)))) * (np.exp(-((x - mean) ** 2) / ((2 * std) ** 2)))


# Finds best focal plane by determining the max of the std deviation vs z curve

# @jit(nopython=True)
def calculate_zrange(idata, parameters):
    """
    Calculates the focal planes based max standard deviation
    """
    zwin = parameters.param["zProject"]["zwindows"]
    numPlanes = parameters.param["zProject"]["zmax"] - parameters.param["zProject"]["zmin"]
    stdMatrix = np.zeros(numPlanes)
    meanMatrix = np.zeros(numPlanes)

    # calculate STD in each plane
    for i in range(0, numPlanes):
        stdMatrix[i] = np.std(idata[i])
        meanMatrix[i] = np.mean(idata[i])

    maxStd = np.max(stdMatrix)
    ifocusPlane = np.where(stdMatrix == maxStd)[0][0]
    # Select a window to avoid being on the edges of the stack

    if ifocusPlane < parameters.param["zProject"]["windowSecurity"] or (
        ifocusPlane > numPlanes - parameters.param["zProject"]["windowSecurity"]
    ):
        focusPlane = ifocusPlane
    else:
        # interpolate zfocus
        axisZ = range(
            max(
                parameters.param["zProject"]["zmin"],
                ifocusPlane - parameters.param["zProject"]["windowSecurity"],
                min(
                    parameters.param["zProject"]["zmax"], ifocusPlane + parameters.param["zProject"]["windowSecurity"],
                ),
            )
        )

        stdMatrix -= np.min(stdMatrix)
        stdMatrix /= np.max(stdMatrix)

        try:
            fitgauss = spo.curve_fit(gaussian, axisZ, stdMatrix[axisZ[0] : axisZ[-1] + 1])
            # printLog("Estimation of focal plane (px): ", int(fitgauss[0][1]))
            focusPlane = int(fitgauss[0][1])
        except RuntimeError:
            printLog("Warning, too many iterations")
            focusPlane = ifocusPlane

    zmin = max(parameters.param["zProject"]["windowSecurity"], focusPlane - parameters.param["zProject"]["zwindows"],)
    zmax = min(
        numPlanes,
        parameters.param["zProject"]["windowSecurity"] + numPlanes,
        focusPlane + parameters.param["zProject"]["zwindows"],
    )
    zrange = range(zmin, zmax + 1)

    return focusPlane, zrange

def find_transform(im_src, im_dst):
    warp = np.eye(3, dtype=np.float32)
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 0.001)
    try:
        _, warp = cv2.findTransformECC(im_src, im_dst, warp, cv2.MOTION_HOMOGRAPHY, criteria)
    except:
        printLog("Warning: find transform failed. Set warp as identity")
    return warp

def variance_of_laplacian(image):
    # compute the Laplacian of the image and then return the focus
    # measure, which is simply the variance of the Laplacian
    return cv2.Laplacian(image, cv2.CV_64F).var()

def scatters3Dimage(client,image):
    """
    splits 3D image plane by plane and scatteres them to a cluster

    Parameters
    ----------
    client : dask CLient()
        result of get_client()
    image : numpy array
        3D image.

    Returns
    -------
    imageListScattered :  List, dict, iterator, or queue of futures matching the type of input.
        scattered image.

    """
    numberPlanes = image.shape[0]
    imageListScattered = [image[z, :, :] for z in range(numberPlanes)]
    # imageListScattered = client.scatter(imageListScattered)
    return imageListScattered

def reassembles3Dimage(client,futures,output_shape):
    """
    waits for futures to arrive
    collects them into a results list
    reassembles 3D image plane by plane

    Parameters
    ----------
    client : dask CLient()
        result of get_client()
    futures : list()
        list of futures
    output_shape : tuple
        result of image.shape

    Returns
    -------
    output : numpy array
        contains reassembled 3D image.

    """
    results = client.gather(futures)
    printLog(" > Retrieving {} results from cluster".format(len(results)))

    output = np.zeros(output_shape)
    for z, result in enumerate(results):
        output[z, :, :] = result

    del results
    return output

# =============================================================================
# CONTRAST and PIXEL INTENSITY NORMALIZATION, INHOMOGENEOUS BACKGROUND
# =============================================================================

def preProcess3DImage(x,lower_threshold, higher_threshold, parallelExecution=True):
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

    # printLog("Removing inhomogeneous background...")
    image = _removesInhomogeneousBackground(image,parallelExecution=parallelExecution)

    # printLog("Rescaling grey levels...")
    image = imageAdjust(image, lower_threshold=lower_threshold, higher_threshold=higher_threshold)[0]

    return image

def imageAdjust(image, lower_threshold=0.3, higher_threshold=0.9999):
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
    # printLog("> Rescaling grey levels...")

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
    image1 = exposure.rescale_intensity(image1, in_range=(lower_cutoff, higher_cutoff), out_range=(0, 1))

    # calculates histogram of intensities of adjusted image
    hist1 = exposure.histogram(image1)

    return image1, hist1_before, hist1, lower_cutoff, higher_cutoff

def saveImageDifferences(I1, I2, I3, I4, outputFileName):
    """
    Overlays two images as R and B and saves them to output file
    """

    I1, I2 = I1 / I1.max(), I2 / I2.max()
    I3, I4 = I3 / I3.max(), I4 / I4.max()

    I1, _, _, _, _ = imageAdjust(I1, lower_threshold=0.5, higher_threshold=0.9999)
    I2, _, _, _, _ = imageAdjust(I2, lower_threshold=0.5, higher_threshold=0.9999)
    I3, _, _, _, _ = imageAdjust(I3, lower_threshold=0.5, higher_threshold=0.9999)
    I4, _, _, _, _ = imageAdjust(I4, lower_threshold=0.5, higher_threshold=0.9999)

    cmap = "seismic"

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((60, 30))

    ax1.imshow(I1 - I2, cmap=cmap)
    ax1.axis("off")
    ax1.set_title("uncorrected")

    ax2.imshow(I3 - I4, cmap=cmap)
    ax2.axis("off")
    ax2.set_title("corrected")

    fig.savefig(outputFileName)

    plt.close(fig)


def _removesInhomogeneousBackground(im, boxSize=(32, 32), filter_size=(3, 3),verbose=True,parallelExecution=True,background=False):
    """
    wrapper to remove inhomogeneous backgrounds for 2D and 3D images

    Parameters
    ----------
    im : numpy array
        input image.
    boxSize : tuple of ints, optional
        size of boxSize used for block decomposition. The default is (32, 32).
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
        output = _removesInhomogeneousBackground2D(im, boxSize=(32, 32), filter_size=(3, 3),verbose=verbose, background=background)
    elif len(im.shape) == 3:
        output = _removesInhomogeneousBackground3D(im, boxSize=(32, 32), filter_size=(3, 3),verbose=verbose,parallelExecution=parallelExecution,background=background)
    else:
        return None
    return output

def _removesInhomogeneousBackground2D(im, boxSize=(32, 32), filter_size=(3, 3), background = False, verbose=True):
    """
    Calls Background2D() from ASTROPY to perform background substraction in a 2D image

    Parameters
    ----------
    im : numpy array
        2D input image.
    boxSize : tuple of ints, optional
        size of boxSize used for block decomposition. The default is (32, 32).
    filter_size : tuple of ints, optional
        Size of gaussian filter used for smoothing results. The default is (3, 3).
    background : boolean, optional
        if True returs the substracted image and the background
        otherwise only the former. The default is False.
    verbose : boolean, optional
        The default is True.

    Returns
    -------
    numpy array
        background substracted 2D image.

    """
    printLog("Removing inhomogeneous background from 2D image...")

    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    bkg = Background2D(im, (64, 64), filter_size=filter_size, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,)

    im1_bkg_substracted = im - bkg.background

    if background:
        return im1_bkg_substracted, bkg
    else:
        return im1_bkg_substracted


def _removesInhomogeneousBackground3D(image3D, boxSize=(64, 64), filter_size=(3, 3),verbose=True, parallelExecution=True, background = False):
    """
    Wrapper to remove inhomogeneous background in a 3D image by recursively calling _removesInhomogeneousBackground2D():
        - addresses output
        - iterates over planes and calls _removesInhomogeneousBackground2D in each plane
        - reassembles results into a 3D image

    Parameters
    ----------
    image3D : numpy array
        input 3D image.
    boxSize : tuple of ints, optional
        size of boxSize used for block decomposition. The default is (32, 32).
    filter_size : tuple of ints, optional
        Size of gaussian filter used for smoothing results. The default is (3, 3).
    verbose : boolean, optional
        The default is True.

    Returns
    -------
    output : numpy array
        processed 3D image.

    """
    if parallelExecution:
        client = try_get_client()
    else:
        client = None

    numberPlanes = image3D.shape[0]
    output = np.zeros(image3D.shape)

    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    if client is not None:
        printLog("> Removing inhomogeneous background from {} planes using {} workers...".format(numberPlanes,len(client.scheduler_info()['workers'])))
        imageList = [image3D[z, :, :] for z in range(numberPlanes)]
        # imageListScattered = client.scatter(imageList)

        futures = [client.submit(Background2D,img,boxSize,
                                    filter_size=filter_size, sigma_clip=sigma_clip,
                                    bkg_estimator=bkg_estimator) for img in imageList]

        results = client.gather(futures)
        printLog(" > Retrieving {} results from cluster".format(len(results)))

        for z, img, bkg in zip(range(numberPlanes),imageList,results):
            output[z, :, :] = img - bkg.background
        del results, futures, imageList
        # del imageListScattered

    else:
        printLog("> Removing inhomogeneous background from {} planes using 1 worker...".format(numberPlanes))
        Zrange=trange(numberPlanes)
        for z in Zrange:
            image2D = image3D[z, :, :]
            bkg = Background2D(
                image2D, boxSize, filter_size=filter_size, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
            )
            output[z, :, :] = image2D - bkg.background

    if background:
        return output, bkg.background
    else:
        return output
# =============================================================================
# IMAGE ALIGNMENT
# =============================================================================

def appliesXYshift3Dimages(image, shift,parallelExecution=True):
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
    if parallelExecution:
        client = try_get_client()
    else:
        client = None

    numberPlanes = image.shape[0]

    if client is None:
        printLog("> Shifting {} planes with 1 thread...".format(numberPlanes))
        shift3D = np.zeros((3))
        shift3D[0], shift3D[1], shift3D[2] = 0, shift[0], shift[1]
        output = shiftImage(image, shift3D)
    else:
        printLog("> Shifting {} planes using {} workers...".format(numberPlanes,len(client.scheduler_info()['workers'])))

        imageListScattered = scatters3Dimage(client,image)

        futures = [client.submit(shiftImage,img, shift) for img in imageListScattered]

        output = reassembles3Dimage(client,futures,image.shape)

        del futures
        del imageListScattered

    printLog("$ Done shifting 3D image.")

    return output

def imageBlockAlignment3D(images, blockSizeXY=256, upsample_factor=100):

    # sanity checks
    if len(images) < 2:
        sys.exit("# Error, number of images must be 2, not {}".format(len(images)))

    # - break in blocks
    numPlanes = images[0].shape[0]
    blockSize = (numPlanes, blockSizeXY, blockSizeXY)

    printLog("$ Breaking images into blocks")
    blocks = [view_as_blocks(x, block_shape=blockSize).squeeze() for x in images]

    block_ref = blocks[0]
    block_target = blocks[1]

    # - loop thru blocks and calculates block shift in xyz:
    shiftMatrices = [np.zeros(block_ref.shape[0:2]) for x in range(3)]

    # printLog("$ Aligning {} blocks".format(len(block_ref.shape[0])))
    for i in trange(block_ref.shape[0]):
        for j in range(block_ref.shape[1]):
            # - cross correlate in 3D to find 3D shift
            shifts_xyz, _, _ = phase_cross_correlation(
                block_ref[i, j], block_target[i, j], upsample_factor=upsample_factor
            )
            for matrix, _shift in zip(shiftMatrices, shifts_xyz):
                matrix[i, j] = _shift

    return shiftMatrices, block_ref, block_target

def combinesBlocksImageByReprojection(block_ref, block_target, shiftMatrices=None, axis1=0):
    """
    This routine will overlap block_ref and block_target images block by block.
    block_ref will be used as a template.
    - block_target will be first translated in ZXY using the corresponding values in shiftMatrices
    to realign each block
    - then an RGB image will be created with block_ref in the red channel, and the reinterpolated
    block_target block in the green channel.
    - the Blue channel is used for the grid to improve visualization of blocks.


    Parameters
    ----------
    block_ref : npy array
        return of view_as_blocks()
    block_target : npy array
        return of view_as_blocks()
    shiftMatrices : list of npy arrays
        index 0 contains Z, index 1 X and index 2 Y
    axis1 : int
        axis used for the reprojection: The default is 0.
        - 0 means an XY projection
        - 1 an ZX projection
        - 2 an ZY projection

    Returns
    -------
    output : NPY array of size imSize x imSize x 3
        RGB image.
    SSIM_as_blocks = NPY array of size numberBlocks x numberBlocks
        Structural similarity index between ref and target blocks
    """
    numberBlocks = block_ref.shape[0]
    blockSizes = list(block_ref.shape[2:])
    blockSizes.pop(axis1)
    imSizes = [x * numberBlocks for x in blockSizes]

    # gets ranges for slicing
    sliceCoordinates = []
    for blockSize in blockSizes:
        sliceCoordinates.append([range(x * blockSize, (x + 1) * blockSize) for x in range(numberBlocks)])

    # creates output images
    output = np.zeros((imSizes[0], imSizes[1], 3))
    SSIM_as_blocks = np.zeros((numberBlocks, numberBlocks))
    MSE_as_blocks = np.zeros((numberBlocks, numberBlocks))
    NRMSE_as_blocks = np.zeros((numberBlocks, numberBlocks))

    # blank image for blue channel to show borders between blocks
    blue = np.zeros(blockSizes)
    blue[0, :], blue[:, 0], blue[:, -1], blue[-1, :] = [0.5] * 4

    # reassembles image
    # takes one plane block
    for i, iSlice in enumerate(tqdm(sliceCoordinates[0])):
        for j, jSlice in enumerate(sliceCoordinates[1]):
            imgs = list()
            imgs.append(block_ref[i, j]) # appends reference image to image list

            if shiftMatrices is not None:
                shift3D = np.array([x[i, j] for x in shiftMatrices]) # gets 3D shift from block decomposition
                imgs.append(shiftImage(block_target[i, j], shift3D)) # realigns and appends to image list
            else:
                imgs.append(block_target[i, j]) # appends original target with no re-alignment

            imgs = [np.sum(x, axis=axis1) for x in imgs] # projects along axis1
            imgs = [exposure.rescale_intensity(x, out_range=(0, 1)) for x in imgs] # rescales intensity values
            imgs = [imageAdjust(x, lower_threshold=0.5, higher_threshold=0.9999)[0] for x in imgs] # adjusts pixel intensities

            NRMSE_as_blocks[i,j] = normalized_root_mse(imgs[0], imgs[1],normalization = 'euclidean')
            MSE_as_blocks[i,j] = mean_squared_error(imgs[0], imgs[1])
            SSIM_as_blocks[i,j] = ssim(imgs[0], imgs[1], data_range=imgs[1].max() - imgs[1].min())

            imgs.append(blue) # appends last channel with grid

            RGB = np.dstack(imgs) # makes block RGB image

            output[iSlice[0] : iSlice[-1] + 1, jSlice[0] : jSlice[-1] + 1, :] = RGB # inserts block into final RGB stack

    return output, SSIM_as_blocks, MSE_as_blocks, NRMSE_as_blocks


def align2ImagesCrossCorrelation(
    image1_uncorrected, image2_uncorrected, lower_threshold=0.999, higher_threshold=0.9999999, upsample_factor=100
):
    """
    Aligns 2 images by contrast adjust and cross correlation
    Parameters
    ----------
    imReference : TYPE
        DESCRIPTION.
    Im2 : TYPE
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
    I_histogram : TYPE
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

    (image1_adjusted, hist1_before, hist1_after, lower_cutoff1, higher_cutoff1,) = imageAdjust(
        image1_uncorrected, lower_threshold=lower_threshold, higher_threshold=higher_threshold
    )
    (image2_adjusted, hist2_before, hist2_after, lower_cutoff2, higher_cutoff2,) = imageAdjust(
        image2_uncorrected, lower_threshold=lower_threshold, higher_threshold=higher_threshold
    )

    # zips histograms
    lower_threshold = {"Im1": lower_cutoff1, "Im2": lower_cutoff2}
    I_histogram = {
        "Im1": (hist1_before, hist1_after),
        "Im2": (hist2_before, hist2_after),
    }

    # calculates shift

    # shift, error, diffphase = register_translation(image1_adjusted, image2_adjusted, upsample_factor=upsample_factor)
    shift, error, diffphase = phase_cross_correlation(image1_adjusted, image2_adjusted, upsample_factor=upsample_factor)

    # corrects image
    # The shift corresponds to the pixel offset relative to the reference image
    image2_corrected = shiftImage(image2_adjusted, shift)
    image2_corrected = exposure.rescale_intensity(image2_corrected, out_range=(0, 1))

    results = (
        shift,
        error,
        diffphase,
        lower_threshold,
        I_histogram,
        image2_corrected,
        image1_adjusted,
        image2_adjusted,
    )

    return results

def alignCV2(im1, im2, warp_mode):

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
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations, termination_eps)

    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(
            im1, im2, warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1
        )
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(im1, im2, warp_matrix, warp_mode, criteria)
    except cv2.error:
        cc = 0
        # printLog('Warning: find transform failed. Set warp as identity')

    return cc, warp_matrix


def applyCorrection(im2, warp_matrix):

    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)

    return im2_aligned


def alignImagesByBlocks(
    I1,
    I2,
    blockSize,
    upsample_factor=100,
    minNumberPollsters=4,
    tolerance=0.1,
    useCV2=False,
    shiftErrorTolerance=5,
):

    Block1 = view_as_blocks(I1, blockSize)
    Block2 = view_as_blocks(I2, blockSize)

    if useCV2:
        warp_matrix = np.eye(2, 3, dtype=np.float32)
        warp_mode = cv2.MOTION_TRANSLATION

    shiftImageNorm = np.zeros((Block1.shape[0], Block1.shape[1]))
    shiftedImage = np.zeros((Block1.shape[0], Block1.shape[1], 2))
    rmsImage = np.zeros((Block1.shape[0], Block1.shape[1]))

    for i in trange(Block1.shape[0]):
        for j in range(Block1.shape[1]):
            if not useCV2:
                # using Scimage registration functions
                shift, error, diffphase = phase_cross_correlation(
                    Block1[i, j], Block2[i, j], upsample_factor=upsample_factor
                )
                shiftImageNorm[i, j] = LA.norm(shift)
                shiftedImage[i, j, 0], shiftedImage[i, j, 1] = shift[0], shift[1]
                I2_aligned = shiftImage(I2, shift)
            else:
                # uses CV2 cause it is 20 times faster than Scimage
                cc, warp_matrix = alignCV2(Block1[i, j], Block2[i, j], warp_mode)
                shiftImageNorm[i, j] = LA.norm(warp_matrix[:, 2])
                shiftedImage[i, j, 0], shiftedImage[i, j, 1] = warp_matrix[:, 2][0], warp_matrix[:, 2][1]
                I2_aligned = applyCorrection(I2, warp_matrix)

            rmsImage[i, j] = np.sum(np.sum(np.abs(I1 - I2_aligned), axis=1))

    # [calculates optimal shifts by polling blocks showing the best RMS]

    # threshold = filters.threshold_otsu(rmsImage)
    threshold = (1 + tolerance) * np.min(rmsImage)
    mask = rmsImage < threshold

    contours = measure.find_contours(rmsImage, threshold)

    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour = np.array([0, 0])

    # [Averages shifts and errors from regions within the tolerated blocks]
    meanShifts = [np.mean(shiftedImage[mask, 0]), np.mean(shiftedImage[mask, 1])]
    stdShifts = [np.std(shiftedImage[mask, 0]), np.std(shiftedImage[mask, 1])]
    meanShiftNorm = np.mean(shiftImageNorm[mask])
    meanError = np.mean(rmsImage[mask])
    relativeShifts = np.abs(shiftImageNorm - meanShiftNorm)

    # [calculates global shift, if it is better than the polled shift, or
    # if we do not have enough pollsters to fall back to then it does a global cross correlation!]
    meanShifts_global, _, _ = phase_cross_correlation(I1, I2, upsample_factor=100)
    I2_aligned_global = shiftImage(I2, shift)
    meanError_global = np.sum(np.sum(np.abs(I1 - I2_aligned_global), axis=1))

    printLog("Block alignment error: {}, global alignment error: {}".format(meanError, meanError_global))

    if np.sum(mask) < minNumberPollsters or meanError_global < meanError or np.max(stdShifts) > shiftErrorTolerance:
        meanShifts = meanShifts_global
        meanError = meanError_global
        printLog("Falling back to global registration")

    printLog("*** Global XY shifts: {:.2f} px | {:.2f} px".format(meanShifts_global[0], meanShifts_global[1]))
    printLog(
        "*** Mean polled XY shifts: {:.2f}({:.2f}) px | {:.2f}({:.2f}) px".format(
            meanShifts[0], stdShifts[0], meanShifts[1], stdShifts[1]
        )
    )

    return np.array(meanShifts), meanError, relativeShifts, rmsImage, contour

# =============================================================================
# FOCAL PLANE INTERPOLATION
# =============================================================================

def focalPlane(data,threshold_fwhm=20, verbose=False):
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
    focalPlane : float
        focal plane: max of fitted z-profile.
    fwhm : float
        full width hald maximum of fitted z-profile.

    """
    # finds focal plane
    rawImages=[data[i,:,:] for i in range(data.shape[0])]
    LaplacianVariance = [cv2.Laplacian(img, cv2.CV_64F).var() for img in rawImages]
    # LaplacianVariance = [0 if np.isnan(x) else x for x in LaplacianVariance] # removes any Nan that will trigger ValueError in spo.curve_fit
    # LaplacianVariance = [0 if np.isinf(x) else x for x in LaplacianVariance] # removes any Inf that will trigger ValueError in spo.curve_fit
    LaplacianVariance  = LaplacianVariance/max(LaplacianVariance)
    xCoord = range(len(LaplacianVariance))
    fitResult, _= fit1DGaussian_scipy(xCoord,LaplacianVariance,title='laplacian variance z-profile',verbose=verbose)

    if len(fitResult)>0:
        focalPlane = fitResult['gauss1d.pos']
        fwhm = fitResult['gauss1d.fwhm']

        if fwhm>threshold_fwhm:
            fwhm, focalPlane = np.nan, np.nan
    else:
        fwhm, focalPlane = np.nan, np.nan

    if verbose:
        printLog("Focal plane found: {} with fwhm: {}".format(focalPlane,fwhm))

    return focalPlane, fwhm

def calculatesFocusPerBlock(data, blockSizeXY=128):
    """
    Calculates the most likely focal plane of an image by breaking into blocks and calculating
    the focal plane in each block

    - breaks image into blocks
    - returns the focal plane + fwhm for each block

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    blockSizeXY : TYPE, optional
        DESCRIPTION. The default is 512.

    Returns
    -------
    focalPlaneMatrix: np array
        matrix containing the maximum of the laplacian variance per block
    fwhm: np array
        matrix containing the fwhm of the laplacian variance per block
    block: np array
        3D block reconstruction of matrix

    """
    nPlanes = data.shape[0]

    blockSize = (nPlanes, blockSizeXY, blockSizeXY)

    block = view_as_blocks(data, blockSize).squeeze()
    focalPlaneMatrix = np.zeros(block.shape[0:2])
    fwhm= np.zeros(block.shape[0:2])

    fwhm = dict()

    for i in trange(block.shape[0]):
        # fwhm[str(i)] = dict()
        for j in range(block.shape[1]):
            focalPlaneMatrix[i, j], fwhm[i,j] = focalPlane(block[i, j])

    return focalPlaneMatrix, fwhm, block

# Program to find most frequent
# element in a list
def most_frequent(List):
    return max(set(List), key = List.count)

def imReassemble(focalPlaneMatrix, block, window=0):
    """
    Makes 2D image from 3D stack by reassembling sub-blocks
    For each sub-block we know the optimal focal plane, which is
    selected for the assembly of the while image

    Parameters
    ----------
    focalPlaneMatrix : numpy 2D array
        matrix containing the focal plane selected for each block.
    block : numpy matrix
        original 3D image sorted by blocks.

    Returns
    -------
    output : numpy 2D array
        output 2D projection

    """
    # gets image size from block image
    numberBlocks = block.shape[0]
    blockSizeXY = block.shape[3]
    imSize = numberBlocks * blockSizeXY

    # gets ranges for slicing
    sliceCoordinates = [range(x * blockSizeXY, (x + 1) * blockSizeXY) for x in range(numberBlocks)]

    # creates output image
    output = np.zeros((imSize, imSize))

    # gets more common plane
    focalPlanes = list()
    for i, iSlice in enumerate(sliceCoordinates):
        for j, jSlice in enumerate(sliceCoordinates):
            focalPlanes.append(int(focalPlaneMatrix[i, j]))
    mostCommonFocalPlane = most_frequent(focalPlanes)

    # reassembles image
    if window == 0:
        # takes one plane block
        for i, iSlice in enumerate(sliceCoordinates):
            for j, jSlice in enumerate(sliceCoordinates):
                focus = int(focalPlaneMatrix[i, j])
                if np.abs(focus-mostCommonFocalPlane)>1:
                    focus = int(mostCommonFocalPlane)
                output[iSlice[0] : iSlice[-1] + 1, jSlice[0] : jSlice[-1] + 1] = block[i, j][focus, :, :]
    else:
        # takes neighboring planes by projecting
        for i, iSlice in enumerate(sliceCoordinates):
            for j, jSlice in enumerate(sliceCoordinates):
                focus = int(focalPlaneMatrix[i, j])
                if np.abs(focus-mostCommonFocalPlane)>1:
                    focus = int(mostCommonFocalPlane)
                zmin = np.max((0, focus - round(window / 2)))
                zmax = np.min((block[i, j].shape[0], focus + round(window / 2)))
                zRange = (focus, range(zmin, zmax))
                # printLog("zrange for ({},{})={}".format(i,j,zRange))
                output[iSlice[0] : iSlice[-1] + 1, jSlice[0] : jSlice[-1] + 1] = projectsImage2D(
                    block[i, j][:, :, :], zRange, "MIP"
                )

    return output

# def findsFocusFromBlocks(focalPlaneMatrix, LaplacianMeans, threshold=0.1):

#     # filters out region with low values of z and of laplacian variances
#     # as planes close to surface and with low variances tend to give problems
#     focalPlaneMatrixWeighted = LaplacianMeans * focalPlaneMatrix
#     mask = focalPlaneMatrixWeighted > threshold * focalPlaneMatrixWeighted.max()

#     # sets problematic blocks to zero
#     mask = focalPlaneMatrixWeighted
#     mask[focalPlaneMatrixWeighted < threshold * focalPlaneMatrixWeighted.max()] = 0
#     mask[focalPlaneMatrix < 10] = 0

#     # calculates focus from the good blocks
#     focus = np.mean(focalPlaneMatrix[mask > 0])

#     return focus, mask

def reinterpolatesFocalPlane(data, param):

    if "blockSize" in param["zProject"]:
        blockSizeXY = param["zProject"]["blockSize"]
    else:
        blockSizeXY = 128

    if "zwindows" in param["zProject"]:
        window = param["zProject"]["zwindows"]
    else:
        window = 0

    focalPlaneMatrix, zRange, block = _reinterpolatesFocalPlane(data, blockSizeXY=blockSizeXY, window=window)

    # reassembles image
    output = imReassemble(focalPlaneMatrix, block, window=window)

    return output,focalPlaneMatrix, zRange


def _reinterpolatesFocalPlane(data, blockSizeXY=256, window=10):
    """
    Reinterpolates the focal plane of a 3D image by breking it into blocks
    - Calculates the focalPlane and fwhm matrices by block
    - removes outliers
    - calculates the focal plane for each block using sigmaClip statistics
    - returns a tuple with focal plane and the range to use

    Parameters
    ----------
    data : numpy array
        input 3D image.
    blockSizeXY : int
        size of blocks in XY, typically 256.
    window : int, optional
        number of planes before and after the focal plane to construct the zRange. The default is 0.

    Returns
    -------
    focalPlaneMatrix : numpy array
        focal plane matrix.
    zRange : tuple
        focusPlane, zRange.
    block : numpy array
        block representation of 3D image.

    """
    printLog("> Reducing 3D image size by slicing around focal plane".format())
    # breaks into subplanes, iterates over them and calculates the focalPlane in each subplane.

    focalPlaneMatrix, fwhm, block= calculatesFocusPerBlock(data, blockSizeXY=blockSizeXY)

    focalPlanes2Process = focalPlaneMatrix[~np.isnan(focalPlaneMatrix)]

    focalPlane,_,_ = sigmaclip(focalPlanes2Process,high = 3, low=3)
    focusPlane = np.mean(focalPlane)
    if np.isnan(focusPlane):
        printLog("# focusPlane detection failed. Using full stack.")
        focusPlane = data.shape[0]//2
        zRange = focusPlane, range(0, data.shape[0])
    else:
        focusPlane = np.mean(focalPlane).astype('int64')
        zmin = np.max([focusPlane-window, 0])
        zmax = np.min([focusPlane+window, data.shape[0]])
        zRange = focusPlane, range(zmin,zmax)

    return focalPlaneMatrix, zRange, block

def _removesZplanes(image3D, Zrange):
    """
    Removes planes in input image.
    For instance, if you provide a Zrange = range(0,image3D.shape[0],2)

    then the routine will remove any other plane. Number of planes skipped
    can be programmed by tuning zRange.

    Parameters
    ----------
    image3D : numpy array
        input 3D image.
    Zrange : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    output = np.zeros((len(Zrange),image3D.shape[1],image3D.shape[2]))
    for i,index in enumerate(Zrange):
        output[i,:,:] = image3D[index,:,:]

    return output

def _averageZplanes(image3D, Zrange):
    """
    Removes z-planes by calculating the average between successive planes

    Parameters
    ----------
    image3D : numpy array
        input 3D image.
    Zrange : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """

    output = np.zeros((len(Zrange),image3D.shape[1],image3D.shape[2]))
    for i,index in enumerate(Zrange):
        av = (image3D[index, :, :].astype(np.float) + image3D[index+1, :, :].astype(np.float))/2
        output[i,:,:] = av.astype(np.uint16)

    return output

def _interpolatesZplanes(image3D, Zrange):
    """
    Removes z planes by reinterpolation

    TO BE CODED!

    Parameters
    ----------
    image3D : numpy array
        input 3D image.
    Zrange : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    from scipy.interpolate import interpn

    output = np.zeros((len(Zrange),image3D.shape[1],image3D.shape[2]))

    # need to code using interpn
    output = image3D

    return output

def reinterpolateZ(image3D, Zrange,mode='average'):
    """
    wrapper function for any kind of z-interpolation
    to reduce the number of planes in an image

    Parameters
    ----------
    image3D : numpy array
        input 3D image.
    Zrange : range
        range of planes for the output image.
    mode : str, optional
        'remove' will remove planes
        'interpolate' will perform an interpolation
        The default is 'remove'.

    Returns
    -------
    output: numpy array

    """
    if 'interpolate' in mode:
        output = _interpolatesZplanes(image3D, Zrange)
    elif 'remove' in mode:
        output = _removesZplanes(image3D, Zrange)
    elif 'average' in mode:
        output = _averageZplanes(image3D, Zrange)

    printLog("$ Reduced Z-planes from {} to {}".format(image3D.shape[0],output.shape[0]))

    return output

# =============================================================================
# SEGMENTATION FUNCTIONS
# =============================================================================

def _segments2DimageByThresholding(image2D,
                             threshold_over_std=10,
                             sigma = 3,
                             boxSize=(32, 32),
                             filter_size=(3, 3),
                             area_min = 3,
                             area_max=1000,
                             nlevels=64,
                             contrast=0.001,
                             deblend3D = False,
                             kernel = []):

    # makes threshold matrix
    threshold = np.zeros(image2D.shape)
    threshold[:]=threshold_over_std*image2D.max()/100

    # segments objects
    segm = detect_sources(image2D, threshold, npixels=area_min, filter_kernel=kernel,)

    if segm.nlabels>0:
        # removes masks too close to border
        segm.remove_border_labels(border_width=10)  # parameter to add to infoList

        if segm.nlabels>0:
            segm_deblend = deblend_sources(image2D,
                                            segm,
                                            npixels=area_min,  # watch out, this is per plane!
                                            filter_kernel=kernel,
                                            nlevels=nlevels,
                                            contrast=contrast,
                                            relabel=True,
                                            mode='exponential')
            if segm_deblend.nlabels>0:
                # removes Masks too big or too small
                for label in segm_deblend.labels:
                    # take regions with large enough areas
                    area = segm_deblend.get_area(label)
                    if area < area_min or area > area_max:
                        segm_deblend.remove_label(label=label)

                # relabel so masks numbers are consecutive
                # segm_deblend.relabel_consecutive()

            # image2DSegmented = segm.data % changed during recoding function
            image2DSegmented = segm_deblend.data

            image2DSegmented[image2DSegmented>0]=1
            return image2DSegmented

        else:
            # returns empty image as no objects were detected
            return segm.data

    else:
        # returns empty image as no objects were detected
        return segm.data



def _segments3Dvolumes_StarDist(image3D,
                                area_min = 3,
                                area_max=1000,
                                nlevels=64,
                                contrast=0.001,
                                deblend3D = False,
                                axis_norm=(0,1,2),
                                model_dir='/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/',
                                model_name='stardist_18032021_single_loci'):

    numberPlanes = image3D.shape[0]

    printLog("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    printLog("> Segmenting {} planes using 1 worker...".format(numberPlanes))
    printLog("> Loading model {} from {}...".format(model_name,model_dir))
    os.environ["CUDA_VISIBLE_DEVICES"]="1"

    model = StarDist3D(None, name=model_name, basedir=model_dir)
    limit_gpu_memory(None, allow_growth=True)

    im = normalize(image3D,1,99.8,axis=axis_norm)
    Lx = im.shape[1]

    if Lx<1000:
        labels, details = model.predict_instances(im)

    else:
        resizer = PadAndCropResizer()
        axes = 'ZYX'

        im = resizer.before(im, axes, model._axes_div_by(axes))
        labels, polys = model.predict_instances(im, n_tiles=(1,8,8))
        labels = resizer.after(labels, axes)

    mask = np.array(labels>0, dtype=int)

    # Now we want to separate objects in 3D using watersheding
    if deblend3D:
        labeledImage = _deblend3Dsegmentation(mask)
    else:
        labeledImage = labels

    return mask, labeledImage

def _segments3DvolumesByThresholding(image3D,
                                     threshold_over_std=10,
                                     sigma = 3,
                                     boxSize=(32, 32),
                                     filter_size=(3, 3),
                                     area_min = 3,
                                     area_max=1000,
                                     nlevels=64,
                                     contrast=0.001,
                                     deblend3D = False,
                                     parallelExecution=True):
    if parallelExecution:
        client = try_get_client()
    else:
        client = None

    numberPlanes = image3D.shape[0]

    kernel = Gaussian2DKernel(sigma, x_size=sigma, y_size=sigma)
    kernel.normalize()

    printLog("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

    parallel = True
    if client is None:
        parallel=False
    else:
        if len(client.scheduler_info()['workers'])<1:
            parallel = False
            printLog("# Failed getting workers. Report of scheduler:")
            for key in client.scheduler_info().keys():
                printLog("{}:{}".format(key, client.scheduler_info()[key]))

    if not parallel:
        printLog("> Segmenting {} planes using 1 worker...".format(numberPlanes))

        output = np.zeros(image3D.shape)

        for z in trange(numberPlanes):
            image2D = image3D[z, :, :]
            image2DSegmented = _segments2DimageByThresholding(image2D,
                                         threshold_over_std=threshold_over_std,
                                         sigma = sigma,
                                         boxSize=boxSize,
                                         filter_size=filter_size,
                                         area_min = area_min,
                                         area_max=area_max,
                                         nlevels=nlevels,
                                         contrast=contrast,
                                         deblend3D = deblend3D,
                                         kernel=kernel)
            output[z,:,:] = image2DSegmented

    else:

        printLog("> Segmenting {} planes using {} workers...".format(numberPlanes,len(client.scheduler_info()['workers'])))

        imageListScattered = scatters3Dimage(client,image3D)

        futures = [client.submit(_segments2DimageByThresholding,
                                         img,
                                         threshold_over_std=threshold_over_std,
                                         sigma = sigma,
                                         boxSize=boxSize,
                                         filter_size=filter_size,
                                         area_min = area_min,
                                         area_max=area_max,
                                         nlevels=nlevels,
                                         contrast=contrast,
                                         deblend3D = deblend3D,
                                         kernel=kernel) for img in imageListScattered]

        output = reassembles3Dimage(client,futures,image3D.shape)

        del futures, imageListScattered

    labels = measure.label(output)

    # Now we want to separate objects in 3D using watersheding
    if deblend3D:
        labels = _deblend3Dsegmentation(output)

    return output, labels


def _deblend3Dsegmentation(binary):
    '''
    Deblends objects in 3D using watersheding

    Parameters
    ----------
    binary : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    binary=binary>0
    printLog(" > Constructing distance matrix from 3D binary mask...")

    distance = apply_parallel(ndi.distance_transform_edt, binary)

    printLog(" > Deblending sources in 3D by watersheding...")
    coords = peak_local_max(distance, footprint=np.ones((10, 10, 25)), labels=binary)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)

    labels = watershed(-distance, markers, mask=binary)
    return labels

def _segments3DMasks(image3D,
                    axis_norm=(0,1,2),
                    pmin=1,
                    pmax=99.8,
                    model_dir="/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks",
                    model_name='stardist_20210625_deconvolved'):

    """
    Parameters
    ----------
    image3D : numpy ndarray (N-dimensional array)
        3D raw image to be segmented

    model_dir : List of strings, optional
        paths of all models directory, the default is ["/mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks"]

    model_name : List of strings, optional
        names of all models, the default is ['stardist_20210625_deconvolved']

    """

    np.random.seed(6)

    numberPlanes = image3D.shape[0]

    printLog("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    printLog("> Segmenting {} planes using 1 worker...".format(numberPlanes))
    printLog("> Loading model {} from {}...".format(model_name,model_dir))
    os.environ["CUDA_VISIBLE_DEVICES"]="1" # why do we need this?

    # Load the model
    # --------------

    model = StarDist3D(None, name=model_name, basedir=model_dir)
    limit_gpu_memory(None, allow_growth=True)

    im = normalize(image3D, pmin=pmin, pmax=pmax, axis=axis_norm)
    Lx = im.shape[1]

    if Lx < 351: # what is this value? should it be a k-arg?
        labels, details = model.predict_instances(im)

    else:
        resizer = PadAndCropResizer()
        axes = 'ZYX'

        im = resizer.before(im, axes, model._axes_div_by(axes))
        labels, polys = model.predict_instances(im, n_tiles=(1, 8, 8))
        labels = resizer.after(labels, axes)

    mask = np.array(labels > 0, dtype=int)
    mask[mask > 0] = 1

    return mask, labels


def plotRawImagesAndLabels(image,label, normalize = False, window = 3):

    """
    Parameters
    ----------
    image : List of numpy ndarray (N-dimensional array)
        3D raw image of format .tif

    label : List of numpy ndarray (N-dimensional array)
        3D labeled image of format .tif
    """
    from stardist import random_label_cmap
    cmap= random_label_cmap()

    moy = np.mean(image,axis=0)
    lbl_moy = np.max(label,axis=0)

    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches((50, 50))
    ax = axes.ravel()
    titles=['raw image','projected labeled image']

    ax[0].imshow(moy, cmap="Greys_r", origin="lower")

    ax[1].imshow(lbl_moy, cmap=cmap, origin="lower")

    for axis,title in zip(ax,titles):
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(title)

    return fig

########################################################
# SAVING ROUTINES
########################################################

def save2imagesRGB(I1, I2, outputFileName):
    """
    Overlays two images as R and B and saves them to output file
    """

    sz = I1.shape
    I1, I2 = I1 / I1.max(), I2 / I2.max()

    I1, _, _, _, _ = imageAdjust(I1, lower_threshold=0.5, higher_threshold=0.9999)
    I2, _, _, _, _ = imageAdjust(I2, lower_threshold=0.5, higher_threshold=0.9999)

    fig, ax1 = plt.subplots()
    fig.set_size_inches((30, 30))

    nullImage = np.zeros(sz)

    RGB = np.dstack([I1, I2, nullImage])
    ax1.imshow(RGB)
    ax1.axis("off")

    fig.savefig(outputFileName)

    plt.close(fig)

def saveImage2Dcmd(image, fileName):
    if image.shape > (1, 1):
        np.save(fileName, image)
        # log.report("Saving 2d projection to disk:{}\n".format(os.path.basename(fileName)),'info')
        printLog("$ Image saved to disk: {}".format(fileName + ".npy"), "info")
    else:
        printLog("# Warning, image is empty", "Warning")

def savesImageAsBlocks(img,fullFileName,blockSizeXY=256,label = 'rawImage'):
    numPlanes = img.shape[0]
    blockSize = (numPlanes, blockSizeXY, blockSizeXY)
    blocks = view_as_blocks(img, block_shape=blockSize).squeeze()
    printLog("\nDecomposing image into {} blocks".format(blocks.shape[0]*blocks.shape[1]))

    folder = fullFileName.split('.')[0]
    fileName = os.path.basename(fullFileName).split('.')[0]

    if not os.path.exists(folder):
        os.mkdir(folder)
        printLog("Folder created: {}".format(folder))

    for i in trange(blocks.shape[0]):
        for j in range(blocks.shape[1]):
            outfile = folder + os.sep+ fileName + "_" + label + "_block_" + str(i) + "_" + str(j) + ".tif"
            imsave(outfile, blocks[i, j])

    # cmap = "seismic"

def imageShowWithValuesSingle(ax, matrix, cbarlabel, fontsize, cbar_kw, valfmt="{x:.0f}", cmap="YlGn"):
    Row = ["".format(x) for x in range(matrix.shape[0])]
    im, cbar = heatmap(matrix, Row, Row, ax=ax, cmap=cmap, cbarlabel=cbarlabel, fontsize=fontsize, cbar_kw=cbar_kw)
    _ = annotate_heatmap(im, valfmt=valfmt, size=fontsize, threshold=None, textcolors=("black", "white")) #, fontsize=fontsize


def imageShowWithValues(matrices, outputName="tmp.png", cbarlabels=["focalPlane"], fontsize=6, verbose=False, title=""):
    """
    Plots a list of matrices with their values in each pixel.

    Parameters
    ----------
    matrices : list
        matrices to plot. Should be 2D numpy arrays
    outputName : TYPE, optional
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
    numberImages = len(matrices)
    fig, axes = plt.subplots(1, numberImages)
    fig.set_size_inches((numberImages*2, 5))
    ax = axes.ravel()
    fig.suptitle(title)
    cbar_kw = {}
    cbar_kw["fraction"] = 0.046
    cbar_kw["pad"] = 0.04

    if len(cbarlabels)!=numberImages:
        cbarlabels=cbarlabels[0]*numberImages

    for matrix,axis,cbarlabel in zip(matrices,ax,cbarlabels):
        imageShowWithValuesSingle(axis, matrix, cbarlabel, fontsize, cbar_kw)

        # imageShowWithValuesSingle(ax[1], matrices[1], "filtered laplacian*focalPlane", fontsize, cbar_kw)

    fig.tight_layout()
    plt.savefig(outputName)

    if not verbose:
        plt.close(fig)


def display3D(image3D = None,labels=None, localizationsList = None,z=40, rangeXY=1000, norm=True, cmap='Greys'):


    if image3D is not None:
        images = list()
        images.append(image3D[z,:,:])
        images.append(image3D[:,rangeXY,:])
        images.append(image3D[:,:,rangeXY])
    else:
        images=[1,1,1]

    if labels is not None:
        segmented = list()
        segmented.append(labels[z,:,:])
        segmented.append(labels[:,rangeXY,:])
        segmented.append(labels[:,:,rangeXY])
    else:
        segmented=[1,1,1]

    if localizationsList is not None:
        localizedList = list()

        for localizations in localizationsList:
            localized = list()
            localized.append(localizations[:,[2,1]])
            localized.append(localizations[:,[2,0]])
            localized.append(localizations[:,[1,0]])
            localizedList.append(localized)

    else:
        localizedList=[1,1,1]

    percent=99.5
    symbols=['+','o','*','^']
    colors=['r','b','g','y']

    fig, axes = plt.subplots(1, len(images))
    fig.set_size_inches(len(images) * 50, 50)
    ax = axes.ravel()

    for image,segm,axis, iPlane in zip(images,segmented, ax, range(len(ax))):
        if image3D is not None:
            if norm:
                norm = simple_norm(image, "sqrt", percent=percent)
                axis.imshow(image, cmap=cmap, origin="lower", norm=norm)
            else:
                axis.imshow(image, cmap=cmap, origin="lower")
        if labels is not None:
            axis.imshow(color.label2rgb(segm, bg_label=0),alpha=.3)
        if localizationsList is not None:
            for iLocList, symbol, Color in zip(range(len(localizedList)),symbols,colors):
                locs =  localizedList[iLocList][iPlane]
                axis.plot(locs[:,0],locs[:,1],symbol,color=Color, alpha=.7)

    return fig

def display3D_assembled(images, localizations = None, plottingRange = None, normalize=True, masks= None):
    wspace=25

    rows_XY, cols_XY = images[0].shape[1], images[0].shape[0]
    rows_YZ, cols_ZX = images[2].shape[0], images[1].shape[0]
    rows = rows_XY + rows_YZ + wspace
    cols = cols_XY + cols_ZX + wspace

    fig = plt.figure(figsize=(50,50), tight_layout=True)
    axis = fig.add_axes([0.1, 0.1, .8, .8])

    displayImage = np.zeros((rows,cols))
    displayImage[0:rows_XY,0:cols_XY] = images[0]
    displayImage[rows_XY+wspace:rows_XY+rows_YZ+wspace,0:cols_XY] = images[2]
    displayImage[0:rows_XY,cols_XY+wspace:cols_XY+cols_ZX+wspace] = images[1].transpose()

    if normalize:
        norm = simple_norm(displayImage[:,:], "sqrt", percent=99)
        axis.imshow(displayImage[:,:],cmap='Greys',alpha=1, norm=norm)
    else:
        axis.imshow(displayImage[:,:],cmap='Greys',alpha=1
                    )

    if masks is not None:
        displayMask = np.zeros((rows,cols))
        displayMask[0:rows_XY,0:cols_XY] = masks[0]
        displayMask[rows_XY+wspace:rows_XY+rows_YZ+wspace,0:cols_XY] = masks[2]
        displayMask[0:rows_XY,cols_XY+wspace:cols_XY+cols_ZX+wspace] = masks[1].transpose()

        axis.imshow(color.label2rgb(displayMask[:,:], bg_label=0),alpha=.3)

    colors = ['r','g','b','y']
    markersizes = [2,1,1,1,1]
    if localizations is not None:
        for i,loc in enumerate(localizations):
            axis.plot(loc[:,2],loc[:,1],'+',color=colors[i],markersize=1)
            if plottingRange is not None:
                selections = [np.abs(loc[:,a]-plottingRange[0])<plottingRange[1] for a in [2,1]]
                axis.plot(loc[selections[0],0]+cols_XY+wspace,loc[selections[0],1],'+',color=colors[i],markersize=markersizes[i],alpha=.9)
                axis.plot(loc[selections[1],2],loc[selections[1],0]+rows_XY+wspace,'+',color=colors[i],markersize=markersizes[i],alpha=.9)
            else:
                axis.plot(loc[:,0]+cols_XY,loc[:,1],'+',color=colors[i],markersize=1)

    axis.axes.yaxis.set_visible(False)
    axis.axes.xaxis.set_visible(False)
    axis.axis("off")

    return fig


def _plotsImage3D(image3D,localizations=None,masks=None,normalize=True,window = 3):
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

    images = list()
    images.append(np.sum(img,axis=0))
    images.append(np.sum(img[:,:,center-window:center+window],axis=2))
    images.append(np.sum(img[:,center-window:center+window,:],axis=1))

    if masks is not None:
        labels = list()
        labels.append(np.max(masks,axis=0))
        labels.append(np.max(masks[:,:,center-window:center+window],axis=2))
        labels.append(np.max(masks[:,center-window:center+window,:],axis=1))
    else:
        labels=None

    fig1 = display3D_assembled(images, localizations = localizations, masks = labels, plottingRange = [center,window],normalize=normalize)

    return fig1


def plots3DshiftMatrices(shiftMatrices, fontsize=8, log=False,valfmt="{x:.1f}"):

    cbar_kw = {}
    cbar_kw["fraction"] = 0.046
    cbar_kw["pad"] = 0.04

    fig, axes = plt.subplots(1, len(shiftMatrices))
    fig.set_size_inches((len(shiftMatrices) * 10, 10))
    ax = axes.ravel()
    titles = ["z shift matrix", "x shift matrix", "y shift matrix"]

    for axis, title, x in zip(ax, titles, shiftMatrices):
        if log:
            x=np.log10(x)
        imageShowWithValuesSingle(axis, x, title, fontsize, cbar_kw, valfmt=valfmt, cmap="YlGn")  # YlGnBu
        axis.set_title(title)

    return fig


def plots4images(allimages, titles=['reference','cycle <i>','processed reference','processed cycle <i>']):

    fig, axes = plt.subplots(2,2)
    fig.set_size_inches((10, 10))
    ax = axes.ravel()

    for axis, img, title in zip(ax, allimages,titles):
        # im = np.sum(img, axis=0)
        axis.imshow(img, cmap="Greys")
        axis.set_title(title)
    fig.tight_layout()

    return fig

def plottingBlockALignmentResults(relativeShifts, rmsImage, contour, fileName="BlockALignmentResults.png"):

    # plotting
    fig, axes = plt.subplots(1, 2)
    ax = axes.ravel()
    fig.set_size_inches((10, 5))

    cbwindow = 3
    p1 = ax[0].imshow(relativeShifts, cmap="terrain", vmin=0, vmax=cbwindow)
    ax[0].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[0].set_title("abs(global-block) shifts, px")
    fig.colorbar(p1, ax=ax[0], fraction=0.046, pad=0.04)

    p2 = ax[1].imshow(rmsImage, cmap="terrain", vmin=np.min(rmsImage), vmax=np.max(rmsImage))
    ax[1].plot(contour.T[1], contour.T[0], linewidth=2, c="k")
    ax[1].set_title("RMS")
    fig.colorbar(p2, ax=ax[1], fraction=0.046, pad=0.04)

    for x in range(len(ax)):
        ax[x].axis("off")

    fig.savefig(fileName)

    plt.close(fig)

def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", fontsize=12, **kwargs):
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
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.1f}", textcolors=("black", "white"), threshold=None, **textkw):
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
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
