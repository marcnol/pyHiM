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
import numpy as np

from skimage import io
import scipy.optimize as spo
import matplotlib.pyplot as plt
import cv2
from numba import jit

from skimage import exposure
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch, simple_norm
# from skimage.feature import register_translation
from scipy.ndimage import shift as shiftImage

from skimage.util.shape import view_as_blocks
from numpy import linalg as LA
from tqdm import trange
from skimage import measure
from scipy.ndimage import shift as shiftImage
from skimage.exposure import match_histograms
from skimage.registration import phase_cross_correlation


from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

# =============================================================================
# CLASSES
# =============================================================================


class Image:
    def __init__(self,param,log1):
        self.param=param
        self.log=log1
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
    def saveImage2D(self, log, rootFolder, tag="_2d"):
        fileName=self.getImageFileName(rootFolder,tag)
        saveImage2Dcmd(self.data_2D, fileName, log)
        
    def getImageFileName(self,rootFolder,tag):
        fileName = rootFolder + os.sep + os.path.basename(self.fileName).split(".")[0] + tag
        return fileName
    
    # read an image as a numpy array
    def loadImage2D(self, fileName, log, masterFolder, tag="_2d"):
        self.fileName = fileName
        fileName=self.getImageFileName(masterFolder,tag)+ ".npy"

        self.data_2D = np.load(fileName)
        log.report(
            "\nLoading from disk:{}".format(os.path.basename(fileName)), "info",
        )

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
        # print("Image Name={}".format(self.fileName))
        self.log.report("Image Size={}".format(self.imageSize))
        # self.log.report("Stage position={}".format(self.stageCoordinates))
        self.log.report("Focal plane={}".format(self.focusPlane))

    # processes sum image in axial direction given range
    # @jit(nopython=True) 
    def zProjectionRange(self):

        # find the correct range for the projection
        if self.param.param["zProject"]["zmax"] > self.imageSize[0]:
            self.log.report("Setting z max to the last plane")
            self.param.param["zProject"]["zmax"] = self.imageSize[0]

        if self.param.param["zProject"]["mode"] == "automatic":
            print("Calculating planes...")
            zRange = calculate_zrange(self.data, self.param)
        elif self.param.param["zProject"]["mode"] == "full":
            (zmin, zmax) = (0, self.imageSize[0])
            zRange = (round((zmin + zmax) / 2), range(zmin, zmax))
        else:
            # Manual: reads from parameters file
            (zmin, zmax) = (
                self.param.param["zProject"]["zmin"],
                self.param.param["zProject"]["zmax"],
            )
            if zmin >= zmax:
                raise SystemExit('zmin is equal or larger than zmax in configuration file. Cannot proceed.')
            zRange = (round((zmin + zmax) / 2), range(zmin, zmax))

        self.log.report("Processing zRange:{}".format(zRange))

        # sums images
        I_collapsed = np.zeros((self.imageSize[1], self.imageSize[2]))
        if self.param.param["zProject"]["zProjectOption"] == "MIP":
            # Max projection of selected planes
            I_collapsed = np.max(self.data[zRange[1][0] : zRange[1][-1]], axis=0)
        else:
            # Sums selected planes
            for i in zRange[1]:
                I_collapsed += self.data[i]

        self.data_2D = I_collapsed
        self.zRange = zRange[1]
        self.focusPlane = zRange[0]

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
            # ax.imshow(self.data_2D, aspect='equal')
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


# =============================================================================
# FUNCTIONS
# =============================================================================


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
    """
    zmin = max(ifocusPlane - parameters.param['process']['windowSecurity'],
                          parameters.param['process']['zmin'])
    zmax = min(ifocusPlane + parameters.param['process']['windowSecurity'],
                          parameters.param['process']['zmax'])
    """
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
            # print("Estimation of focal plane (px): ", int(fitgauss[0][1]))
            focusPlane = int(fitgauss[0][1])
        except RuntimeError:
            print("Warning, too many iterations")
            focusPlane = ifocusPlane

    zmin = max(parameters.param["zProject"]["windowSecurity"], focusPlane - parameters.param["zProject"]["zwindows"],)
    zmax = min(
        numPlanes,
        parameters.param["zProject"]["windowSecurity"] + numPlanes,
        focusPlane + parameters.param["zProject"]["zwindows"],
    )
    zrange = range(zmin, zmax + 1)

    return focusPlane, zrange


def imageAdjust(image, lower_threshold=0.3, higher_threshold=0.9999):
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


def save2imagesRGB(I1, I2, outputFileName):
    """
    Overlays two images as R and B and saves them to output file
    """

    sz = I1.shape
    I1, I2 = I1 / I1.max(), I2 / I2.max()   

    I1,_,_,_,_ = imageAdjust(I1, lower_threshold=0.5, higher_threshold=0.9999)
    I2,_,_,_,_ = imageAdjust(I2, lower_threshold=0.5, higher_threshold=0.9999)
    
    fig, ax1 = plt.subplots()
    fig.set_size_inches((30, 30))
   
    nullImage = np.zeros(sz)
   
    RGB = np.dstack([I1, I2, nullImage])
    ax1.imshow(RGB)
    ax1.axis("off")
   
    fig.savefig(outputFileName)

    plt.close(fig)

def saveImageDifferences(I1, I2, I3, I4, outputFileName):
    """
    Overlays two images as R and B and saves them to output file
    """

    I1, I2 = I1 / I1.max(), I2 / I2.max()   
    I3, I4 = I3 / I3.max(), I4 / I4.max()   

    I1,_,_,_,_ = imageAdjust(I1, lower_threshold=0.5, higher_threshold=0.9999)
    I2,_,_,_,_ = imageAdjust(I2, lower_threshold=0.5, higher_threshold=0.9999)
    I3,_,_,_,_ = imageAdjust(I3, lower_threshold=0.5, higher_threshold=0.9999)
    I4,_,_,_,_ = imageAdjust(I4, lower_threshold=0.5, higher_threshold=0.9999)

    cmap = 'seismic'
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    fig.set_size_inches((60, 30))
    
    ax1.imshow(I1-I2, cmap=cmap)
    ax1.axis("off")
    ax1.set_title("uncorrected")
    
    ax2.imshow(I3-I4, cmap=cmap)
    ax2.axis("off")
    ax2.set_title("corrected")
    
    fig.savefig(outputFileName)

    plt.close(fig)



def saveImage2Dcmd(image, fileName, log):
    if image.shape > (1, 1):
        np.save(fileName, image)
        # log.report("Saving 2d projection to disk:{}\n".format(os.path.basename(fileName)),'info')
        log.report("Image saved to disk: {}".format(fileName + ".npy"), "info")
    else:
        log.report("Warning, image is empty", "Warning")


def align2ImagesCrossCorrelation(image1_uncorrected, 
                                 image2_uncorrected,
                                 lower_threshold=0.999, 
                                 higher_threshold=0.9999999,
                                 upsample_factor=100):
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

    results=(
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

def find_transform(im_src, im_dst):
    warp = np.eye(3, dtype=np.float32)
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 0.001)
    try:
        _, warp = cv2.findTransformECC(im_src, im_dst, warp, cv2.MOTION_HOMOGRAPHY, criteria)
    except:
        print('Warning: find transform failed. Set warp as identity')
    return warp
    
def alignCV2(im1,im2,warp_mode):
    
   
    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY :
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else :
        warp_matrix = np.eye(2, 3, dtype=np.float32)
    
    # Specify the number of iterations.
    number_of_iterations = 1000 # 5000
    
    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-10;
    
    # Define termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations,  termination_eps)
    
    # Run the ECC algorithm. The results are stored in warp_matrix.
    try:
        cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria, inputMask=None, gaussFiltSize=1)
    except TypeError:
        cc, warp_matrix = cv2.findTransformECC(im1,im2,warp_matrix, warp_mode, criteria)
    except cv2.error:
        cc=0
        # print('Warning: find transform failed. Set warp as identity')

    return cc, warp_matrix

def applyCorrection(im2,warp_matrix):
    
    sz = im2.shape

    # Use warpAffine for Translation, Euclidean and Affine
    im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);

    return im2_aligned 

   
def alignImagesByBlocks(I1, I2, blockSize, log1, upsample_factor=100, minNumberPollsters=4, tolerance=0.1, useCV2=False,shiftErrorTolerance = 5):

    Block1=view_as_blocks(I1,blockSize)
    Block2=view_as_blocks(I2,blockSize)
    
    if useCV2:
        warp_matrix = np.eye(2, 3, dtype=np.float32)
        warp_mode = cv2.MOTION_TRANSLATION
    
    shiftImageNorm= np.zeros((Block1.shape[0],Block1.shape[1]))
    shiftedImage = np.zeros((Block1.shape[0],Block1.shape[1],2))
    rmsImage= np.zeros((Block1.shape[0],Block1.shape[1]))
    
    for i in trange(Block1.shape[0]):
        for j in range(Block1.shape[1]):
            if not useCV2:
                # using Scimage registration functions
                shift, error, diffphase = phase_cross_correlation(Block1[i,j], Block2[i,j],upsample_factor=upsample_factor)
                shiftImageNorm[i,j] = LA.norm(shift)
                shiftedImage[i,j,0], shiftedImage[i,j,1] = shift[0], shift[1]
                I2_aligned = shiftImage(I2, shift)
            else:              
                # uses CV2 cause it is 20 times faster than Scimage
                cc, warp_matrix = alignCV2(Block1[i,j], Block2[i,j], warp_mode)
                shiftImageNorm[i,j] = LA.norm(warp_matrix[:,2])
                shiftedImage[i,j,0],shiftedImage[i,j,1] = warp_matrix[:,2][0], warp_matrix[:,2][1]
                I2_aligned = applyCorrection(I2,warp_matrix)
            
            rmsImage[i,j] =np.sum(np.sum(np.abs(I1-I2_aligned),axis=1)) 
            
    # [calculates optimal shifts by polling blocks showing the best RMS]

    # threshold = filters.threshold_otsu(rmsImage)
    threshold = (1+tolerance)*np.min(rmsImage)
    mask = rmsImage < threshold 
    
    contours = measure.find_contours(rmsImage, threshold)
    
    try:
        contour = sorted(contours, key=lambda x: len(x))[-1]
    except IndexError:
        contour=np.array([0,0])
            
    # [Averages shifts and errors from regions within the tolerated blocks]
    meanShifts = [np.mean(shiftedImage[mask,0]), np.mean(shiftedImage[mask,1])]
    stdShifts =[np.std(shiftedImage[mask,0]), np.std(shiftedImage[mask,1])]
    meanShiftNorm = np.mean(shiftImageNorm[mask])
    meanError = np.mean(rmsImage[mask])
    relativeShifts= np.abs(shiftImageNorm-meanShiftNorm)

    # [calculates global shift, if it is better than the polled shift, or
    # if we do not have enough pollsters to fall back to then it does a global cross correlation!]
    meanShifts_global, _, _= phase_cross_correlation(I1,I2,upsample_factor=100)
    I2_aligned_global  = shiftImage(I2, shift)
    meanError_global = np.sum(np.sum(np.abs(I1-I2_aligned_global),axis=1)) 

    log1.info("Block alignment error: {}, global alignment error: {}".format(meanError,meanError_global))

    if np.sum(mask)<minNumberPollsters or meanError_global<meanError or np.max(stdShifts)>shiftErrorTolerance:
        meanShifts = meanShifts_global             
        meanError=meanError_global
        log1.info("Falling back to global registration")

    log1.info("*** Global XY shifts: {:.2f} px | {:.2f} px".format(meanShifts_global[0], meanShifts_global[1]))            
    log1.info("*** Mean polled XY shifts: {:.2f}({:.2f}) px | {:.2f}({:.2f}) px".format(meanShifts[0], stdShifts[0], meanShifts[1], stdShifts[1]))
    
    return np.array(meanShifts), meanError, relativeShifts, rmsImage, contour

def plottingBlockALignmentResults(relativeShifts, rmsImage, contour, fileName='BlockALignmentResults.png'):    

    # plotting
    fig, axes = plt.subplots(1,2)
    ax=axes.ravel()
    fig.set_size_inches((10,5))
    
    cbwindow=3
    p1 = ax[0].imshow(relativeShifts,cmap='terrain',vmin=0, vmax=cbwindow)
    ax[0].plot(contour.T[1], contour.T[0], linewidth=2, c='k')
    ax[0].set_title("abs(global-block) shifts, px")
    fig.colorbar(p1,ax=ax[0],fraction=0.046, pad=0.04)
    
    p2 = ax[1].imshow(rmsImage,cmap='terrain',vmin=np.min(rmsImage), vmax=np.max(rmsImage))
    ax[1].plot(contour.T[1], contour.T[0], linewidth=2, c='k')
    ax[1].set_title("RMS")
    fig.colorbar(p2,ax=ax[1],fraction=0.046, pad=0.04)
    
    for x in range(len(ax)):
        ax[x].axis('off')
        
    fig.savefig(fileName)

    plt.close(fig)