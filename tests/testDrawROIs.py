#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:21:14 2020

@author: marcnol

test ROI drawing functionalities

This was written for a user to manually define RNA masks

"""

#%%


import logging
import numpy as np
from matplotlib import pyplot as plt
from roipoly import MultiRoi
from imageProcessing import Image

"""
logging.basicConfig(format='%(levelname)s ''%(processName)-10s : %(asctime)s '
                            '%(module)s.%(funcName)s:%(lineno)s %(message)s',
                    level=logging.INFO)
"""


def createsUserMask(filename):

    # loads image
    img = np.load(filename).squeeze()
    outputFileName = "masks"

    # displays image
    Im = Image()
    Im.data_2D = img
    Im.imageShow(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple ROIs
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    numberROIs = len(multiroi_named.rois)
    print("Number of ROIs drawn: {}".format(numberROIs))

    masks = np.zeros((img.shape[0], img.shape[1], numberROIs))

    # Display image with all ROIs
    Im.imageShow(show=True, normalization="simple")

    roi_names = []
    i = 0
    for name, roi in multiroi_named.rois.items():
        roi.display_roi()
        roi.display_mean(img)
        roi_names.append(name)
        masks[:, :, i] = roi.get_mask(img)
        i += 1
    plt.legend(roi_names, bbox_to_anchor=(1.2, 1.05))
    plt.show()

    # saves result
    np.save(outputFileName, masks)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # Create image
    rootFolder = "/home/marcnol/data/Experiment_4/0_Embryo/alignImages/"
    fileNameRNA = rootFolder + "scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy"

    createsUserMask(fileNameRNA)
