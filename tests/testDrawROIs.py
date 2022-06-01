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
    output_filename = "masks"

    # displays image
    im_obj = Image()
    im_obj.data_2d = img
    im_obj.show_image(show=True, normalization="simple")
    print("Click on the button to add a new ROI")

    # Draw multiple rois
    multiroi_named = MultiRoi(roi_names=["My first ROI", "My second ROI"])

    number_rois = len(multiroi_named.rois)
    print("Number of rois drawn: {}".format(number_rois))

    masks = np.zeros((img.shape[0], img.shape[1], number_rois))

    # Display image with all rois
    im_obj.show_image(show=True, normalization="simple")

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
    np.save(output_filename, masks)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # Create image
    root_folder = "/home/marcnol/data/Experiment_4/0_Embryo/alignImages/"
    fileNameRNA = root_folder + "scan_002_DAPI_001_ROI_converted_decon_ch01_2d_registered.npy"

    createsUserMask(fileNameRNA)
