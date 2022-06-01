#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:17:04 2020

@author: marcnol
"""


from dcimg import DCIMGFile
import os
from matplotlib.pylab import plt
import numpy as np

root_folder = "/mnt/grey/DATA/rawData_2020/Experiment_6_sara/dcimg_raw_data"
file_name = "185__FTL.dcimg"
fileID = DCIMGFile(root_folder + os.sep + file_name)

fileID.open()
image = fileID.mma

# image_2d=np.sum(image,axis=0)

plt.imshow(image[100, :, :])
# fileID.close()
