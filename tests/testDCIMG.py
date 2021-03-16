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

rootFolder = "/mnt/grey/DATA/rawData_2020/Experiment_6_sara/dcimg_raw_data"
fileName = "185__FTL.dcimg"
fileID = DCIMGFile(rootFolder + os.sep + fileName)

fileID.open()
image = fileID.mma

# image2D=np.sum(image,axis=0)

plt.imshow(image[100, :, :])
# fileID.close()
