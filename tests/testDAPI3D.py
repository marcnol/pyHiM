#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:10:17 2020

@author: marcnol
"""

import numpy as np
import os

from mayavi.mlab import *

from imageProcessing.imageProcessing import Image
from fileProcessing.fileManagement import Parameters, log
from skimage import io

x, y, z = np.ogrid[-10:10:20j, -10:10:20j, -10:10:20j]

s = np.sin(x*y*z)/(x*y*z)

# contour3d(s)

rootFolder = '/home/marcnol/Downloads'
fileName = 'scan_006_DAPI_001_ROI_converted_decon_ch00.tif'

fullFileName=rootFolder+os.sep+fileName

data = io.imread(fileName).squeeze()
subdata=data[22:45,0:2000,0:2000]
contour3d(subdata, vmin=5000, vmax=40000)
