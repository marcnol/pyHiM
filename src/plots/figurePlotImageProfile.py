#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:11:33 2020

@author: marcnol


plot line profile from npy image

"""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


rootFolder = "/home/marcnol/data/Embryo_debug_dataset/Experiment_18/zProject"
file = "scan_001_RT29_001_ROI_converted_decon_ch01_2d.npy"
filename = rootFolder + os.sep + file

pixelSize = 0.1


xmin = 802
xmax = xmin
ymin = 935
ymax = 947
window = 5

im = np.load(filename)

# shows full image
figure, ax = plt.subplots(1)
im1 = im.copy()
im1 = im1 / 5
# im1=im1-im1.min()
# im1=1-im1
ax.imshow(im1, cmap="Greys", vmax=1000)
rect = Rectangle((xmax, ymin), 2 * window, ymax - ymin, edgecolor="b", facecolor="none")
ax.add_patch(rect)
# show crop

fig1 = plt.figure(constrained_layout=True)
crop = im[xmin - window : xmin + window, ymin:ymax]
crop = crop / crop.max()
crop = crop - crop.min()
crop = 1 - crop
plt.imshow(crop, cmap="Greys", vmax=1)

# shows profile
fig2 = plt.figure(constrained_layout=True)

profile = im[xmin, ymin:ymax]
profile = profile / profile.max()
# x = pixelSize*np.linspace(0, profile.shape[0], num=profile.shape[0], endpoint=True)
x = pixelSize * np.arange(profile.shape[0])
plt.plot(x, profile)
plt.plot(x, profile, "o")
