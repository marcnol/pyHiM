#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 14:31:58 2022

@author: marcnol

this verifies the results of a 3D mask segmentation using stardist 3D networks

- loads 3D NPY array with labeled image and plots it alongside original image

"""

import os
import numpy as np
import napari
import matplotlib.pyplot as plt
import glob

folder_segmentedObjects = "segmentedObjects"
folder = "/home/marcnol/data/pancreas_Exp6_2022_002_ROI"+os.sep+folder_segmentedObjects

files = glob.glob(folder+os.sep+"*mask*3Dmasks*npy")

print(f"\n{len(files)} files found matching criteria:\n\n{files}")

#%% loads data

data = [np.load(x) for x in files]


#%% projects labeled image

numberObjects = [np.amax(x) for x in data]

data2D = [np.max(x,axis=0) for x in data]

plt.imshow(data2D[0], cmap="gist_earth")

for numberObject,file in zip(numberObjects,files):
    print(f"\nFile: {os.path.basename(file)}")
    print(f"Number of objects detected: {numberObject}")

#%% saves projections

output_files = [x.split('.')[0]+"_Masks.npy" for x in files]

print(f"output files: {output_files}\n\n")

for output_file,_data2D in zip(output_files,data2D):
    if os.path.exists(output_file):
        print(f"----Warning!----\nRenaming {output_file} as it exists already!\n")
        os.rename(output_file, output_file+"._2Dmasks.npy")

    print(f"> Saving {output_file}\n")
    np.save(output_file,_data2D)


#%% views results with Napari

data2D_reloaded = [np.load(x) for x in output_files]

viewer = napari.view_image(data2D_reloaded[0])



