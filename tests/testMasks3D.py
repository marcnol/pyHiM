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
from tifffile import imsave

folder_segmentedObjects = "segmentedObjects"
folder = "/home/marcnol/data/pancreas_Exp6_2022_002_ROI"+os.sep+folder_segmentedObjects

files = glob.glob(folder+os.sep+"*mask*3Dmasks*npy")

print(f"\n{len(files)} files found matching criteria:\n\n{files}")

#%% loads data


data = [np.load(x) for x in files]

for file,im in zip(files,data):
    output_TIFF = file.split(".npy")[0]+".tif"
    print(f"Saving original NPY as TIFF:\n--> {output_TIFF}")
    imsave(output_TIFF,im)

#%% projects labeled image

numberObjects = [np.amax(x) for x in data]

data2D = [np.max(x,axis=0) for x in data]

for numberObject,file in zip(numberObjects,files):
    print(f"\nFile: {os.path.basename(file)}")
    print(f"Number of objects detected: {numberObject}")

#%% saves projections

output_files = [x.split('.npy')[0]+"_Masks.npy" for x in files]
output_files_TIFF = [x.split('.npy')[0]+"_Masks.tif" for x in files]

print(f"output files: {output_files}\n\n")

for output_file,_data2D, output_file_TIFF in zip(output_files,data2D,output_files_TIFF):
    if os.path.exists(output_file):
        print(f"----Warning!----\nRenaming {output_file} as it exists already!\n")
        os.rename(output_file, output_file+"._2Dmasks.npy")

    print(f"> Saving:\n--> {output_file}\n--> {output_file_TIFF} \n")
    np.save(output_file,_data2D)
    imsave(output_file_TIFF,_data2D)

#%% views results with Napari

data2D_reloaded = [np.load(x) for x in output_files]

viewer = napari.view_image(data2D_reloaded[0])

#%% converts original TIF by reducing planes by half

from tifffile import imread

def reinterpolateZ(image3D, Zrange):
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

    output = np.zeros((len(Zrange),image3D.shape[1],image3D.shape[2]))
    for i,index in enumerate(Zrange):
        output[i,:,:] = image3D[index,:,:]

    print("$ Reduced Z-planes from {} to {}".format(image3D.shape[0],output.shape[0]))

    return output


folder_raw = "/home/marcnol/data/pancreas_Exp6_2022_002_ROI"+os.sep

files_raw = glob.glob(folder_raw + os.sep + "*_ch00.tif")

data_raw = [imread(x) for x in files_raw]

for file,im in zip(files_raw,data_raw):
    output_TIFF = file.split(".npy")[0]+"_reduced_planes.tif"
    output_TIFF_2D = file.split(".npy")[0]+"_reduced_planes_2D.tif"

    new_im = reinterpolateZ(im,range(0, im.shape[0],2))
    new_im_2D = np.max(new_im,axis=0)

    imsave(output_TIFF,new_im)
    imsave(output_TIFF_2D, new_im_2D)
    print(f"Saving original NPY as TIFF:\n--> {output_TIFF}")
    print(f"N planes in new image: {new_im.shape}")

#%% represents segmentations and original RAW images
from astropy.visualization import SqrtStretch, simple_norm
from stardist import random_label_cmap
lbl_cmap = random_label_cmap()

percent=99.5

fig, axes = plt.subplots(3,1, figsize=(12,24), sharex=True, sharey=True)
axes = axes.ravel()


for im in data_raw:
    norm = simple_norm(im, "sqrt", percent=percent)
    axes[0].imshow(np.max(im,axis=0), cmap="Greys", origin="lower", norm=norm, alpha=1)
    axes[2].imshow(np.max(im,axis=0), cmap="Greys", origin="lower", norm=norm, alpha=1)

plot_flag = [True, False]


for _data2D,cmap, flag in zip(data2D,cmaps, plot_flag):
    if flag:
        axes[1].imshow(_data2D, cmap=lbl_cmap, alpha=1)
        axes[2].imshow(_data2D, cmap=lbl_cmap, alpha=0.4)
