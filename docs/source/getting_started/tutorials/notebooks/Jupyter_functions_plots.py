#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:55:21 2022

@author: Olivier Messina
"""

import glob
import os
import shutil

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from IPython.core.display import HTML
from IPython.display import Image
from matplotlib import rcParams
from tifffile import TiffWriter


def copy_localization_table(Input_folder):
    source = os.path.join(Input_folder, "localizations_3D_barcode_example.dat")
    dest = os.path.join(
        Input_folder,
        "analysis_localizations",
        "localize_3d",
        "data",
        "localizations_3D_barcode.dat",
    )
    folder_path = os.path.join(
        Input_folder, "analysis_localizations", "localize_3d", "data"
    )
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    shutil.copyfile(source, dest)

    # copy the parameters.json in new root
    source = os.path.join(Input_folder, "parameters.json")
    dest = os.path.join(Input_folder, "analysis_localizations", "parameters.json")
    shutil.copyfile(source, dest)


def plot_zprojection(Input_folder, RTs_references, titles, datatype):
    # Figure size in inches optional
    rcParams["figure.figsize"] = 15, 15

    # Enter in the folder contaning the zProjected output images
    os.chdir(Input_folder + "project")

    if datatype == "DAPI":
        # Create list of images
        img_A = glob.glob("*DAPI*" + "*ch01*" + ".png")
        img_B = glob.glob("*DAPI*" + "*ch00*" + ".png")
        Concat_images = [img_A, img_B]

    if datatype == "RT":
        # Create list of images
        img_A = glob.glob("*" + RTs_references + "*" + "*ch00*" + ".png")
        img_B = glob.glob("*" + RTs_references + "*" + "*ch01*" + ".png")
        Concat_images = [img_A, img_B]

    # plots images
    fig, ax = plt.subplots(1, 2)
    for x, file in enumerate(Concat_images):
        # Read images
        Image = mpimg.imread(Input_folder + "project" + os.sep + file[0])[
            0:1000, 0:1000
        ]
        # Display images
        ax[x].set_title(titles[x])
        ax[x].imshow(Image)


def plot_alignment(Input_folder, RTs_references, titles):
    # Figure size in inches optional
    rcParams["figure.figsize"] = 15, 10

    # Create list of images
    jpgFilenamesList_RTs_alignement_Difference = glob.glob(
        Input_folder
        + "register_global"
        + os.sep
        + "*"
        + RTs_references
        + "*"
        + "_referenceDifference.png"
    )
    jpgFilenamesList_RTs_alignement_overlay = glob.glob(
        Input_folder
        + "register_global"
        + os.sep
        + "*"
        + RTs_references
        + "*"
        + "_overlay"
        + "*"
        + ".png"
    )

    img_A = mpimg.imread(jpgFilenamesList_RTs_alignement_Difference[0])[
        1000:2000, 1000:2000
    ]  # Zoom in the region of interest
    img_B = mpimg.imread(jpgFilenamesList_RTs_alignement_Difference[0])[
        1000:2000, 3500:4500
    ]  # Zoom in the region of interest
    img_C = (
        mpimg.imread(jpgFilenamesList_RTs_alignement_overlay[0])[1000:2000, 700:1700]
        * 5
    )  # Zoom in the region of interest

    # Concatenates images
    Concat_images = [img_A, img_B, img_C]

    # Plot
    fig, ax = plt.subplots(1, 3)
    for x, file in enumerate(titles):
        # Read images
        Image = Concat_images[x]
        # Display images
        ax[x].set_title(titles[x])
        ax[x].imshow(Image)


def show_plot(titles, imgs, files):
    if len(imgs) > 1:
        fig, ax = plt.subplots(1, 2)
        for title, img, axis in zip(titles, imgs, ax):
            axis.set_title(title)
            axis.imshow(img)
    else:
        plt.imshow(imgs[0])
        plt.axis("off")


def show_plot_traces(titles, imgs, scales):
    fig, ax = plt.subplots(1, 2)
    for title, img, axis, sca in zip(titles, imgs, ax, scales):
        axis.set_title(title + " (μm)")
        axis.imshow(img, extent=(0, sca, 0, sca))


def plot_segment_object(Input_folder, titles, datatype, RTs_references=""):
    # Figure size in inches optional
    rcParams["figure.figsize"] = 15, 10

    if datatype == "DAPI" or datatype == "RT":
        if datatype == "RT":
            # Create list of images
            files = glob.glob(
                Input_folder
                + "localize_3d"
                + os.sep
                + "*"
                + RTs_references
                + "*"
                + "_3DimageNlocalizations.png"
            )
            imgs = [mpimg.imread(files[0])[500:4500, 500:4500]]

        if datatype == "DAPI":
            # Create list of images
            files = glob.glob(
                Input_folder + "mask_3d" + os.sep + "*" + "DAPI" + "*" + "_3Dmasks.png"
            )
            imgs = [mpimg.imread(files[0])[1500:3600, 500:4500]]

        print(f"$ Will plot: {files}")

        # Plots images
        show_plot(titles, imgs, files)

    elif datatype == "TRACES":
        # Create list of images
        tracing_folder = os.path.join(Input_folder, "tracing")
        files = glob.glob(
            tracing_folder + os.sep + "*" + "KDtree" + "*" + "_XYZ" + "*" + ".png"
        )
        imgs = [mpimg.imread(x) for x in files]
        # crops png to get full and zoomed images
        # this img represent 200*200 micron-meters
        imgs = [x[70:1840, 170:1940] for x in imgs]
        scales = [200]
        imgs.append(imgs[0][1470:1770, 1390:1690])
        scales.append(200 * (1770 - 1470) / (1840 - 70))

        print(f"$ Will plot: {files}")

        # Plots images
        show_plot_traces(titles, imgs, scales)


def plot_matrix(input_folder, data_type="proximity"):
    # figure size in inches optional
    rcParams["figure.figsize"] = 15, 15

    tracing_folder = os.path.join(input_folder, "tracing")
    # Create list of images
    if data_type == "PWD":
        files = glob.glob(tracing_folder + os.sep + "*_PWDmatrixKDE.png")
    elif data_type == "proximity":
        files = glob.glob(tracing_folder + os.sep + "*HiMmatrix.png")
    elif data_type.split(",")[0] == "3D_alignments":
        files = glob.glob(
            input_folder
            + "register_local"
            + os.sep
            + "*"
            + data_type.split(",")[1]
            + "*3Dalignments.png"
        )

    titles = [data_type]

    if len(files) > 0:
        img = mpimg.imread(files[0])
        show_plot(titles, [img], "-")
