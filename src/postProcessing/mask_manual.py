#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 17:47:13 2023

@author: marcnol
"""
import sys
import argparse
import os
import numpy as np

from matplotlib import pyplot as plt
from skimage import exposure, io
from roipoly import MultiRoi

from astropy.visualization import simple_norm

def show_image(data_2d, normalization="simple", size=(10, 10)):
    fig = plt.figure()
    fig.set_size_inches(size)

    axes = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    axes.set_axis_off()

    if normalization == "simple":
        norm = simple_norm(data_2d, "sqrt", percent=99.9)

    fig.add_axes(axes)
    axes.imshow(data_2d, origin="lower", cmap="Greys_r", norm=norm)
    axes.set_title("2D Data")
    return fig, axes


def creates_user_mask(file_name,label):
    # loads image

    # displays image
    data = io.imread(file_name).squeeze()
    print(f'$ Image size: {data.shape}')
    data_2d = np.max(data,axis=0)
    
    fig = plt.figure()
    fig.set_size_inches((10,10))
    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    fig.add_axes(ax)
    ax.set_axis_off()    
    norm = simple_norm(data_2d, "sqrt", percent=99.9)
    ax.imshow(data_2d, origin="lower", cmap="Greys_r", norm=norm)

    print("Click on the button to add a new ROI")

    # Draw multiple rois
    multiroi_named = MultiRoi(roi_names=["ROI1", "ROI2"])

    number_rois = len(multiroi_named.rois)
    print(f"Number of rois drawn: {number_rois}")

    masks = np.zeros((data_2d.shape[0], data_2d.shape[1], number_rois))

    # Display image with all rois
    fig, _ = show_image(data_2d)

    roi_names = []
    i = 0
    for name, roi in multiroi_named.rois.items():
        roi.display_roi()
        roi.display_mean(data_2d)
        roi_names.append(name)
        masks[:, :, i] = roi.get_mask(data_2d)
        i += 1
    plt.legend(roi_names, bbox_to_anchor=(1.2, 1.05))
    # plt.show()

    print("Saving and closing image with rois...")
    fig.savefig(file_name.split('.')[0] + '_'+label+ ".png")
    plt.close(fig)

    # saves result
    output = os.path.basename(file_name).split('.')[0]+'_'+label+'.npy'
    np.save(output, masks)
    print(f"$ saved output image mask as {output}")
    
def main():

    print("="*10+"Started execution"+"="*10)

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input file name (TIF assumed)")
    parser.add_argument("--label", help="Label to add to image file name. Default=none")

    args = parser.parse_args()

    processing_list = {}

    if args.input:
        input_file = args.input
    else:
        print(
            ">> ERROR: you must provide a filename with the image file"
        )
        sys.exit(-1)

    if args.label:
        label = args.label
    else:
        label = ''
        
    print(f"parameters> input_file: {input_file}\n")

    creates_user_mask(input_file,label)
    
    print("="*9+"Finished execution"+"="*9)


if __name__ == "__main__":
    main()
    
