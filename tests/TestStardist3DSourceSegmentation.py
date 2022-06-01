#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:31:50 2021

Script tested on LOPEVI. Before launching the segmentation, in a terminal :

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/cuda/lib64
    conda activate StarDIst3D

@author: fiche, marcnol
"""

from __future__ import print_function, unicode_literals, absolute_import, division

import numpy as np
import os
from datetime import datetime

from skimage import io, measure

from csbdeep.utils.tf import limit_gpu_memory

from stardist.models import StarDist3D

from imageProcessing.imageProcessing  import (
    image_adjust,
    reinterpolate_z,
    _segment_3d_volumes_stardist,
    _plot_image_3d,
    )

import matplotlib.pyplot as plt

# selects GPU
os.environ["CUDA_VISIBLE_DEVICES"]="1"

# Load the model
# --------------
if "atlantis" in os.uname()[1]:
    model_dir = '/home/marcnol/data/AImodels/StarDist/barcodes/'
    root_folder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
    file = root_folder+'scan_001_RT27_001_ROI_converted_decon_ch01.tif'
else:
    model_dir = '/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/'
    root_folder="/mnt/grey/DATA/users/marcnol/models/StarDist3D/training3Dbarcodes/dataset1/"
    file = root_folder+'scan_001_RT25_001_ROI_converted_decon_ch01_preProcessed_index0.tif'

model_name = 'stardist_18032021_single_loci'
model = StarDist3D(None, name=model_name, basedir=model_dir)
limit_gpu_memory(None, allow_growth=True)


#%% Segments image using Stardist

output = root_folder+"segmentedStardist3D_rebinned.npy"
im= io.imread(file).squeeze()

z_binning = 2

img = reinterpolate_z(im, range(0,im.shape[0],z_binning),mode='remove')

axis_norm = (0,1,2)
begin_time = datetime.now()

mask, labeled = _segment_3d_volumes_stardist(img,
                                area_min = 3,
                                area_max=1000,
                                nlevels=64,
                                contrast=0.001,
                                deblend_3d = True,
                                axis_norm=(0,1,2),
                                model_dir=model_dir,
                                model_name=model_name)

print("Elapsed time: {}".format(datetime.now() - begin_time))

np.save(output,mask)

#%% labels image and displays

labels = measure.label(mask)

img2distplay=img.copy()
img2distplay= image_adjust(img2distplay, lower_threshold=0.999, higher_threshold=0.99999)[0]

fig = _plot_image_3d(img2distplay,normalize=True, masks=labels)

