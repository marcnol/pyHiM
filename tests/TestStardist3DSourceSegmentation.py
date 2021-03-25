#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:31:50 2021

Script tested on LOPEVI. Before launching the segmentation, in a terminal :
    
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/cuda/lib64
    conda activate StarDIst3D

@author: fiche
"""

from __future__ import print_function, unicode_literals, absolute_import, division

import numpy as np
import os
import time

from glob import glob
from tqdm import tqdm
from tifffile import imread
from csbdeep.utils import Path, normalize
from csbdeep.io import save_tiff_imagej_compatible
from csbdeep.data import PadAndCropResizer
from csbdeep.utils.tf import limit_gpu_memory

from stardist.models import StarDist3D

# Define the folders where the trained model and the test data are saved
# ----------------------------------------------------------------------

model_dir = '/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Annotated_data/data_loci_small/models/'
data_dir = '/mnt/PALM_dataserv/DATA/JB/2021/Data_single_loci/Data_test_small/'
model_name = 'stardist_18032021_single_loci'

# Load the model
# --------------

model = StarDist3D(None, name=model_name, basedir=model_dir)
limit_gpu_memory(None, allow_growth=True)

# Look for the data. Depending whether the images are large (>1000x1000pix) or 
# small, the procedure for the image reconstruction is not the same. 
# ------------------------------------------------------------------

os.chdir(data_dir)
X = sorted(glob('*raw.tif'))
axis_norm = (0,1,2)

for n_im in tqdm(range(len(X))):
    im_name = Path(X[n_im]).stem
    print(im_name)
    
    im = imread(X[n_im])
    im = normalize(im,1,99.8,axis=axis_norm)
    Lx = im.shape[1]

    t0= time.time()
    
    if Lx<1000:
        labels, details = model.predict_instances(im)
        
    else:
        resizer = PadAndCropResizer()
        axes = 'ZYX'
        
        im = resizer.before(im, axes, model._axes_div_by(axes))
        labels, polys = model.predict_instances(im, n_tiles=(1,8,8))
        labels = resizer.after(labels, axes)
        
    t1 = time.time() - t0
        
    print("Computation time : " + str(t1))

    mask = np.array(labels>0, dtype=int)
    
    label_name = im_name + '_label.tif'
    save_tiff_imagej_compatible(label_name, labels, axes='ZYX')
    
    mask_name = im_name + '_mask.tif'
    save_tiff_imagej_compatible(mask_name, mask, axes='ZYX')