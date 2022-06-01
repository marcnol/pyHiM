#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_almost_equal

from skimage.filters import threshold_local, gaussian
from skimage.util.apply_parallel import apply_parallel
from datetime import datetime
from skimage import io
from skimage import exposure
from imageProcessing.imageProcessing  import (
    _remove_inhomogeneous_background_2d,
    _remove_inhomogeneous_background,
    image_adjust,
    _segment_3d_volumes_by_thresholding,
    save_image_as_blocks,
    display_3d,
    combine_blocks_image_by_reprojection,
    display_3d_assembled,
    apply_xy_shift_3d_images,
    )
# a = np.arange(144).reshape(12, 12).astype(float)
root_folder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
file = root_folder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
imageRaw = io.imread(file).squeeze()


#%%
# begin_time = datetime.now()

# image_3d = exposure.rescale_intensity(imageRaw, out_range=(0, 1))
# image_3d= _remove_inhomogeneous_background(image_3d)

# print("Elapsed time: {}".format(datetime.now() - begin_time))

# #%%
# begin_time = datetime.now()

# image_3d = apply_parallel(exposure.rescale_intensity, imageRaw, #chunks=(60, 256, 256),
# #                         extra_arguments=(),
#                          extra_keywords={'out_range': (0, 1)})

# image_3d = apply_parallel(_remove_inhomogeneous_background, image_3d, chunks=(1, 2048, 2048),
# #                         extra_arguments=(),
#                           extra_keywords={'verbose': False})

# print("Elapsed time: {}".format(datetime.now() - begin_time))

#%%
from dask.distributed import Client, LocalCluster

parallel=True
image3D_small=imageRaw[0:4,:,:]
if parallel:
    n_threads = 10

    print("Go to http://localhost:8787/status for information on progress...")

    # cluster = LocalCluster(n_workers=n_threads)

    # client = Client(cluster)

    image_3d_aligned = apply_xy_shift_3d_images(image3D_small, np.array([1,1]))

    # client.close()
    # cluster.close()

else:
    image_3d_aligned = apply_xy_shift_3d_images(image3D_small, np.array([1,1]))