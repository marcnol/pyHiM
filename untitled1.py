

from __future__ import absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal

from skimage.filters import threshold_local, gaussian
from skimage.util.apply_parallel import apply_parallel
from datetime import datetime
from skimage import io
from skimage import exposure
from imageProcessing.imageProcessing  import (
    _removesInhomogeneousBackground2D,
    _removesInhomogeneousBackground,
    imageAdjust,
    _segments3DvolumesByThresholding,
    savesImageAsBlocks,
    display3D,
    combinesBlocksImageByReprojection,
    display3D_assembled,
    appliesXYshift3Dimages,
    )
# a = np.arange(144).reshape(12, 12).astype(float)
rootFolder="/home/marcnol/data/Embryo_debug_dataset/test_dataset/"
file = rootFolder+'scan_001_RT31_001_ROI_converted_decon_ch01.tif'
imageRaw = io.imread(file).squeeze()


#%%
begin_time = datetime.now()

image3D = exposure.rescale_intensity(imageRaw, out_range=(0, 1))
image3D= _removesInhomogeneousBackground(image3D)

print("Elapsed time: {}".format(datetime.now() - begin_time))

begin_time = datetime.now()

image3D = apply_parallel(exposure.rescale_intensity, imageRaw, #chunks=(60, 256, 256),
#                         extra_arguments=(),
                         extra_keywords={'out_range': (0, 1)})

image3D = apply_parallel(_removesInhomogeneousBackground, image3D, #chunks=(60, 256, 256),
#                         extra_arguments=(),
                         extra_keywords={'out_range': (0, 1)})

print("Elapsed time: {}".format(datetime.now() - begin_time))


