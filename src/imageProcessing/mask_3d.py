#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Segments mask in 3D"""

from core.parameters import SegmentationParams
from imageProcessing.makeProjections import Feature
# from imageProcessing.alignImages3D import Class3D
# import apply DF function 

class Mask3D(Feature):
    def __init__(self, params: SegmentationParams):
        super().__init__(params)
        self.out_folder = self.params.mask_3d_folder
        self.name = "Mask3D"
        #  if DF exists :
            #    warped_new_np = applyDF(new_image_np,displacement_fields_np, fixed_image_np.shape)
            #    warped_new_sitk = to_sitk(warped_new_np, ref_img=fixed_image)
            #    write_image(warped_new_sitk, "name_image_corrected.tif")
