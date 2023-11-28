#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Segments mask in 3D"""

from core.parameters import SegmentationParams
from imageProcessing.makeProjections import Feature


class Mask3D(Feature):
    def __init__(self, params: SegmentationParams):
        super().__init__(params)
        self.out_folder = self.params.mask_3d_folder
        self.name = "Mask3D"
