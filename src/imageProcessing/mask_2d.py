#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Segments mask in 2D"""

from core.parameters import SegmentationParams
from imageProcessing.makeProjections import Feature


class Mask2D(Feature):
    def __init__(self, params: SegmentationParams):
        super().__init__(params)
        self.out_folder = self.params.mask_2d_folder
        self.name = "Mask2D"
