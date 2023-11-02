#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Localize sources (fluorescent spots) in 3D"""

from core.parameters import SegmentationParams
from imageProcessing.makeProjections import Feature


class Localize3D(Feature):
    def __init__(self, params: SegmentationParams):
        super().__init__(params)
        self.out_folder = self.params.localize_3d_folder
        self.name = "Localize3D"
