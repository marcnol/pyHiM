#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Localize sources (fluorescent spots) in 2D"""

from core.parameters import SegmentationParams
from imageProcessing.makeProjections import Feature


class Localize2D(Feature):
    def __init__(self, params: SegmentationParams):
        super().__init__(params)
        self.out_folder = self.params.localize_2d_folder
        self.name = "Localize2D"
