#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
class used to define the order of execution inside the pipeline.
"""

from core.parameters import MatrixParams
from imageProcessing.makeProjections import Feature


class BuildMatrixTempo(Feature):
    def __init__(self, params: MatrixParams):
        super().__init__(params)
        self.out_folder = self.params.folder
        self.name = "BuildMatrix"
