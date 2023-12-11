#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrects drift in 3D"""

from core.parameters import RegistrationParams
from imageProcessing.makeProjections import Feature


class RegisterLocal(Feature):
    def __init__(self, params: RegistrationParams):
        super().__init__(params)
        self.out_folder = self.params.register_local_folder
        self.name = "RegisterLocal"
