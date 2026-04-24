#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.ndimage import shift as shift_image

from imageProcessing.alignImages import apply_shift_3d_images


def test_apply_shift_3d_images_from_xy_shift():
    image = np.zeros((5, 7, 9), dtype=float)
    image[2, 3, 4] = 1.0

    shifted = apply_shift_3d_images(image, [1, -2])
    expected = shift_image(image, [0, 1, -2])

    assert np.allclose(shifted, expected)


def test_apply_shift_3d_images_from_zxy_shift():
    image = np.zeros((5, 7, 9), dtype=float)
    image[2, 3, 4] = 1.0

    shifted = apply_shift_3d_images(image, [1, -2, 3])
    expected = shift_image(image, [1, -2, 3])

    assert np.allclose(shifted, expected)
