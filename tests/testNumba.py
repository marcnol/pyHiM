#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:38:15 2020

@author: marcnol
"""

# import numba
# import numpy as np

# from datetime import datetime

#%% mandelbrot


from timeit import default_timer as timer

try:
    from matplotlib.pylab import imshow, show

    have_mpl = True
except ImportError:
    have_mpl = False
import numpy as np
from numba import jit


@jit(nopython=True)
# @jit(nopython=True, parallel=True)
def mandel(x, y, max_iters):
    """
    Given the real and imaginary parts of a complex number,
    determine if it is a candidate for membership in the Mandelbrot
    set given a fixed number of iterations.
    """
    i = 0
    c = complex(x, y)
    z = 0.0j
    for i in range(max_iters):
        z = z * z + c
        if (z.real * z.real + z.imag * z.imag) >= 4:
            return i

    return 255


@jit(nopython=True)
# @jit(nopython=True, parallel=True)
def create_fractal(min_x, max_x, min_y, max_y, image, iters):
    height = image.shape[0]
    width = image.shape[1]

    pixel_size_x = (max_x - min_x) / width
    pixel_size_y = (max_y - min_y) / height
    for x in range(width):
        real = min_x + x * pixel_size_x
        for y in range(height):
            imag = min_y + y * pixel_size_y
            color = mandel(real, imag, iters)
            image[y, x] = color

    return image


image = np.zeros((500 * 2, 750 * 2), dtype=np.uint8)
s = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20)
e = timer()

s = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20)
e = timer()

print("Time elapsed:{}".format(e - s))
if have_mpl:
    imshow(image)
    show()
