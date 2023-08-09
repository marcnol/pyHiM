#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for projection
"""

import numpy as np
import scipy.optimize as spo
from apifish.stack import projection

from core.pyhim_logging import print_log


# @jit(nopython=True)
def calculate_zrange(idata, parameters):
    """
    Calculates the focal planes based max standard deviation
    Finds best focal plane by determining the max of the std deviation vs z curve
    """
    num_planes = (
        parameters.param_dict["zProject"]["zmax"]
        - parameters.param_dict["zProject"]["zmin"]
    )
    std_matrix = np.zeros(num_planes)
    mean_matrix = np.zeros(num_planes)

    # calculate STD in each plane
    for i in range(num_planes):
        std_matrix[i] = np.std(idata[i])
        mean_matrix[i] = np.mean(idata[i])

    max_std = np.max(std_matrix)
    i_focus_plane = np.where(std_matrix == max_std)[0][0]
    # Select a window to avoid being on the edges of the stack

    if i_focus_plane < parameters.param_dict["zProject"]["windowSecurity"] or (
        i_focus_plane > num_planes - parameters.param_dict["zProject"]["windowSecurity"]
    ):
        focus_plane = i_focus_plane
    else:
        # interpolate zfocus
        axis_z = range(
            max(
                parameters.param_dict["zProject"]["zmin"],
                i_focus_plane - parameters.param_dict["zProject"]["windowSecurity"],
                min(
                    parameters.param_dict["zProject"]["zmax"],
                    i_focus_plane + parameters.param_dict["zProject"]["windowSecurity"],
                ),
            )
        )

        std_matrix -= np.min(std_matrix)
        std_matrix /= np.max(std_matrix)

        try:
            fitgauss = spo.curve_fit(
                projection.gaussian, axis_z, std_matrix[axis_z[0] : axis_z[-1] + 1]
            )
            focus_plane = int(fitgauss[0][1])
        except RuntimeError:
            print_log("Warning, too many iterations")
            focus_plane = i_focus_plane

    zmin = max(
        parameters.param_dict["zProject"]["windowSecurity"],
        focus_plane - parameters.param_dict["zProject"]["zwindows"],
    )
    zmax = min(
        num_planes,
        parameters.param_dict["zProject"]["windowSecurity"] + num_planes,
        focus_plane + parameters.param_dict["zProject"]["zwindows"],
    )
    zrange = range(zmin, zmax + 1)

    return focus_plane, zrange


def project_image_2d(img, z_range, mode):
    # sums images
    i_collapsed = None

    if mode == "MIP":
        # Max projection of selected planes
        i_collapsed = projection.maximum_projection(
            img[z_range[1][0] : z_range[1][-1] + 1]
        )
    elif mode == "sum":
        # Sums selected planes
        i_collapsed = projection.sum_projection(
            img[z_range[1][0] : (z_range[1][-1] + 1)]
        )
    else:
        print_log(f"ERROR: mode not recognized. Expected: MIP or sum. Read: {mode}")

    return i_collapsed


def reinterpolate_focal_plane(data, param_dict):
    if "blockSize" in param_dict["zProject"]:
        block_size_xy = param_dict["zProject"]["blockSize"]
    else:
        block_size_xy = 128

    if "zwindows" in param_dict["zProject"]:
        window = param_dict["zProject"]["zwindows"]
    else:
        window = 0

    focal_plane_matrix, z_range, block = projection.reinterpolate_focal_plane(
        data, block_size_xy=block_size_xy, window=window
    )

    # reassembles image
    output = projection.reassemble_images(focal_plane_matrix, block, window=window)

    return output, focal_plane_matrix, z_range
