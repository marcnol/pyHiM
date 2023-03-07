"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/plugins/guides.html?#widgets

Replace code below according to your needs.
"""

import numpy as np
from napari.types import ImageData
from napari.utils.notifications import show_info
from napari.layers import Image
from magicgui import magic_factory
from apifish.stack import projection


def do_focal_plane_projection(layer: ImageData) -> ImageData:
    focal_plane_matrix, _, block = projection.reinterpolate_focal_plane(layer)
    # reassembles image
    projected_image = projection.reassemble_images(focal_plane_matrix, block)
    return projected_image


@magic_factory(
    call_button="Run",
    radio_option={
        "widget_type": "RadioButtons",
        "orientation": "vertical",
        "choices": [
            ("Max", 1),
            ("Mean", 2),
            ("Median", 3),
            ("Sum", 4),
            ("Focus", 5),
            ("Focal", 6),
        ],
    },
)
def do_projection(layer: ImageData, radio_option=1) -> ImageData:
    if radio_option == 1:
        show_info("Max !")
        im_proj = projection.maximum_projection(layer)
        return (im_proj, {"shape": 0}, 'image')
    if radio_option == 2:
        show_info("Mean !")
        return projection.mean_projection(layer)
    if radio_option == 3:
        show_info("Median !")
        return projection.median_projection(layer)
    if radio_option == 4:
        show_info("Sum !")
        return projection.sum_projection(layer)
    if radio_option == 5:
        show_info("Focus !")
        return projection.focus_projection(layer)
    if radio_option == 6:
        show_info("Focal !")
        return do_focal_plane_projection(layer)
