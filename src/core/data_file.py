#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data files module

Manage files operations, depending of DataManager.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from skimage import io

from core.pyhim_logging import print_log
from core.saving import image_show_with_values


class DataFile:
    def __init__(self, data=None):
        self.data = data

    def save(self, folder_path: str, basename: str):
        pass

    def delete_data(self):
        self.data = None


class TifFile:
    def __init__(self, path, basename, ext, label):
        self.all_path = path
        self.basename = basename
        self.extension = ext
        self.root = self.get_root()
        self.label = label

    def get_root(self):
        length = len(self.all_path) - len(self.basename) - 1 - len(self.extension)
        return self.all_path[:length]

    def load(self):
        return io.imread(self.all_path).squeeze()


class NpyFile(DataFile):
    def __init__(self, npy_data, dimension: int):
        super().__init__(npy_data)
        self.dim = dimension
        self.extension = "npy"
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def save(self, folder_path: str, basename: str):
        self.folder_path = folder_path + os.sep + "data"
        self.basename = f"{basename}_{str(self.dim)}d"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )
        if self.data.shape <= (1, 1):
            raise ValueError(f"Image is empty! Original file: {basename}.tif")
        np.save(self.path_name, self.data)
        print_log(f"$ Image saved to disk: {self.basename}.npy")


class Png2DFile(DataFile):
    def __init__(self, data):
        super().__init__(data)
        self.extension = "png"
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def save(self, folder_path: str, basename: str):
        self.folder_path = folder_path
        self.basename = f"{basename}_2d"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )
        fig = plt.figure()
        size = (10, 10)
        fig.set_size_inches(size)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()
        norm = ImageNormalize(stretch=SqrtStretch())

        ax.set_title("2D Data")

        fig.add_axes(ax)
        ax.imshow(self.data, origin="lower", cmap="Greys_r", norm=norm)
        fig.savefig(self.path_name)
        plt.close(fig)


class FocalPlaneMatrixFile(DataFile):
    def __init__(self, data, title):
        super().__init__(data)
        self.title = title
        self.extension = "png"
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def save(self, folder_path: str, basename: str):
        self.folder_path = folder_path
        self.basename = f"{basename}_focalPlaneMatrix"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )
        image_show_with_values(
            [self.data],
            title=self.title,
            output_name=self.path_name,
        )


class BlockAlignmentFile(DataFile):
    def __init__(self, relative_shifts, rms_image, contour):
        super().__init__()
        self.extension = "png"
        self.relative_shifts = relative_shifts
        self.rms_image = rms_image
        self.contour = contour
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def save(self, folder_path, basename):
        self.folder_path = folder_path
        self.basename = f"{basename}_block_alignments"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )

        # plotting
        fig, axes = plt.subplots(1, 2)
        ax = axes.ravel()
        fig.set_size_inches((10, 5))

        cbwindow = 3
        p_1 = ax[0].imshow(self.relative_shifts, cmap="terrain", vmin=0, vmax=cbwindow)
        ax[0].plot(self.contour.T[1], self.contour.T[0], linewidth=2, c="k")
        ax[0].set_title("abs(global-block) shifts, px")
        fig.colorbar(p_1, ax=ax[0], fraction=0.046, pad=0.04)

        p_2 = ax[1].imshow(
            self.rms_image,
            cmap="terrain",
            vmin=np.min(self.rms_image),
            vmax=np.max(self.rms_image),
        )
        ax[1].plot(self.contour.T[1], self.contour.T[0], linewidth=2, c="k")
        ax[1].set_title("RMS")
        fig.colorbar(p_2, ax=ax[1], fraction=0.046, pad=0.04)

        for axe in ax:
            axe.axis("off")

        fig.savefig(self.path_name)

        plt.close(fig)

    def delete_data(self):
        self.relative_shifts = None
        self.rms_image = None
        self.contour = None
