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
from imageProcessing.imageProcessing import image_adjust


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
    def __init__(self, npy_data, status: str):
        super().__init__(npy_data)
        self.status = status
        self.extension = "npy"
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def save(self, folder_path: str, basename: str):
        self.folder_path = folder_path + os.sep + "data"
        self.basename = basename + self.status
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


class BothImgRbgFile(DataFile):
    def __init__(self, image1_uncorrected, image2_corrected_raw):
        super().__init__()
        self.extension = "png"
        self.image1_uncorrected = image1_uncorrected
        self.image2_corrected_raw = image2_corrected_raw
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def delete_data(self):
        self.image1_uncorrected = None
        self.image2_corrected_raw = None

    def save(self, folder_path, basename):
        """
        Overlays two images as R and B and saves them to output file
        """
        self.folder_path = folder_path
        self.basename = f"{basename}_overlay_corrected"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )
        sz = self.image1_uncorrected.shape
        img_1, img_2 = (
            self.image1_uncorrected / self.image1_uncorrected.max(),
            self.image2_corrected_raw / self.image2_corrected_raw.max(),
        )
        img_1, _, _, _, _ = image_adjust(
            img_1, lower_threshold=0.5, higher_threshold=0.9999
        )
        img_2, _, _, _, _ = image_adjust(
            img_2, lower_threshold=0.5, higher_threshold=0.9999
        )
        fig, ax1 = plt.subplots()
        fig.set_size_inches((30, 30))
        null_image = np.zeros(sz)
        rgb = np.dstack([img_1, img_2, null_image])
        ax1.imshow(rgb)
        ax1.axis("off")
        fig.savefig(self.path_name)
        plt.close(fig)


class RefDiffFile(DataFile):
    def __init__(self, preprocessed_ref, shifted_img, preprocessed_img):
        super().__init__()
        self.extension = "png"
        self.preprocessed_ref = preprocessed_ref
        self.shifted_img = shifted_img
        self.preprocessed_img = preprocessed_img
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def delete_data(self):
        self.preprocessed_ref = None
        self.shifted_img = None
        self.preprocessed_img = None

    def save(self, folder_path, basename):
        """
        Overlays two images as R and B and saves them to output file
        """
        self.folder_path = folder_path
        self.basename = f"{basename}_referenceDifference"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )

        img_1, img_2 = (
            self.preprocessed_ref / self.preprocessed_ref.max(),
            self.preprocessed_img / self.preprocessed_img.max(),
        )
        img_3, img_4 = (
            self.preprocessed_ref / self.preprocessed_ref.max(),
            self.shifted_img / self.shifted_img.max(),
        )

        img_1, _, _, _, _ = image_adjust(
            img_1, lower_threshold=0.5, higher_threshold=0.9999
        )
        img_2, _, _, _, _ = image_adjust(
            img_2, lower_threshold=0.5, higher_threshold=0.9999
        )
        img_3, _, _, _, _ = image_adjust(
            img_3, lower_threshold=0.5, higher_threshold=0.9999
        )
        img_4, _, _, _, _ = image_adjust(
            img_4, lower_threshold=0.5, higher_threshold=0.9999
        )

        cmap = "seismic"

        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches((60, 30))

        ax1.imshow(img_1 - img_2, cmap=cmap)
        ax1.axis("off")
        ax1.set_title("uncorrected")

        ax2.imshow(img_3 - img_4, cmap=cmap)
        ax2.axis("off")
        ax2.set_title("corrected")

        fig.savefig(self.path_name)

        plt.close(fig)


class EqualizationHistogramsFile(DataFile):
    def __init__(self, i_histogram, lower_threshold):
        super().__init__()
        self.extension = "png"
        self.i_histogram = i_histogram
        self.lower_threshold = lower_threshold
        self.folder_path = ""
        self.basename = ""
        self.path_name = ""

    def delete_data(self):
        self.i_histogram = None
        self.lower_threshold = None

    def save(self, folder_path, basename):
        """
        Overlays two images as R and B and saves them to output file
        """
        self.folder_path = folder_path
        self.basename = f"{basename}_intensityHist"
        self.path_name = (
            self.folder_path + os.sep + self.basename + "." + self.extension
        )

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

        ax1.plot(self.i_histogram["Im1"][0][1], self.i_histogram["Im1"][0][0])
        ax2.plot(self.i_histogram["Im2"][0][1], self.i_histogram["Im2"][0][0])
        ax3.plot(self.i_histogram["Im1"][1][1], self.i_histogram["Im1"][1][0])
        ax4.plot(self.i_histogram["Im2"][1][1], self.i_histogram["Im2"][1][0])
        ax3.set_yscale("log")
        ax4.set_yscale("log")
        ax1.vlines(
            self.lower_threshold["Im1"],
            0,
            self.i_histogram["Im1"][0][0].max(),
            colors="r",
        )
        ax2.vlines(
            self.lower_threshold["Im2"],
            0,
            self.i_histogram["Im2"][0][0].max(),
            colors="r",
        )
        plt.savefig(self.path_name)
        plt.close(fig)
