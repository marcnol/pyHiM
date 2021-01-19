#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:39:42 2020

@author: marcnol
"""
import sys, os
import numpy as np
import matplotlib

matplotlib.rcParams["image.interpolation"] = None
import matplotlib.pyplot as plt

from tqdm import tqdm
from csbdeep.utils import Path, download_and_extract_zip_file
from tifffile import imread
from glob import glob
from csbdeep.utils import Path, download_and_extract_zip_file
from csbdeep.utils import normalize
from stardist import fill_label_holes, relabel_image_stardist, random_label_cmap
from stardist.matching import matching_dataset
from skimage import measure
from stardist import calculate_extents, gputools_available, _draw_polygons
from stardist.matching import matching
from stardist.models import Config2D, StarDist2D, StarDistData2D
import scipy.io as spio

np.random.seed(42)
lbl_cmap = random_label_cmap()

#%% functions


def plot_img_label(img, lbl, img_title="image", lbl_title="label", **kwargs):
    fig, (ai, al) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw=dict(width_ratios=(1.25, 1)))
    im = ai.imshow(img, cmap="gray", clim=(0, 1))
    ai.set_title(img_title)
    fig.colorbar(im, ax=ai)
    al.imshow(lbl, cmap=lbl_cmap)
    al.set_title(lbl_title)
    plt.tight_layout()


def random_fliprot(img, mask):
    assert img.ndim >= mask.ndim
    axes = tuple(range(mask.ndim))
    perm = tuple(np.random.permutation(axes))
    img = img.transpose(perm + tuple(range(mask.ndim, img.ndim)))
    mask = mask.transpose(perm)
    for ax in axes:
        if np.random.rand() > 0.5:
            img = np.flip(img, axis=ax)
            mask = np.flip(mask, axis=ax)
    return img, mask


def random_intensity_change(img):
    img = img * np.random.uniform(0.6, 2) + np.random.uniform(-0.2, 0.2)
    return img


def augmenter(x, y):
    """Augmentation of a single input/label image pair.
    x is an input image
    y is the corresponding ground-truth label image
    """
    x, y = random_fliprot(x, y)
    x = random_intensity_change(x)
    # add some gaussian noise
    sig = 0.02 * np.random.uniform(0, 1)
    x = x + sig * np.random.normal(0, 1, x.shape)
    return x, y


def loadsTrainingDataJB(rootFolder):

    folderMasks = "Labeled_images"
    folderImages = "Original_images"

    ListMasks, ListImages = [], []
    ListMasks = sorted(glob(rootFolder + os.sep + folderMasks + os.sep + "*.tif"))

    baseNameMasks = [os.path.basename(basename) for basename in ListMasks]

    for target in baseNameMasks:

        expectedFolder = rootFolder + os.sep + folderImages + os.sep + target.split(".tif")[0]
        expectedFolder45 = rootFolder + os.sep + folderImages + os.sep + target.split("_45.tif")[0]

        if os.path.exists(expectedFolder):
            fileName = expectedFolder + os.sep + "00_Raw_Embryo_segmentation.mat"
            if os.path.exists(fileName):
                ListImages.append(fileName)

        elif os.path.exists(expectedFolder45):
            fileName = expectedFolder45 + os.sep + "00_Raw_Embryo_segmentation_45.mat"
            if os.path.exists(fileName):
                ListImages.append(fileName)

    print("Number of Masks: {}".format(len(ListMasks)))
    print("Number of Images: {}".format(len(ListImages)))

    if len(ListMasks) == len(ListImages):

        Y = list(map(imread, ListMasks))
        # measure.label
        Y = [measure.label(y) for y in Y]
        Xmat = list(map(spio.loadmat, ListImages))
        # for y in Ymat:
        #     if 'im_raw' in y.keys():
        #         Y.append(y['im_raw'])
        #     elif 'im_raw_45' in y.keys():
        #         Y.append(y['im_raw_45'])

        X = [x[list(x.keys())[-1]] for x in Xmat]

    else:
        print("Warning, something is wrong...")
    # mat = spio.loadmat(rootFolder+fileName, squeeze_me=True)
    return X, Y


def example(model, i, show_dist=True):
    img = normalize(X[i], 1, 99.8, axis=axis_norm)
    labels, details = model.predict_instances(img)

    plt.figure(figsize=(13, 10))
    img_show = img if img.ndim == 2 else img[..., 0]
    coord, points, prob = details["coord"], details["points"], details["prob"]
    plt.subplot(121)
    plt.imshow(img_show, cmap="gray")
    plt.axis("off")
    a = plt.axis()
    _draw_polygons(coord, points, prob, show_dist=show_dist)
    plt.axis(a)
    plt.subplot(122)
    plt.imshow(img_show, cmap="gray")
    plt.axis("off")
    plt.imshow(labels, cmap=lbl_cmap, alpha=0.5)
    plt.tight_layout()
    plt.show()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    #%% sets parameters
    ownTrainingSet = True
    sampleName = "nc14"
    sep = "_"
    # 32 is a good default choice (see 1_data.ipynb)
    n_rays = 32
    quick_demo = False
    train_epochs = 20
    # Predict on subsampled grid for increased efficiency and larger field of view
    grid = (2, 2)
    modelName = "".join(
        [
            "stardist",
            sep,
            sampleName,
            sep,
            "nrays:",
            str(n_rays),
            sep,
            "epochs:",
            str(train_epochs),
            sep,
            "grid:",
            str(grid[0]),
        ]
    )

    #%% loads data

    if ownTrainingSet:
        rootFolder = "/mnt/PALM_dataserv/DATA/JB/JB/Sara/Deep_Learning/Training_data/Embryo/Marcelo_embryo_data/DAPI_nuclei_segmentation/stage_14"
        X, Y = loadsTrainingDataJB(rootFolder)
        print("Loading data from: {}".format(rootFolder))

    else:

        download_and_extract_zip_file(
            url="https://github.com/mpicbg-csbd/stardist/releases/download/0.1.0/dsb2018.zip",
            targetdir="data",
            verbose=1,
        )

        X = sorted(glob("data/dsb2018/train/images/*.tif"))
        Y = sorted(glob("data/dsb2018/train/masks/*.tif"))
        assert all(Path(x).name == Path(y).name for x, y in zip(X, Y))

        X, Y = X[:10], Y[:10]

        X = list(map(imread, X))
        Y = list(map(imread, Y))

        print("Loading Example data from StarDist")

    #%% preprocesses datasets
    i = min(4, len(X) - 1)
    img, lbl = X[i], fill_label_holes(Y[i])
    assert img.ndim in (2, 3)
    img = img if img.ndim == 2 else img[..., :3]
    # assumed axes ordering of img and lbl is: YX(C)

    plt.figure(figsize=(16, 10))
    plt.subplot(121)
    plt.imshow(img, cmap="gray")
    plt.axis("off")
    plt.title("Raw image")
    plt.subplot(122)
    plt.imshow(lbl, cmap=lbl_cmap)
    plt.axis("off")
    plt.title("GT labels")
    None

    n_channel = 1 if X[0].ndim == 2 else X[0].shape[-1]

    #%% Normalize images and fill small label holes.

    axis_norm = (0, 1)  # normalize channels independently
    # axis_norm = (0,1,2) # normalize channels jointly
    if n_channel > 1:
        print(
            "Normalizing image channels %s." % ("jointly" if axis_norm is None or 2 in axis_norm else "independently")
        )
        sys.stdout.flush()

    X = [normalize(x, 1, 99.8, axis=axis_norm) for x in tqdm(X)]
    Y = [fill_label_holes(y) for y in tqdm(Y)]

    #%% splits into training and validation datasets

    assert len(X) > 1, "not enough training data"
    rng = np.random.RandomState(42)
    ind = rng.permutation(len(X))
    n_val = max(1, int(round(0.15 * len(ind))))
    ind_train, ind_val = ind[:-n_val], ind[-n_val:]
    X_val, Y_val = [X[i] for i in ind_val], [Y[i] for i in ind_val]
    X_trn, Y_trn = [X[i] for i in ind_train], [Y[i] for i in ind_train]
    print("number of images: %3d" % len(X))
    print("- training:       %3d" % len(X_trn))
    print("- validation:     %3d" % len(X_val))

    #%%  shows example

    i = min(9, len(X) - 1)
    # i=3
    img, lbl = X[i], Y[i]
    assert img.ndim in (2, 3)
    img = img if (img.ndim == 2 or img.shape[-1] == 3) else img[..., 0]
    plot_img_label(img, lbl)
    None

    #%%  Show parameters for the network
    print(Config2D.__doc__)

    #%% Setting paramters doe TRAINING network and DATA augmentation

    # Use OpenCL-based computations for data generator during training (requires 'gputools')
    use_gpu = False and gputools_available()

    conf = Config2D(n_rays=n_rays, grid=grid, use_gpu=use_gpu, n_channel_in=n_channel, train_epochs=train_epochs,)
    # print(conf)
    vars(conf)

    if use_gpu:
        from csbdeep.utils.tf import limit_gpu_memory

        # adjust as necessary: limit GPU memory to be used by TensorFlow to leave some to OpenCL-based computations
        limit_gpu_memory(0.8)

    # THIS SETS up the model to be trained, including name.
    # Using default values: prob_thresh=0.5, nms_thresh=0.4.
    model = StarDist2D(conf, name=modelName, basedir="models")

    # Check if the neural network has a large enough field of view to see up to the boundary of most objects.
    median_size = calculate_extents(list(Y), np.median)
    fov = np.array(model._axes_tile_overlap("YX"))
    print(f"median object size:      {median_size}")
    print(f"network field of view :  {fov}")

    if any(median_size > fov):
        print("WARNING: median object size larger than field of view of the neural network.")

    # plot some augmented examples
    img, lbl = X[0], Y[0]
    plot_img_label(img, lbl)
    for _ in range(3):
        img_aug, lbl_aug = augmenter(img, lbl)
        plot_img_label(img_aug, lbl_aug, img_title="image augmented", lbl_title="label augmented")

    #%% TRAINING

    if quick_demo:
        print(
            "NOTE: This is only for a quick demonstration!\n"
            "      Please set the variable 'quick_demo = False' for proper (long) training.",
            file=sys.stderr,
            flush=True,
        )
        model.train(X_trn, Y_trn, validation_data=(X_val, Y_val), augmenter=augmenter, epochs=2, steps_per_epoch=10)

        print(
            "====> Stopping training and loading previously trained demo model from disk.", file=sys.stderr, flush=True
        )
        model = StarDist2D.from_pretrained("2D_demo")
    else:
        model.train(X_trn, Y_trn, validation_data=(X_val, Y_val), augmenter=augmenter)
    None

    #%% Optimize thresholds
    # While the default values for the probability and non-maximum suppression thresholds already yield good results in
    # many cases, we still recommend to adapt the thresholds to your data. The optimized threshold values are saved to
    # disk and will be automatically loaded with the model.

    if quick_demo:
        model.optimize_thresholds(X_val[:2], Y_val[:2])
    else:
        model.optimize_thresholds(X_val, Y_val)

    #%% Evaluation and Detection Performance

    # First predict the labels for all validation images:

    for index in range(3):

        # example(model, index)
        Y_val_pred = [
            model.predict_instances(x, n_tiles=model._guess_n_tiles(x), show_tile_progress=False)[0]
            for x in tqdm(X_val)
        ]

        plot_img_label(X_val[index], Y_val[index], lbl_title="label GT")
        plot_img_label(X_val[index], Y_val_pred[index], lbl_title="label Pred")

    #%% Choose several IoU thresholds $\tau$ that might be of interest and for each compute matching statistics for the validation data.

    taus = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    stats = [matching_dataset(Y_val, Y_val_pred, thresh=t, show_progress=False) for t in tqdm(taus)]
    print("\nMatching stats for 0.5: {}\n".format(stats[taus.index(0.5)]))

    # Plot the matching statistics and the number of true/false positives/negatives as a function of the IoU threshold $\tau$.

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

    for m in ("precision", "recall", "accuracy", "f1", "mean_true_score", "mean_matched_score", "panoptic_quality"):
        ax1.plot(taus, [s._asdict()[m] for s in stats], ".-", lw=2, label=m)
    ax1.set_xlabel(r"IoU threshold $\tau$")
    ax1.set_ylabel("Metric value")
    ax1.grid()
    ax1.legend()

    for m in ("fp", "tp", "fn"):
        ax2.plot(taus, [s._asdict()[m] for s in stats], ".-", lw=2, label=m)
    ax2.set_xlabel(r"IoU threshold $\tau$")
    ax2.set_ylabel("Number #")
    ax2.grid()
    ax2.legend()
