#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Replot 3-way co-localization matrix from saved .npy file
"""
import argparse
import os
import re
import select
import sys

import matplotlib.pyplot as plt
import numpy as np


def plot_threeway_matrix(
    pair_means,
    anchor_barcode,
    output_file="3way",
    distance_cutoff=0.2,
    vmin=None,
    vmax=None,
    cmap="RdBu",
):
    """
    Creates a heatmap of three-way co-localization frequencies using matplotlib.

    Parameters:
    ----------
    pair_means : dict
        Dictionary with (barcode1, barcode2) tuples as keys and mean frequencies as values
    anchor_barcode : int
        The anchor barcode number
    output_file : str
        Output file name for the plot
    """
    # Get all unique barcodes from the pairs
    all_barcodes = set()
    for b1, b2 in pair_means.keys():
        all_barcodes.add(b1)
        all_barcodes.add(b2)

    # Sort barcodes for consistent matrix indexing
    sorted_barcodes = sorted(all_barcodes)
    n_barcodes = len(sorted_barcodes)

    # Create empty matrices for means
    mean_matrix = np.zeros((n_barcodes, n_barcodes))

    # Create a mapping from barcode to matrix index
    barcode_to_idx = {b: i for i, b in enumerate(sorted_barcodes)}

    # Fill the matrices with the computed values
    for (b1, b2), mean_val in pair_means.items():
        i, j = barcode_to_idx[b1], barcode_to_idx[b2]
        mean_matrix[i, j] = mean_val
        mean_matrix[j, i] = mean_val  # Mirror the matrix (symmetric)

    # Create the figure and subplots for the mean frequencies
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the mean heatmap using matplotlib
    im = ax.imshow(
        mean_matrix,
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin if vmin is not None else 0,
        vmax=(
            vmax
            if vmax is not None
            else (0.9 * mean_matrix.max() if mean_matrix.max() > 0 else 1)
        ),
    )

    # Set up the axes with the correct labels
    ax.set_xticks(np.arange(len(sorted_barcodes)))
    ax.set_yticks(np.arange(len(sorted_barcodes)))
    ax.set_xticklabels(sorted_barcodes, fontsize=10)
    ax.set_yticklabels(sorted_barcodes, fontsize=10)

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)

    # Add perpendicular lines for the anchor barcode NOT WORKING
    if anchor_barcode in barcode_to_idx:
        anchor_idx = barcode_to_idx[anchor_barcode]

        # Horizontal line across the anchor barcode row
        ax.axhline(
            y=anchor_idx - 0.5, color="black", linestyle="-", linewidth=2, alpha=0.7
        )

        # Vertical line across the anchor barcode column
        ax.axvline(
            x=anchor_idx - 0.5, color="black", linestyle="-", linewidth=2, alpha=0.7
        )

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Co-localization frequency", fontsize=12)

    # Add title and labels
    ax.set_title(
        f"3-way co-localization with anchor {anchor_barcode}\n(distance cutoff: {distance_cutoff} µm)",
        fontsize=14,
    )
    ax.set_xlabel("Barcode #", fontsize=14)
    ax.set_ylabel("Barcode #", fontsize=14)

    # Adjust layout and saves npy matrix and image
    plt.tight_layout()
    output_filename = f"{output_file.split('.')[0]}_anchor_{anchor_barcode}_replot"

    np.save(f"{output_filename}.npy", mean_matrix)

    plt.savefig(f"{output_filename}.png", dpi=300)
    print(f"Saved three-way co-localization heatmap to: {output_filename}")
    plt.close()


def extract_anchor(filename):
    match = re.search(r"anchor_(\d+)", filename)
    return int(match.group(1)) if match else -1


def args_parser():
    parser = argparse.ArgumentParser(
        description="Replot 3-way co-localization matrices from .npy files"
    )
    parser.add_argument(
        "--input", nargs="*", help="List of input .npy matrix files.", default=[]
    )
    parser.add_argument("--pipe", action="store_true", help="Read filenames from stdin")
    parser.add_argument(
        "--vmin", type=float, default=None, help="Minimum color scale value"
    )
    parser.add_argument(
        "--vmax", type=float, default=None, help="Maximum color scale value"
    )
    parser.add_argument(
        "--cmap", default="RdBu", help="Matplotlib colormap (default: RdBu)"
    )
    return parser.parse_args()


def main():
    args = args_parser()

    matrix_files = list(args.input)

    if args.pipe and select.select([sys.stdin], [], [], 0.0)[0]:
        matrix_files += [
            line.strip() for line in sys.stdin if line.strip().endswith(".npy")
        ]

    if not matrix_files:
        print("No matrix files provided.")
        return

    for npy_file in matrix_files:
        if not os.path.exists(npy_file):
            print(f"File not found: {npy_file}")
            continue

        print(f"Loading matrix: {npy_file}")
        matrix = np.load(npy_file)
        anchor = extract_anchor(npy_file)

        if anchor == -1:
            print(f"Could not determine anchor from filename: {npy_file}")
            continue

        # Simulate a fake pair_means dictionary from the symmetric matrix
        n = matrix.shape[0]
        barcodes = list(range(n))
        pair_means = {}
        for i in range(n):
            for j in range(i + 1, n):
                pair_means[(barcodes[i], barcodes[j])] = matrix[i, j]

        output_base = os.path.splitext(npy_file)[0]

        plot_threeway_matrix(
            pair_means,
            anchor,
            output_file=output_base,
            distance_cutoff=0.2,
            vmin=args.vmin,
            vmax=args.vmax,
            cmap=args.cmap,
        )


if __name__ == "__main__":
    main()
