#!/usr/bin/env python3
"""
# Localization file copy and rename utility

At the moment the localization files form multiple ROIs have the same name. Thus we cannot copy them
onto the same folder to merge.

This script copies the files to a destination folder while renaming them to include
the unique folder name from their original path.

## Usage Examples

Basic usage with shell wildcard expansion:
```bash
./localization_cp_files.py --files /path/to/*/data/file.dat --destination_folder /output/path
```

For instance:
```bash
localization_cp_files.py --files /home/marcnol/grey/ProcessedData_2024/Experiment_76_David_DNAFISH_HiM_HRKit_G1E_20240711/*/localize_3d/data/localizations_3D_barcode.dat --destination_folder .
```

With explicit file list:
```bash
./localization_cp_files.py --files /path/to/001/data/file.dat /path/to/002/data/file.dat --destination_folder /output/path
```

## Dependencies
- Python 3.6+
- Standard library modules: os, shutil, re, argparse

## Notes
- The script automatically detects which part of the path varies (usually a folder number)
- Original file metadata (timestamps, permissions) is preserved
- Destination folder is created if it doesn't exist

## Author
marcnol
Created on: March 24, 2025
"""

import argparse
import os
import re
import shutil


def copy_and_rename_files(file_list, destination_folder):
    """
    Copy files from the list to the destination folder,
    renaming them to include the unique folder name.

    Args:
        file_list (list): List of file paths to copy
        destination_folder (str): Destination folder path
    """
    # Create destination folder if it doesn't exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        print(f"Created destination folder: {destination_folder}")

    if not file_list:
        print("No files provided to copy")
        return

    print(f"Processing {len(file_list)} files")

    # Find the common pattern and the variable part (folder name)
    # We'll analyze the first file to determine the pattern
    sample_paths = [file_list[0]]
    if len(file_list) > 1:
        sample_paths.append(file_list[1])

    # Extract the common path structure
    path_parts = [p.split("/") for p in sample_paths]

    # Find the varying segment (usually a number like 001, 002, etc.)
    varying_index = None
    for i in range(min(len(path_parts[0]), len(path_parts[-1]))):
        if len(sample_paths) == 1 or path_parts[0][i] != path_parts[-1][i]:
            varying_index = i
            break

    if varying_index is None and len(file_list) > 1:
        print("Could not identify the varying folder part in the file paths")
        varying_index = -3  # Default to 3rd from last segment if we can't determine
    elif varying_index is None:
        # For a single file, we'll assume the folder structure is /path/to/NUMBER/rest/of/path
        # Let's find a segment that looks like a number (e.g., 001, 002)
        for i, part in enumerate(path_parts[0]):
            if re.match(r"^\d+$", part):
                varying_index = i
                break

        if varying_index is None:
            varying_index = -3  # Default to 3rd from last segment if we can't determine

    # Process each file
    for file_path in file_list:
        # Extract the variable folder name
        path_parts = file_path.split("/")
        folder_name = (
            path_parts[varying_index] if varying_index < len(path_parts) else "unknown"
        )

        # Get the original filename
        original_filename = os.path.basename(file_path)

        # Create the new filename with the folder name appended
        new_filename = f"{original_filename.split('.')[0]}_{folder_name}.{original_filename.split('.')[-1]}"

        # Full path for the destination file
        dest_path = os.path.join(destination_folder, new_filename)

        # Copy the file
        shutil.copy2(file_path, dest_path)
        print(f"Copied: {file_path} -> {dest_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Copy and rename files based on folder names"
    )
    parser.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="List of files to copy (will be expanded by shell)",
    )
    parser.add_argument(
        "--destination_folder", required=True, help="Destination folder path"
    )

    args = parser.parse_args()

    copy_and_rename_files(args.files, args.destination_folder)
    print("Operation completed!")
