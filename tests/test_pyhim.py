import sys
import os
import tempfile
import shutil
import pytest
import filecmp
import numpy as np
from astropy.table import Table

# sys.path.append("..")
from pyHiM import main


"""Check the non regression of modules called inside the jupyter notebook on pyHiM documentation"""


################################# FUNCTIONS #################################

def extract_files(root: str):
    """Extract recursively file informations of all files into a given directory.
    Note: filename is the name without extension

    Parameters
    ----------
    root : str
        The name of root directory

    Returns
    -------
    List[Tuple(str,str,str)]
        List of file informations: (filepath, filename, extension)
    """
    files = []
    # Iterate into dirname and each subdirectories dirpath, dirnames, filenames
    for dirpath, dirnames, filenames in os.walk(root):
        for filename in filenames:
            split_filename = filename.split(".")
            if len(split_filename) > 1:
                extension = split_filename.pop()
            else:
                extension = None
            short_filename = ".".join(split_filename)
            filepath = os.path.join(dirpath, filename)
            files.append((filepath, short_filename, extension))

        if len(dirnames) > 0:
            print("[INFO] Inside:")
            print(dirpath)
            print("[INFO] Subdirectories detected:")
            print(dirnames)

    return files

def compare_npy_files(first_file, second_file, shuffled_plans=False):
    """Load both files as numpy array and compare them.

    Args:
        first_file (string): first complete file path
        second_file (string): second complete file path

    Returns:
        bool: True if they are the same array with the same value
    """
    first_npy = np.load(first_file)
    second_npy = np.load(second_file)
    is_same = np.array_equal(first_npy, second_npy, equal_nan=True)
    if not is_same and shuffled_plans:
        reverse_first = np.swapaxes(first_npy,0,2)
        reverse_second = np.swapaxes(second_npy,0,2)
        is_same = True
        for plan in reverse_first:
            is_inside = (np.equal(plan,reverse_second) | np.isnan(reverse_second)).all((1,2)).any()
            is_same = is_same and is_inside
    return is_same

def compare_ecsv_files(first_file, second_file, columns_to_remove=[], shuffled_lines=False):
    first_ecsv = Table.read(first_file, format="ascii.ecsv")
    second_ecsv = Table.read(second_file, format="ascii.ecsv")

    for col in columns_to_remove:
        first_ecsv.remove_column(col)
        second_ecsv.remove_column(col)
    first_npy = first_ecsv.as_array()
    second_npy = second_ecsv.as_array()
    is_same = True
    if shuffled_lines:
        for line in first_npy:
            if not line in second_npy:
                print(f"SHUFFLE: At line {line}\n from {first_file}\n\n")
                is_same = False
                break
    else:
        comparison = first_npy == second_npy
        is_same = comparison.all()
    return is_same

def compare_line_by_line(first_file, second_file, shuffled_lines=False):
    with open(first_file) as f1:
        with open(second_file) as f2:
            f1_lines = f1.read().splitlines()
            f2_lines = f2.read().splitlines()
            f1_length = len(f1_lines)
            f2_length = len(f2_lines)
            line_index = 0
            is_same = f1_length == f2_length
            while is_same and (line_index < f1_length):
                if shuffled_lines:
                    is_same = f1_lines[line_index] in f2_lines
                    if not is_same:
                        print(f"SHUFFLE: At line number {line_index}\n from {first_file}\n{f1_lines[line_index]}\n")
                else:
                    is_same = f1_lines[line_index] == f2_lines[line_index]
                    if not is_same:
                        print(f"At line number {line_index}\n from {first_file}\n{f1_lines[line_index]}\n from {second_file}\n{f2_lines[line_index]}")
                line_index += 1
            return is_same


################################# TESTS #################################

# build a temporary directory
tmp_dir = tempfile.TemporaryDirectory()
tmp_resources = os.path.join(tmp_dir.name, "resources")
shutil.copytree("pyhim-small-dataset/resources", tmp_resources)
tmp_small_inputs = os.path.join(tmp_resources, "small_dataset/IN")
tmp_stardist_basename = os.path.join(tmp_resources, "stardist_models")
tmp_traces_inputs = os.path.join(tmp_resources, "traces_dataset/IN")


def test_make_projections():
    """Check makeProjections"""
    main(["-F", tmp_small_inputs, "-C", "makeProjections", "-S", tmp_stardist_basename])
    tmp_z_project = os.path.join(tmp_small_inputs, "zProject")
    out_z_project = "pyhim-small-dataset/resources/small_dataset/OUT/makeProjections/"
    out_files = extract_files(out_z_project)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_z_project, filename)
        out_file = os.path.join(out_z_project, filename)
        assert compare_npy_files(tmp_file, out_file)

def test_align_images():
    """Check alignImages"""
    main(["-F", tmp_small_inputs, "-C", "alignImages", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    out_align_images = "pyhim-small-dataset/resources/small_dataset/OUT/alignImages/"
    out_files = extract_files(out_align_images)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        if extension:
            filename = short_filename + "." + extension
        else:
            filename = short_filename
        tmp_file = os.path.join(tmp_align_images, filename)
        out_file = os.path.join(out_align_images, filename)
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file)
        else:
            assert compare_line_by_line(tmp_file, out_file, shuffled_lines=True)


def test_applies_registrations():
    """Check appliesRegistrations"""
    main(["-F", tmp_small_inputs, "-C", "appliesRegistrations", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    out_align_images = "pyhim-small-dataset/resources/small_dataset/OUT/appliesRegistrations/"
    out_files = extract_files(out_align_images)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_align_images, filename)
        out_file = os.path.join(out_align_images, filename)
        assert compare_npy_files(tmp_file, out_file)

def test_align_images_3d():
    """Check alignImages3D"""
    main(["-F", tmp_small_inputs, "-C", "alignImages3D", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    out_align_images = "pyhim-small-dataset/resources/small_dataset/OUT/alignImages3D/"
    out_files = extract_files(out_align_images)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_align_images, filename)
        out_file = os.path.join(out_align_images, filename)
        assert compare_line_by_line(tmp_file, out_file, shuffled_lines=True)

def test_segment_masks_3d():
    """Check segmentMasks3D"""
    main(["-F", tmp_small_inputs, "-C", "segmentMasks3D", "-S", tmp_stardist_basename])
    tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects")
    out_segmented_objects = "pyhim-small-dataset/resources/small_dataset/OUT/segmentMasks3D/"
    out_files = extract_files(out_segmented_objects)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_segmented_objects, filename)
        out_file = os.path.join(out_segmented_objects, filename)
        assert compare_npy_files(tmp_file, out_file)


# def test_segment_sources_3d():
#     """Check segmentSources3D"""
#     main(["-F", tmp_small_inputs, "-C", "segmentSources3D", "-S", tmp_stardist_basename])
#     tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects")
#     out_segmented_objects = "pyhim-small-dataset/resources/small_dataset/OUT/segmentSources3D/"
#     out_files = extract_files(out_segmented_objects)
#     assert len(out_files) > 0
#     for filepath, short_filename, extension in out_files:
#         filename = short_filename + "." + extension
#         tmp_file = os.path.join(tmp_segmented_objects, filename)
#         out_file = os.path.join(out_segmented_objects, filename)
#         assert compare_ecsv_files(tmp_file, out_file,columns_to_remove= ["Buid","id"], shuffled_lines=True)


def test_build_traces():
    """Check build_traces"""
    main(["-F", tmp_traces_inputs, "-C", "build_traces"])
    tmp_builds_pwd_matrix = os.path.join(tmp_traces_inputs, "buildsPWDmatrix")
    out_builds_pwd_matrix = "pyhim-small-dataset/resources/traces_dataset/OUT/build_traces/"
    out_files = extract_files(out_builds_pwd_matrix)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_builds_pwd_matrix, filename)
        out_file = os.path.join(out_builds_pwd_matrix, filename)
        assert compare_ecsv_files(tmp_file, out_file, columns_to_remove=["Trace_ID"])

def test_build_matrix():
    """Check build_matrix"""
    main(["-F", tmp_traces_inputs, "-C", "build_matrix"])
    tmp_builds_pwd_matrix = os.path.join(tmp_traces_inputs, "buildsPWDmatrix")
    out_builds_pwd_matrix = "pyhim-small-dataset/resources/traces_dataset/OUT/build_matrix/"
    out_files = extract_files(out_builds_pwd_matrix)
    assert len(out_files) > 0
    for filepath, short_filename, extension in out_files:
        filename = short_filename + "." + extension
        tmp_file = os.path.join(tmp_builds_pwd_matrix, filename)
        out_file = os.path.join(out_builds_pwd_matrix, filename)
        if extension == "npy":
            assert compare_npy_files(tmp_file, out_file, shuffled_plans=True)
        elif extension == "ecsv":
            assert compare_line_by_line(tmp_file, out_file)