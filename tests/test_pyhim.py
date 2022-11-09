import sys
import os
import tempfile
import shutil
import pytest

# sys.path.append("..")
from pyHiM import main


"""Check the non regression of modules called inside the jupyter notebook on pyHiM documentation"""


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
    assert os.path.isdir(tmp_z_project)

def test_align_images():
    """Check alignImages"""
    main(["-F", tmp_small_inputs, "-C", "alignImages", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    assert os.path.isdir(tmp_align_images)

def test_applies_registrations():
    """Check appliesRegistrations"""
    main(["-F", tmp_small_inputs, "-C", "appliesRegistrations", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    assert os.path.isdir(tmp_align_images)

def test_align_images_3d():
    """Check alignImages3D"""
    main(["-F", tmp_small_inputs, "-C", "alignImages3D", "-S", tmp_stardist_basename])
    tmp_align_images = os.path.join(tmp_small_inputs, "alignImages")
    assert os.path.isdir(tmp_align_images)

def test_segment_masks_3d():
    """Check segmentMasks3D"""
    main(["-F", tmp_small_inputs, "-C", "segmentMasks3D", "-S", tmp_stardist_basename])
    tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects")
    assert os.path.isdir(tmp_segmented_objects)


def test_segment_sources_3d():
    """Check segmentSources3D"""
    main(["-F", tmp_small_inputs, "-C", "segmentSources3D", "-S", tmp_stardist_basename])
    tmp_segmented_objects = os.path.join(tmp_small_inputs, "segmentedObjects")
    assert os.path.isdir(tmp_segmented_objects)


def test_build_traces():
    """Check build_traces"""
    main(["-F", tmp_traces_inputs, "-C", "build_traces"])
    tmp_builds_pwd_matrix = os.path.join(tmp_small_inputs, "buildsPWDmatrix")
    assert os.path.isdir(tmp_builds_pwd_matrix)

def test_build_matrix():
    """Check build_matrix"""
    main(["-F", tmp_traces_inputs, "-C", "build_matrix"])
    tmp_builds_pwd_matrix = os.path.join(tmp_small_inputs, "buildsPWDmatrix")
    assert os.path.isdir(tmp_builds_pwd_matrix)
