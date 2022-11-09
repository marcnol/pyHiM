import sys
import os
import tempfile
import shutil
sys.path.append("..")
from pyHiM import main


def test_non_regression_for_tuto_rtd():
    """Check the non regression of modules called inside the jupyter notebook on pyHiM documentation"""

    # build a temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_resources = os.path.join(tmp_dir, "resources")
        shutil.copytree("pyhim-small-dataset/resources", tmp_resources)

        tmp_small_inputs = os.path.join(tmp_resources, "small_dataset/IN")
        cmd_list = "makeProjections,alignImages,appliesRegistrations,alignImages3D,segmentMasks3D,segmentSources3D"
        tmp_stardist_basename = os.path.join(tmp_resources, "stardist_models")
        main(["-F", tmp_small_inputs, "-C", cmd_list, "-S", tmp_stardist_basename])
        
        tmp_zProject = os.path.join(tmp_small_inputs, "zProject")
        assert os.path.isdir(tmp_zProject)

        tmp_traces_inputs = os.path.join(tmp_resources, "traces_dataset/IN")
        main(["-F", tmp_traces_inputs, "-C", "build_traces,build_matrix"])

        tmp_buildsPWDmatrix = os.path.join(tmp_traces_inputs, "buildsPWDmatrix")
        assert os.path.isdir(tmp_buildsPWDmatrix)



# def test_mocking_cmd(mocker):

    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.make_projections',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.align_images',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.apply_registrations',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.align_images_3d',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.segment_masks',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.segment_masks_3d',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.segment_sources_3d',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.filter_localizations',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.register_localizations',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.build_traces',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.build_matrix',
    # )
    # mocker.patch(
    #     'fileProcessing.functionCaller.HiMFunctionCaller.process_pwd_matrices',
    # )


    # main(["-F", "resources/small_dataset", "-C", "makeProjections,alignImages,appliesRegistrations,alignImages3D,segmentMasks3D,segmentSources3D"])

    # main(["-F", "resources/traces_dataset", "-C", "build_traces,build_matrix"])