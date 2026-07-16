from types import SimpleNamespace

import pytest

pytest.importorskip("astropy")
from astropy.table import Table

from matrixOperations.merge_inputs import MergeInputs


def _write_dummy_table(path, value):
    table = Table({"ROI #": [value]})
    table.write(path, format="ascii.ecsv", overwrite=True)


def test_merge_inputs_creates_expected_outputs(tmp_path):
    data_path = tmp_path
    seg_params = SimpleNamespace(
        localize_2d_folder="localize_2d",
        localize_3d_folder="localize_3d",
        outputFile="localizations",
    )
    reg_params = SimpleNamespace(
        register_local_folder="register_local", outputFile="shifts"
    )

    for folder, suffix, values in (
        (seg_params.localize_3d_folder, "_3D_barcode.ecsv", [1, 2]),
        (reg_params.register_local_folder, "_block3D.ecsv", [3, 4]),
    ):
        for idx, val in enumerate(values):
            file_path = (
                data_path
                / folder
                / "data"
                / f"{seg_params.outputFile if 'barcode' in suffix else reg_params.outputFile}_{idx}{suffix}"
            )
            file_path.parent.mkdir(parents=True, exist_ok=True)
            _write_dummy_table(file_path, val)

    merger = MergeInputs(None)
    merged_local_shifts = merger.merge_all(data_path, seg_params, reg_params)

    merged_localizations = (
        data_path
        / seg_params.localize_3d_folder
        / "data"
        / "localizations_3D_barcode.ecsv"
    )
    merged_registration = (
        data_path / reg_params.register_local_folder / "data" / "shifts_block3D.ecsv"
    )

    assert merged_localizations.exists()
    assert merged_registration.exists()
    assert merged_local_shifts == str(merged_registration)

    stacked_localizations = Table.read(merged_localizations, format="ascii.ecsv")
    stacked_registration = Table.read(merged_registration, format="ascii.ecsv")
    assert len(stacked_localizations) == 2
    assert len(stacked_registration) == 2


def test_merge_inputs_skips_when_outputs_present(tmp_path):
    data_path = tmp_path
    seg_params = SimpleNamespace(
        localize_2d_folder="localize_2d",
        localize_3d_folder="localize_3d",
        outputFile="localizations",
    )
    reg_params = SimpleNamespace(
        register_local_folder="register_local", outputFile="shifts"
    )

    existing_localization = (
        data_path
        / seg_params.localize_3d_folder
        / "data"
        / "localizations_3D_barcode.ecsv"
    )
    existing_localization.parent.mkdir(parents=True, exist_ok=True)
    _write_dummy_table(existing_localization, 10)

    existing_registration = (
        data_path / reg_params.register_local_folder / "data" / "shifts_block3D.ecsv"
    )
    existing_registration.parent.mkdir(parents=True, exist_ok=True)
    _write_dummy_table(existing_registration, 20)

    merger = MergeInputs(None)
    merged_local_shifts = merger.merge_all(data_path, seg_params, reg_params)

    assert merged_local_shifts == str(existing_registration)
    assert Table.read(existing_localization, format="ascii.ecsv")["ROI #"][0] == 10
    assert Table.read(existing_registration, format="ascii.ecsv")["ROI #"][0] == 20
