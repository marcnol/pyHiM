"""Utilities to merge localization and registration tables prior to trace building.

This module consolidates multiple per-image outputs from ``SegmentSources3D``
and ``register_local`` into the single-file convention expected by
``build_traces``. When the aggregated file already exists, the merger exits
without modifying the data.
"""

from __future__ import annotations

import glob
import os
from typing import Iterable, Optional

from astropy.table import Table, vstack

from core.parameters import RegistrationParams, SegmentationParams
from core.pyhim_logging import print_log


class MergeInputs:
    """Merge localization and registration tables before trace building."""

    def __init__(self, param):
        self.current_param = param

    def _merge_tables(
        self,
        files: Iterable[str],
        output_path: str,
        description: str,
    ) -> bool:
        """Merge a collection of ECSV tables into a single file.

        Parameters
        ----------
        files
            Iterable of table paths to merge.
        output_path
            Destination file path following the build_traces naming convention.
        description
            Human-friendly label used in log messages.

        Returns
        -------
        bool
            ``True`` when a merged file is written, ``False`` otherwise.
        """

        tables = [Table.read(path, format="ascii.ecsv") for path in files]
        if not tables:
            print_log(f"! No {description} tables found to merge.")
            return False

        merged_table = vstack(tables, metadata_conflicts="silent")

        merged_table.meta["comments"] = [""]

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        merged_table.write(
            output_path,
            format="ascii.ecsv",
            overwrite=True,
        )

        print_log(f"$ Merged {len(tables)} {description} table(s) into: {output_path}")
        return True

    def merge_localization_tables(
        self, data_path: str, seg_params: SegmentationParams
    ) -> None:
        """Merge per-file localization tables into the expected filenames."""

        for folder, suffix in (
            (seg_params.localize_2d_folder, "_2D_barcode.dat"),
            (seg_params.localize_3d_folder, "_3D_barcode.dat"),
        ):
            output_path = os.path.join(
                data_path, folder, "data", f"{seg_params.outputFile}{suffix}"
            )
            if os.path.exists(output_path):
                print_log(
                    f"$ Localization file already present: {output_path}\n> Nothing to merge."
                )
                continue
            pattern = os.path.join(
                data_path, folder, "data", f"{seg_params.outputFile}_*{suffix}"
            )
            files = sorted(glob.glob(pattern))
            self._merge_tables(files, output_path, "localization")

    def merge_registration_tables(
        self, data_path: str, reg_params: RegistrationParams
    ) -> Optional[str]:
        """Merge per-file local registration tables.

        Returns the path of the merged registration file when created, otherwise
        ``None``.
        """

        output_path = os.path.join(
            data_path,
            reg_params.register_local_folder,
            "data",
            f"{reg_params.outputFile}_block3D.dat",
        )
        if os.path.exists(output_path):
            print_log(
                f"$ Registration file already present: {output_path}\n> Nothing to merge."
            )
            return output_path

        pattern = os.path.join(
            data_path,
            reg_params.register_local_folder,
            "data",
            f"{reg_params.outputFile}_*_block3D.dat",
        )
        files = sorted(glob.glob(pattern))
        merged = self._merge_tables(files, output_path, "registration")
        return output_path if merged else None

    def merge_all(
        self,
        data_path: str,
        seg_params: SegmentationParams,
        reg_params: RegistrationParams,
    ) -> Optional[str]:
        """Merge both localization and registration tables.

        Returns the merged registration file path when one is created.
        """

        self.merge_localization_tables(data_path, seg_params)
        return self.merge_registration_tables(data_path, reg_params)
