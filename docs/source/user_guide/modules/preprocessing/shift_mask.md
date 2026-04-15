# shift_mask

*Applies XY registration shifts to previously segmented 3D mask volumes.*

## Invoke

Inside the folder with your input data, run:
```shell
pyhim -C shift_mask
```

`shift_mask` is intended to be run after `mask_3d` has produced unregistered masks and after `register_global` has produced XY shifts.

## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|parameters.json|1|Yes|Parameter file.|
|register_global.ecsv|1|Yes|XY alignment table used to shift masks.|
|mask_3d/data/*_3Dmasks_unregistered.npy|1..n|Yes|Unregistered segmented 3D mask arrays to be shifted.|

## Outputs

|Name shape|Quantity|Description|
|---|---|---|
|mask_3d/data/*_3Dmasks.npy|1..n|Registered segmented 3D mask arrays (same NPY 3D labeled format, shifted in XY).|

## Relevant options

- `segmentedObjects/operation` must include `3D`.
- The routine is executed for mask labels (`DAPI` and `mask`).
- XY shifts are read from `register_global.ecsv` for each ROI and cycle label.

## Description

This routine loads segmented 3D mask arrays that were saved as unregistered files and applies the XY drift correction associated with each ROI and cycle.

### Segmented input format expected by `shift_mask`

`shift_mask` searches in:

- `mask_3d/data/`

with this filename pattern:

- `*_3Dmasks_unregistered.npy`

Each file is expected to be a 3D labeled NumPy array (`.npy`) where background is `0` and mask IDs are positive integers.

### Output filename and format

For each input file, the suffix `_unregistered.npy` is removed and replaced by `.npy`.

Example:

- input: `scan_006_DAPI_001_ROI_converted_decon_ch00_3Dmasks_unregistered.npy`
- output: `scan_006_DAPI_001_ROI_converted_decon_ch00_3Dmasks.npy`

Output files are written to the same folder:

- `mask_3d/data/`

and keep the same file format (`.npy`) and dimensionality (3D labeled mask volume).

## Step by step

- Finds candidate files matching `mask_3d/data/*_3Dmasks_unregistered.npy`.
- Filters files by current cycle label.
- Loads the corresponding XY shift for ROI + label.
- Rounds shift values to integer pixels (mask alignment).
- Applies XY shift to all Z planes.
- Saves shifted mask as `*_3Dmasks.npy` in `mask_3d/data/`.
