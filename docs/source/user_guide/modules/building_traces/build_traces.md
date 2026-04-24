# build_traces

## Invoke
Inside the folder with your input data, run:
```shell
pyhim -C build_traces
```

![build_traces](../../../_static/from_tuto/build_traces.png)

## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|parameters.json|1|Yes|Parameter file.|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
||||


## Relevant options

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```parameters.json```.

```
"tracing_method": ["masking","clustering"], # list of methods it will use
"mask_expansion": 8,# number of pixels masks will be expanded to assign localizations
"masks2process":{"nuclei":"DAPI","mask1":"mask0"}, # masks identities to process
"KDtree_distance_threshold_mum": 1, # threshold distance used for KDtree clustering
"barcode_coordinates_BEDfile": "", # optional BED file used to annotate barcodes with genomic coordinates
```

### Optional BED-driven barcode annotation

If `barcode_coordinates_BEDfile` is set to a valid BED file path, `build_traces` will load that file and use it to enrich each saved trace entry with genomic metadata.

Example in `parameters.json` (`buildsPWDmatrix` section):

```json
"buildsPWDmatrix": {
  "tracing_method": ["masking", "clustering"],
  "mask_expansion": 8,
  "masks2process": {"nuclei": "DAPI"},
  "KDtree_distance_threshold_mum": 1,
  "barcode_coordinates_BEDfile": "resources/barcode_coordinates.bed"
}
```

When a BED file is provided:

- Barcode identities are validated against the BED mapping.
- If a 5th BED column is present, barcode IDs can be remapped before writing the trace table.
- `Chrom`, `Chrom_Start`, and `Chrom_End` columns in the output trace table are populated from BED coordinates.

### BED file format expected by `build_traces`

Use a tab-separated BED-like file (no header), with at least 4 columns:

1. `chrom` (e.g. `chr2L`, `chrX`)
2. `start` (genomic start, basepairs)
3. `end` (genomic end, basepairs)
4. `barcode` (barcode ID expected in localizations)
5. *(optional)* `barcode_remap` (barcode ID written to output table after remapping)

Example:

```text
chr2L	150000	151000	1
chr2L	220000	221000	2	102
chr2L	300000	301000	3
```

In this example, barcode `2` is remapped to `102` when traces are saved, while all entries receive genomic coordinates from columns 1 to 3.

Output images:

- `_XYZ_ROI*.png`

|  | full image | zoomed images |
| --- |   ---- | --- |
| 3D **mask** | ![image-20220210221402082](../../../_static/user_guide/image-20220210221402082.png) |![image-20220210221430543](../../../_static/user_guide/image-20220210221430543.png)|
| 3D **mask** | ![image-20220210222233148](../../../_static/user_guide/image-20220210222233148.png) |![image-20220210222354093](../../../_static/user_guide/image-20220210222354093.png)|
| 3D **KDtree** |  ||
