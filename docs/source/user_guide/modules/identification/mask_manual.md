# mask_manual
*This script is used to mask images manually*.

Often, we need to select traces from a specific region of a FOV, for instance from cells with a specific expression pattern, or from a well-defined morphological region of a tissue. For this, we use `mask_manual`. The user inputs a 2D image in TIF format and defines by hand multiple ROIs. Each will be saved with a different index in an output NUMPY array 2D image file. This image file can then be used with `trace_label` to automatically label the traces that fall within the pattern defined by hand by the user.



## Invoke
To process 
```shell
$ mask_manual.py --input my_image.tif --label mylabel
```



## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|label|1|Yes|string. Do not use special characters|
|input|1|Yes|2D image in tif format|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
|NUMPY 2D array|1|The output labeled image|

