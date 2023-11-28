# mask_cellpose
This script is provided for the user to be able to segment nuclei or other objects (e.g. mask images) using cellpose.





## Invoke
To process 
```shell
$ mask_cellpose.py --input my_image.tif 
```



## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|input|1|Yes|3D image in tif format|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
|NUMPY 3D array|1|The output labeled image|

