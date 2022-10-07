# register_localizations

## Invoke
Inside the folder with your input data, run:
```shell
pyhim -C alignImages
```

## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|infoList.json|1|Yes|Parameter file.|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
||||

## Relevant options

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

```
"toleranceDrift": 1 # tolerance drift in pixels. Above this value localizations will not be locally registered
```




outputs images:

- `_registered_barcode_stats.png`
- `_registered_barcode_localizations_ROI*.png`

| statistics of registration                                   | localization map                                             |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20220210221852444](../../../_static/user_guide/image-20220210221852444.png) | ![image-20220210221942291](../../../_static/user_guide/image-20220210221942291.png) |
| ![image-20220210222028835](../../../_static/user_guide/image-20220210222028835.png) | ![image-20220210222006297](../../../_static/user_guide/image-20220210222006297.png) |