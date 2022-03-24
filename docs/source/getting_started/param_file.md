# WIP - What is infoList.json ?

This file contains all input parameters necessary for pyHiM to run as you want.

To know what and how modify in this file for a specific feature go to [user guide section](../user_guide/fundamental.md) of this feature.

You can find a global description of each parameters in [reference guide](../reference/infoList_comprehension.md).

Find below an example extract:
```json
"labels": {
    "DAPI": {
        "order": 3, 
        "acquisition": {
            "label_channel": "ch00",
            "label_channel_fiducial": "ch01"            
        },            
        "segmentedObjects": {
            "area_max": 3000,
            "area_min": 150,
            "background_method": "inhomogeneous",
            "centroidDifference_max": 5,
            "operation": "2D,3D",
            "stardist_network": "DAPI_2D_stardist_nc14_nrays:64_epochs:40_grid:2",
            "stardist_network3D": "DAPI_3D_stardist_20210720_deconvolved",
            "tesselation": true,
            "threshold_over_std": 1.0
        },
        "zProject": {
            "mode": "full",
            "zProjectOption": "sum",
            "zmax": 59,
            "zmin": 1,
            "zwindows": 15
        }
    },
```


**TODO:**
- *a example of a possible change to customize a parameter just for a label*