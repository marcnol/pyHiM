# filter_localizations



**Invoke**

To run this function exclusively, run *pyHiM* using the ```-C filter_localizations``` argument. This function will find and process all the localization files in the `segmentedObjects` folder. To avoid overwriting data, existing files will be renamed with the extension `_version_n` where `n`will be incremented from run to run. The output of `filter_localizations` will be saved with the original localizations filename. A comment in the header will be added to indicate that a *filtering* operation was run on this file.



**Relevant options**

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

```
"flux_min": 10, # maximum flux allowed for 2D
"flux_min_3D": 4000,# maximum flux allowed for 3D
```



**Output images**

- `_filtered_barcode_localizations_ROI*.png`

