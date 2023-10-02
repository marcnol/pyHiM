# Building chromatin traces

## Build traces: new method

The new method requires executing several modules:

- `filter_localizations`
- `register_localizations`
- `build_traces`
- `build_matrices`

```{toctree}
:maxdepth: 1

building_traces/filter_localizations
building_traces/register_localizations
building_traces/build_traces
building_traces/trace_selector
building_traces/trace_combinator
building_traces/trace_filter
building_traces/trace_analyser
building_traces/build_matrices
```


## build traces: old method

The old module `buildPWDmatrix `  does all operations at once: filtering, local registration, tracing by masking, and construction of PWD matrix.

**Operation**

The ```processesPWDmatrices``` script:

- iterates over ROIs
  - assigns barcode localizations to DAPI masks
  - applies local drift correction, if available
  - removes localizations using flux and driftTolerance
  - calculates the pair-wise distances for each single-cell mask
  - outputs are:
    - Table with #cell #PWD #coordinates (e.g. ```buildsPWDmatrix_3D_order:0_ROI:1.ecsv```)
    - NPY array with single cell PWD single cell matrices (e.g. ```buildsPWDmatrix_3D_HiMscMatrix.npy```)
    - NPY array with barcode identities (e.g. ```buildsPWDmatrix_3D_uniqueBarcodes.ecsv```)
    - the files with no ```3D``` tag contain data analyzed using 2D localizations.
- Single-cell results are combined together to calculate:
     - Distribution of pairwise distance for each barcode combination
     - Ensemble mean pairwise distance matrix using mean of distribution
     - Ensemble mean pairwise distance matrix using Kernel density estimation
     - Ensemble Hi-M matrix using a predefined threshold
          - For each of these files, there is an image in ```PNG``` format saved. Images containing ```3D``` are for 3D other are for 2D.



**Invoke**

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```parameters.json```.

If you want to run this function exclusively, run *pyHiM* using the ```-C buildHiMmatrix``` argument.

```sh
usage: pyhim [-h] [-F ROOTFOLDER] [-C CMD] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  -C CMD, --cmd CMD     Comma-separated list of routines to run: 
  						makeProjections, appliesRegistrations,
                        alignImages,alignImages3D, segmentMasks,
                        segmentSources3D,buildHiMmatrix
  --threads THREADS     Number of threads to run in parallel mode. If none,
                        then it will run with one thread.
```



**Options**

"buildsPWDmatrix"

```
"folder": "buildsPWDmatrix",  # output folder
"flux_min": 1000, *Description:* minimum flux per spot. If flux is smaller, localization will be discarded
"flux_min_3D": 0.1, *Description:* minimum flux per spot for 3D localizations. If flux is smaller, localization will be discarded
"toleranceDrift":1, *Description*: tolerance used for block drift correction, in px
```


**Outputs**

Mean pairwise distance matrix. By default means are calculated using Kernel Density Estimators of the PWD distributions.

<img src="../../_static/user_guide/buildsPWDmatrix_HiMmatrix.png" alt="buildsPWDmatrix_HiMmatrix" style="zoom:50%;" />

In addition, the function outputs the distribution of distances for each combination of barcodes:

<img src="../../_static/user_guide/buildsPWDmatrix_PWDhistograms.png" alt="buildsPWDmatrix_PWDhistograms" style="zoom: 25%;" />


### Filtering barcode localizations

There are several filters:
1. Properties of 2D localization algorithm (e.g. brightness)

2. Accuracy of 3D localization: sigma of fit, correlation between z-position from weighted moment and from gaussian fit, etc

3. Accuracy of drift correction in the region where the barcode was localized. 

   This is only applied if LocalDrift correction was **not** run. 


*Examples.*

| Filtering | Matrix |
| --- |  ---- |
| Unfiltered matrix. Total barcode localizations: 18700 | <img src="../../_static/user_guide/buildsPWDmatrix.png" alt="buildsPWDmatrix" style="zoom:25%;" /> |
|```toleranceDrift = 1px```. Barcode localizations kept: 12377 of a total: 18700.| <img src="../../_static/user_guide/buildsPWDmatrixFilterBlockDrift.png" alt="buildsPWDmatrixFilterBlockDrift" style="zoom:25%;" />|
| ```toleranceDrift = 1px```  ```Flux = 100```. Barcode localizations kept: 5562 of a total: 18700. | <img src="../../_static/user_guide/buildsPWDmatrix_filterFlux100.png" alt="buildsPWDmatrix_filterFlux100" style="zoom:25%;" /> |
|```toleranceDrift = 1px```  ```Flux = 200```. Barcode localizations kept: 4528 of a total: 18700. | <img src="../../_static/user_guide/buildsPWDmatrix_filterFlux.png" alt="buildsPWDmatrix_filterFlux" style="zoom:25%;" />|
|```toleranceDrift = 1px``` ```Flux = 1000```. Barcode localizations kept: 1923 of a total: 18700.| <img src="../../_static/user_guide/buildsPWDmatrix_filterFlux1000.png" alt="buildsPWDmatrix_filterFlux1000" style="zoom:25%;" />|



### Other outputs

In addition to the PWD matrix, we now also have a map of the alignment accuracy  and scatter plots showing the flux of each barcode, its sharpness, magnitude and roundness. These are used in order to validate the segmentation process and help with the selection of the ```flux``` threshold used in this filtering step.

*Alignment accuracy*

This provides a map of all barcode localizations in an ROI, colorcoded by the accuracy of localization. Colorbar scale is in pixels.

<img src="../../_static/user_guide/BarcodeAlignmentAccuracy_ROI1_2D2.png" alt="BarcodeAlignmentAccuracy_ROI:1_2D2" style="zoom: 50%;" />

*Barcode localization statistics*

This provides  localization statistics from ASTROPY. The main use of these plots is to determine if the threshold ```flux``` used is correct. Default is *200*.

<img src="../../_static/user_guide/BarcodeStats_ROI1_2D.png" alt="BarcodeStats_ROI:1_2D" style="zoom: 67%;" />
