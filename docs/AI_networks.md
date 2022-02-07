# AI networks

[TOC]

### DAPI 2D

`Folder: /mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks/`

```sh
Model: stardist 2D
Training set: 3 ROIs, DAPI 2D embryos nc14, deconvolved, presegmented using image analysis.
Name: DAPI_2D_stardist_nc14_nrays:64_epochs:40_grid:2
Comment: third network trial. Excellent performance.
```



### DAPI 3D

`Folder: /mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks/`

```sh
Model: stardist 3D
Training set: 3 ROIs, DAPI 3D embryos nc14, deconvolved
Name: DAPI_3D_stardist_20210617_deconvolved
Comment: first network trial.

Model: stardist 3D
Training set: 5 ROIs, DAPI 3D embryos nc14, deconvolved
Name: DAPI_3D_stardist_20210625_deconvolved
Comment: Second network, trained with a larger number of ROIs. Good performance in validation tests.

Model: stardist 3D
Training set: 17 ROIs, DAPI 3D embryos nc14, nc13, nc12, deconvolved
Name: DAPI_3D_stardist_20210720_deconvolved
Comment: Third network, trained with a larger number of ROIs. Better performance than 0625.

Model: stardist 3D
Training set: 17 ROIs, DAPI 3D embryos nc14, nc13, nc12, RAW and deconvolved images.
Name: DAPI_3D_stardist_20210805_mixed
Comment: Third network, trained with a larger number of ROIs. Better performance than 0625.

Model: stardist 3D
Training set: pancreas, deconvolved images.
Name: DAPI_3D_stardist_17032021
Comment: First network for tissues.
```



### Barcodes 2D

`Folder: /mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks/`

```sh
Model: stardist 2D
Training set: 6 ROIs, 3 barcodes, embryos nc14, deconvolved
Name: PSF_2D_stardist_18032021_single_loci
Comment: first network trial. Good performance in validation tests.

Model: stardist 2D
Training set: 6 ROIs, 3 barcodes, embryos nc14, deconvolved
Name: PSF_2D_stardist_19032021_single_loci
Comment: Second network. Good performance in validation tests.
```



### Barcodes 3D

`Folder: /mnt/grey/DATA/users/marcnol/pyHiM_AI_models/networks/`

```sh
Model: stardist 3D
Training set: 178 barcode images 256x256, deconvolved
Name: PSF_3D_stardist_18032021_single_loci
Comment: first network trial. Good performance in validation tests.

Model: stardist 3D
Training set: 178 barcode images 256x256, deconvolved
Name: PSF_3D_stardist_19032021_single_loci
Comment: Second network. Good performance in validation tests. Bright spots are sometimes not segmented properly.

Model: stardist 3D
Training set: 132 barcode images 256x256, simulations of PSF.
Name: PSF_3D_stardist_20210618_simu_deconvolved_thresh_0_01
Comment: Third network. Excellent performance in validation tests. Some bright spots are still not perfectly segmented (very low percentage).
```





