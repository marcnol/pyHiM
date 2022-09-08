# segmentMasks
*Segments DAPI and sources in 2D*

## Invoke
To run a 2D segmentation exclusively, run *pyHiM* using the ``` -C segmentMasks ``` argument. This fucnction will be applied when the parameter ```operation``` is set to ```2D```, in the section ```segmentedObjects``` of ```infoList.json```. 


## Relevant options
|Name|Option|Description|
|:-:|:-:|:-:|
|operation|2D| Select 2D mask segmentation|
||3D| Select 3D mask segmentation|
|background_method|inhomogeneous| |
||flat| | 
||stardist| |
|stardist_network| | Name of the network used for segmentation|
|stardist_basename| | Folder containing AI models|
|background_sigma| | Used to remove inhomogeneous background. Default: 3.0|
|threshold_over_std| | Threshold used to detect sources. Default: 1.0|
|area_min| | Minimal area to keep object|
|area_max| | Maximal area to keep object|
|residual_max| | Maximum difference between axial spot intensity and gaussian fit| 

## Outputs

A 2D mask segmentation produces two outputs saved in the `segmentedObjects` folder:

```
scan_002_mask0_002_ROI_converted_decon_ch01_segmentedMasks.png
scan_002_mask0_002_ROI_converted_decon_ch01_Masks.npy
```

The PNG file is a representation of the raw image and the segmented objects.

The NPY file is a 2D labeled numpy array containing the segmented objects with a size identical to the original image. Background has the value _0_ and then each mask contains a different integer. The maximum value in this matrix is identical to the number of masks detected. The file name is constructed using the original root filename with the tag `_Masks`.

_Warning_: This mode operates in 2D, therefore the Startdist network provided **must be** in 2D.
```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([aligned_2D_DAPI])
			A2([aligned_2D_Barcodes])
		end
		subgraph C[OUTPUT]
			C1([Table with barcode coordinates])
			C2([Numpy array with segmented masks])
		end
		A --> B1
		B1[["segmentMasks()"]] --> B2
		B2[["makesSegmentations()"]] --> B3
		B2 --> B4
		B3[barcode] --> B5
		B3 --> B6
		B4[DAPI] --> B9
		B4 --> B10
		B4 --> B11
		B4 ----> B14
		B5[flat] --> B7
		B6[inhomogeneous] --> B8
		B7[["segmentSourceFlatBackground()"]] ---> C
		B8[["segmentSourceInhomogBackground()"]] ---> C
		B9[flat] --> B12
		B10[inhomogeneous] --> B12
		B11[stardist] --> B13
		B12[["segmentMaskInhomogBackground()"]] ---> C
		B13[["segmentMaskStardist()"]] ---> C
		B14[tesselation] --> C
		
	
```