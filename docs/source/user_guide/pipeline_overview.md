# Pipeline overview
## Default *pyHiM* flow

To run default pipeline, *pyHiM* need two kinds of data:
- A dictionary of initialization parameters, named `infoList.json`
- 3D images with TIFF format. Four types of images are accepted and will be processed in the following order:
	1. Fiducial
	2. Barcode
	3. Mask (like DAPI)
	4. RNA (optional)

These types of images are called labels. **Note that labels 1,2 & 3 are mandatory for running the default analysis pipeline.**

The default pipeline consists in a sequence of seven basic routines:

1. **makeProjections**: Project all 3D images in 2D 
2. **alignImages**: Compute the best shift to align all 2D fiducials
3. **appliesRegistrations**: Shift 2D barcodes, masks and RNA according to the transformation computed at the alignImages step
4. **alignImages3D**: Take 2D aligned fiducial images and find the best shift along the Z-axis. This shift will be applied on the 3D segmented barcodes at buildHiMmatrix step.
5. **segmentMasks**: Segments 2D aligned barcodes and masks
6. **segmentSources3D**: Applies 2D shift, computed at alignImages step, to 3D barcodes. Then, segments them in 3D.
7. **buildHiMmatrix**: Filter the segmentation results, associate barcode coordinates with the right mask and calculate the pairwise distance (PWD) matrix for each mask.

The detailed pipeline organization is summarized in the table below:

|Command|Fiducial|Barcode|DAPI|RNA|
|:-:|:-:|:-:|:-:|:-:|
|makeProjections|1|4|8|12|
|alignImages|2||||
|appliesRegistrations||5|9|13|
|alignImages3D|3||||
|segmentMasks||6|10||
|segmentSources3D||7|||
|buildsPWDmatrix|||11||

## Input / Output data

Here is a table summarizing the type of input and output data for each routine:

|Routine|Input|Output|
|:-:|---|---|
|**makeProjections**|3D_raw.tif|2D_raw.npy|
|**alignImages**|2D_raw.npy + 2D_reference.npy|alignImages.ecsv|
|**AppliesRegistrations**|2D_raw.npy + alignImages.ecsv|2D_registered.npy|
|**alignImages3D**|3D_raw.tif + 3D_reference.tif + alignImages.ecsv|alignImages_block3D.ecsv|
|**segmentMasks**|2D_registered.npy|segmented_barecode.ecsv + segmented_mask.npy|
|**segmentSources3D**|3D_raw.tif + alignImages.ecsv|3D_segmented_barcode.ecsv|
|**buildHiMmatrix**|segmented_barcode.ecsv + segmented_mask.npy + alignImages_block3D.ecsv + 3D_segmented_barcode.ecsv|PWDMatrix.ecsv + 3D_PWDMatrix.ecsv|

## Processing flowchart for different data sources

Different data sources exist and need to processed differently by *pyHiM*:
- fiducials
- DNA-FISH barcodes
- masks (e.g. DAPI, cell, genomic region)
- RNA-FISH images

Each data source is handled by a different workflow in *pyHiM*. In each workflow, you will see different symbols for `modules`, `input data` and `I/O data`, as follows:

```{mermaid}
flowchart TD
	subgraph graph legend
		Z1((module))
		Z2([initial data])
		Z3[I/O data]
	end 
```

### fiducial flow
This scheme shows the steps involved in the analysis of fiducial images.

```{mermaid}
flowchart
	
	Fidu([3D_Fiducials])

	routine1bis((makeProjections))
	routine2((alignImages))
	routine4((alignImages3D))
	
	zProj0bis[3D_raw.tif]
	zProj1bis[2D_raw.npy]
	zProj2[2D_reference.npy]
	zProj3[3D_reference.npy]
	
	align1[alignImages.ecsv]
	align2[alignImages_block3D.ecsv]

	Fidu --> zProj0bis
	zProj0bis -->
	routine1bis --> 
	zProj1bis --> zProj2
	zProj1bis & zProj2 --> 
	routine2  --> align1
	
	zProj0bis --> zProj3
	zProj3 --> routine4
	zProj0bis --> routine4
	align1 -.-> routine4
	routine4 --> align2

```

### barcode flow
This scheme shows the steps involved in the analysis of DNA-FISH barcode images.

```{mermaid}
flowchart
	
	Barc([3D_Barcodes])

	routine1((makeProjections))
	routine3((appliesRegistrations))
	routine5((segmentMasks))
	routine6((segmentSources3D))
	
	zProj0[3D_raw.tif]
	zProj1[2D_raw.npy]
	
	align1[alignImages.ecsv]
	align3[2D_registered.npy]
	
	seg2[segmentedObjects_barcode.ecsv]
	seg3[segmentedObjects_3D_barcode.ecsv]

	align1 --> routine3
	
	Barc -->
	zProj0 -->
	routine1 -->
	zProj1 --> routine3
	
	routine3 --> align3
	
	
	align3 --> 
	routine5 --> seg2
	
	align1 --> routine6
	Barc --> routine6
	routine6 --> seg3
	
	


```
### masks flow
This scheme shows the steps involved in the analysis of masks images.

```{mermaid}
flowchart
	
	DAPI([3D_Masks])

	routine1((makeProjections))
	routine3((appliesRegistrations))
	routine5((segmentMasks))
	routine7((buildHiMmatrix))
	
	zProj0[3D_raw.tif]
	zProj1[2D_raw.npy]
	
	align1[alignImages.ecsv]
	align2[alignImages_block3D.ecsv]
	align3[2D_registered.npy]
	
	seg1[Masks.npy]
	seg2[segmentedObjects_barcode.ecsv]
	seg3[segmentedObjects_3D_barcode.ecsv]
	
	build1[PWDMatrix.ecsv]
	build2[3D_PWDMatrix.ecsv]


	align1 --> routine3
	
	DAPI -->
	zProj0 -->
	routine1 -->
	zProj1 --> routine3
	
	routine3 --> align3
	
	
	align3 --> 
	routine5 --> seg1
	
	
	align2 -.-> routine7
	seg1 & seg2 & seg3 --> routine7
	
	routine7 --> build1 & build2
	routine7 <--> matching
	subgraph matching
		matching1[HiMscMatrix.npy]
		matching2[3D_HiMscMatrix.npy]
		matching3[Nmatrix.npy]
		matching4[3D_Nmatrix.npy]
		matching7[uniqueBarcodes.ecsv]
		matching8[3D_uniqueBarcodes.ecsv]
	end


```
### RNA-FISH image flow
This scheme shows the steps involved in the analysis of RNA-FISH images.

```{mermaid}
flowchart
	
	RNA([3D_RNA])

	routine1((makeProjections))
	routine3((appliesRegistrations))
	
	zProj0[3D_raw.tif]
	zProj1[2D_raw.npy]
	
	align1[alignImages.ecsv]
	align3[2D_registered.npy]

	align1 --> routine3
	
	RNA -->
	zProj0 -->
	routine1 -->
	zProj1 --> routine3
	
	routine3 --> 
	align3 --> 
	rout4((processSNDchannel.py)) -->
	proc1[SNDassignedCells.ecsv]
```
