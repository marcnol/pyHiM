# segmentSources3D
*Segments sources in 3D*


## Invoke
Inside the folder with your input data, run:
```shell
pyhim -C segmentSources3D
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

## Description

## Graph

```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([3D_Barcodes])
			A2([shift table])
		end
		subgraph C[OUTPUT]
			C1([ECSV table with barcode coordinates in 3D])
		end
		A1 --> B1
		A2 --> B5
		B1[["segmentSources3D_file()"]] --> B2
		B1 --> B3
		B2[stardist] ---> B5
		B3[not stardist] --> B4
		B4[["preProcess3DImage()"]] --> B5
		B5[["appliesXYshift3Dimages()"]] --> B6
		B6([aligned_3D_Barcodes]) --> B7
		B6 --> B8
		B7[["_segments3Dvolumes_StarDist()"]] --> B9
		B8[["_segments3DvolumesByThresholding()"]] --> B9
		B9([segmentedImage3D]) --> B10
		B10[["getMaskProperties()"]] --> B11
		B11([spots]) --> B12
		B12[["bigfish.fit_subpixel()"]] --> B13
		B13([spots_subpixel]) --> C
		
	
```