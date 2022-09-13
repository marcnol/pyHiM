# makeProjections
*Projects 3D images in 2D*
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
Parameters to run this scropt will be read from the ```zProject``` field of ```infoList.json```


|Name|Option|Description|
|:-:|:-:|:-:|
|mode|manual|Assign plans between *zmin* and *zmax* to "zRange"|
||automatic|Estimates the focal plane using the maximum standard deviation plane by plane. Use *zwindows* to set "zRange" arround this focal plane.|
||full|Assign all plans to "zRange"|
||laplacian|Split 3D image into blocks of the size given by *blockSize*. Find Laplacian Variance maximum (blur estimation) for each block in order to estimate the focal plane. Rebuild block-by-block 2D image with optimal focal plane of each block. if *zwindows* option is activated, project each block with MIP option.|
|windowSecurity||Used for *automatic* mode, removes the lowest and highest Z-plans.|
|zwindows| | In automatic mode, selects the number of planes below and above the focal plane to be used for making the projection.
|display| | Saves output 2D projections as png files
|zProjectOption|sum|Sum plans in "zRange"|
||MIP|Maximum Intensity Projection of plans in "zRange"|    
|zmax| | Select ending plane to use for projection
|zmin| | Select starting plane to use for projection

## Description


This function will take 3D stacks and project them into 2D.

There are many choices of how to do this:

-   `manual`: indicate the planes in zmin and zmax and set to manual.
    
-   `automatic`: the function estimates focal plane using the maximum of the std deviation from plane to plane, then projects around `zwindows` of the focal plane. Set to automatic.
    
-   `full`: projects all planes into a 2D image. Set to full.
    
    There are some additional options that can be provided to indicate how projections are made:
    
-   `laplacian`: breaks the image into blocks of size `blockSize`. Then calculates the laplacian variance in each block, and estimates the focal position per block by maximizing the laplacian variance. The overall focal plane for the image will be outputed to the terminal and to the block image (see title in image below). The 2D image is reconstructed block by block by using the optimal focal plane for each block. If the parameter `zwindows` is set to zero, then only the image at the focal point will be used. Otherwise we will do an MIP in the subvolume: `focalPlane-zwindows/2:focalPlane+zwindows/2`.Set to laplacian.
    
    There are some additional options that can be provided to indicate how projections are made:
    
-   `windowSecurity`: during automatic focal plane search, it will discard maxima located this number of planes away from the border.
    
-   `zProjectOption`: how it converts a 3D stack into a 2D projection:
    
    -   sum: sums all planes
    -   MIP: maximum intensity projection

## (Invoke)
To run this function exclusively, run *pyHiM* using the ``` -C makeProjections ``` argument. This routine take all 3D images and project its in 2D. Depending on the chosen *mode*, this module start to find the good set of Z-plans, where there is the least noise. This step give a range centered on a focal plan, named *zRange*. After, projection is done on this range either by sum or by maximum intensity projection.

## Graph


```{mermaid}
flowchart TD

		A1[["zProjectionRange()"]] --> A2
		A1 --> A3
		A1 --> A4
		A1 --> A5
		A2[manual] ---> A7
		A3[automatic] --> A6
		A4[full] ---> A7
		A5[laplacian] --> A8
		A6[["calculate_zrange()"]] --> A7
		A7([zRange]) --> A9
		A8[["reinterpolatesFocalPlane()"]] ----> A10
		A8 ---> A9
		A9[["projectsImage2D()"]] --> A14
		A9 --> A15
		subgraph OUTPUT
			A10([focalPlaneMatrix])
			subgraph common
				A11([data_2D])
				A12(["focusPlane"])
				A13([zRange])
			end
		end
		A14[MIP] --> common
		A15[sum] --> common
		subgraph INPUT
			K([3D_DAPI])
			M([3D_Fiducials])
			N([3D_RNA])
			L([3D_Barcodes])
		end
		INPUT --> A1

```