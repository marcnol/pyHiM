
# alignImages3D
*Aligns fiducials in 3D.*

## Invoke
Inside the folder with your input data, run:
```shell
pyhim -C alignImages3D
```
#image
## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|infoList.json|1|Yes|Parameter file.|
|<image_name>.tif|2..n|Yes|Images with a fiducial channel to align in 3D.|
|alignImages.table|1|No|XY alignment resulting from the XY alignment produced while running `alignImages`.|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
|alignImages_block3D.dat|1|Shift block maps for X, Y and Z and quality matrices.|
|?||Corrected blocks in XY, ZX, ZY|

## Relevant options

To run, the value for ```localAlignment``` key should be ```block3D```. The other parameters are shared with ```alignImages```.

## Description

None of the methods above takes into account the drift of the sample in the z-plane. While this is typically very small given the good performance of autofocus, it could be an issue in some instances. This method will first apply the 2D drift obtained using methods 1 or 2 to the 3D stack of cycle _. Then it will background-substract and level-normalize the reference and cycle _fiducial images and will break them into 3D blocks (somewhat similar to method 2, which was breaking images into 2D blocks). Next, it will x-correlate every single 3D block in the reference image to the corresponding, pre-aligned block in the cycle _image to obtain a local 3D drift correction. The results are outputted as 3 matrices that indicate the correction applied to each block in z, x and y. In addition, a reassembled image made of XY, XZ and YZ projections is outputted to evaluate performance. Needless to say, this is the slowest but most performant method in the stack.

## Step by step

To run this function exclusively, run *pyHiM* using the ``` -C alignImages3D ``` argument.
The following steps are implemented:
- Iterate over reference fiducials available for each ROI
- Iterate over all cycles for a given ROI
- Load 3D reference fiducial image the current imaging cycle
- Re-align 3D fiducial image using XY alignment resulting from the XY alignment produced while running `alignImages`. If this is not available, it will XY project the 3D stack of reference and cycle _fiducial to get an XY realignment. Beware, this will be global and will not use blockAlignment._
- Breaks 3D images for both reference and cycle _fiducials in blocks (defined by `blockSizeXY`)
- Cross-correlates each block to get an XYZ shift. This provides a 3D local drift correction for each block
- Store shifts in output Table that contains values for each block (columns shift_z, shift_x and shift_y).
- Store quality of image superposition based on the normalized root mean square matrix for each block in output Table (columns quality_xy, quality_zy, quality_zx).

## Graph (useless?)

```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([3D_Fiducials])
			A2([reference list])
		end
		subgraph C[OUTPUT]
			C1([ECSV table with all shifts])
			C2([aligned_3D_Fiducials])
		end
		A --> B1
		B1[["alignFiducials3Dfile()"]] --> B3
		B1 --> B2
		B2([shift]) --> B3
		B3[["appliesXYshift3Dimages()"]] --> B4
		B4[["imageBlockAlignment3D()"]] --> C
	
```