# alignImages
*Registers fiducials using a barcode as reference*

There are several ways of correcting for drift within *pyHiM*:

2.1 **Global drift correction by cross-correlation.** This option just runs a x-correlation between the 2D projected images for the reference and cycle _fiducials. It is the fastest, but will ignore local deformations in the sample and, sometimes, can get fooled by bright dirt in the image that will drive the x-correlation to the wrong place. If your sample is clean and does not show much deformation, this is the way to go. The method will output overlap images that should be used whether the method worked as expected, and to what extent a local correction is needed._

2.2 **Block drift correlation.** This option will also use the 2D projection images of reference and cycle _fiducials, but it will first break them up into blocks and will perform a block-by-block optimization of XY drift. This means that this method is very robust and is not easily fooled by dirt in the sample. However, the method will find a consensus global drift that will be applied and therefore local drift issues are not solved. An additional advantage to method 1 is that it can estimate how much local drift is present in each block and will use this to discard blocks where the local drift is higher than a user-provided tolerance (see below). After you run this method, you will get the uncorrected and corrected images so you can evaluate whether it worked properly and whether local drift correction methods need to be applied. 

2.3 **2D Local drift correction.** This method will be applied after methods 2.1 and 2.2. It will iterate over the DAPI masks detected in the segmentation function (see below), extract a 2D region around each mask, and x-correlate the reference and cycle _fiducials in this 2D sub-region. Thus, this method is slower than methods 1 and 2, but provides for local corrections that account for deformations of the sample. The method will output images with the uncorrected and corrected overlaps for each DAPI mask sub-region so you can evaluate its performance. 

## Invoke

To run this function exclusively, run *pyHiM* using the ``` -C alignImages ``` argument. 
In the set of *fiducial* images, one is chosen by initialization parameters to be the reference. 
The algorithm takes images one by one and align with the reference.
There are several ways to compute the shift:
- Global alignement make simple cross-correlation with two images
- Split image in block and make cross-correlation block by block. The ```alignByBlock``` parameter in the ```alignImages``` field of ```infoList.json``` should be set to ```True```. It calculates the optimal shift between fiducial and reference in each block. It estimates the root mean squared error (RMS) between the reference and the shifted image for each block, and uses the blocks in which the RMS is within ```tolerance```. Mean and standar deviation of the XY shifhs are calcualted, and mean shifts are used for shifting the image and getting the final RMS error. This method is more robust against a bright noise spot.
- Local drift correction in 2D using a bounding box that is ```bezel``` pixels larger than the mask for both the reference fiducial and the fiducial of each cycle. It applies the same cross-correlation algotrithm as before to find aditional local shift. If this shift is larger than ```localShiftTolerance``` in any direction, it will not apply it. 


## Relevant options
Parameters for this script will be read from the  ```alignImages``` field of ```infoList.json```

|Name|Option|Description|
|:-:|:-:|:-:|
|referenceFiducial| |Selects reference barcode image|
|alignByBlock| | Set to false if a block correction is not needed. Default: True|
|bezel| |Selects number of pixels around the fiducial mask for local shift correction|
|localShiftTolerance | | Maximal tolerance in pixels to apply local correction |


```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([2D_Fiducials])
			A2([imReference])
		end
		subgraph C[OUTPUT]
			C1([alignmentResultsTable.ecsv])
			C2(["shift[X][Y]"])
			C3([image2_corrected_raw])
		end
		A1 --> B1
		A2 --> B1
		B1[["align2Files()"]] --> B2
		B1 --> B3
		B2[global alignment] --> B5
		B3[alignByBlock] --> B4
		B4[["alignImagesByBlocks()"]] --> B6
		B5[["align2ImagesCrossCorrelation()"]] --> B6
		B6([shift]) --> B7
		B6 --> C1
		B6 --> C2
		B7[["shiftImage()"]] --> C3

```