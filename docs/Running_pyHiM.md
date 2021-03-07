

# Running pyHiM

[TOC]

## 1- Basic run

### Run pyHiM

Ensure you followed the steps described previously during installation when you did the test run:

1. Identify a ```destination_directory``` where your data are stored. The raw deconvolved files can be in your ```destination_directory``` or within a sub-folder.
2. Be aware of not putting more than ONE sub-folder with TIFF files in the ```destination_directory```. If your ```destination_directory``` already has the raw deconvolved TIFFs then remvove any other directory with TIFFs from ```destination_directory```
3. copy files to your ```destination_directory```  (names are self-explanatory)
   1. infoList_DAPI.json
   2. infoList_RNA.json
   3. infoList_fiducial.json
   4. infoList_barcode.json
4. Change the fiducial RT by running ```changeRT_infoList.py``` at the command line in the ```destination_directory```. The input arguments are the RT currently present in the infoList files and the RT that you want to change it for. For instance: ```changeRT_infoList.py RT33 RT95```changes RT33 to RT95 in all the infoList files.
5. Run pyHiM by the following command at the command line:

```bash
pyHiM.py
```

This assumes that you are running it from the ```destination_directory```. If it is not the case, use the ``-F`` flag with the directory with your data.



The arguments are

```sh
usage: pyHiM.py [-h] [-F ROOTFOLDER] [--parallel] [--localAlignment] [--refit]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --parallel            Runs in parallel mode
  --localAlignment      Runs localAlignment function
  --refit               Refits barcode spots using a Gaussian axial fitting function.

```



```-F ``` indicates the rootFolder where pyHiM expects to find the dataset.

```--parallel``` flag will make it run in parallel mode. Be ready to open your browser in ```http://localhost:8787``` and make sure you connect by ```ssh -L 8787:localhost:8787 marcnol@lopevi```. Change your username and server of course...

```---localAlignment``` will run the local alignment function. See below

```---refit```  will refit barcode spots using a Gaussian axial fitting function.



#### infoList parameters files

a typical file (DAPI example) looks like:

```bash
{
    "acquisition": {
        "DAPI_channel": "ch00",
        "RNA_channel": "ch01",
        "fileNameRegExp": "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\w|-]+).tif",
		"barcode_channel": "ch01",
        "fiducialBarcode_channel": "ch00",
        "fiducialDAPI_channel": "ch02",
        "label": "DAPI",
        "pixelSizeXY": 0.1,
        "pixelSizeZ": 0.25,        
        "positionROIinformation": 3
    },
    "alignImages": {
        "folder": "alignImages",
        "operation": "overwrite",
        "outputFile": "alignImages.bed",
        "alignByBlock": true,
        "tolerance": 0.1,
        "localAlignment": "overwrite",
        "lower_threshold": 0.999, 
        "higher_threshold": 0.9999999, 
		"localShiftTolerance": 1,
        "background_sigma": 3.0,  
        "bezel": 20,               
        "referenceFiducial": "RT27"
    },
    "projectsBarcodes": {
        "folder": "projectsBarcodes",
        "operation": "overwrite",
        "outputFile": "projectsBarcodes"
    },
    "buildsPWDmatrix": {
        "folder": "buildsPWDmatrix",  # output folder
        "flux_min": 200  # min flux to keeep object                
        "toleranceDrift":1,
    },    
    "segmentedObjects": {
        "area_max": 3000,
        "area_min": 150,
        "background_method": "stardist",
        "stardist_network": "stardist_nc14_nrays:128_epochs:400_grid:2",
        "stardist_basename": "/mnt/grey/DATA/users/marcnol/models",
        "background_sigma": 3.0,
        "folder": "segmentedObjects",
        "fwhm": 3.0,
        "brightest": 1100,
        "intensity_max": 59,
        "intensity_min": 0,
        "operation": "overwrite",
        "outputFile": "segmentedObjects",
        "residual_max": 2.5,
        "sigma_max": 5,
        "centroidDifference_max": 5,       
        "3Dmethod":"zASTROPY",
        "3DGaussianfitWindow": 3,
        "threshold_over_std": 1.0,
        "3dAP_window": 5,
        "3dAP_flux_min": 2,
        "3dAP_brightest": 100,
        "3dAP_distTolerance": 1
    },
    "zProject": {
        "display": true,
        "folder": "zProject",
        "mode": "full",
        "operation": "skip",
        "saveImage": true,
        "windowSecurity": 2,
        "zProjectOption": "sum",
        "zmax": 59,
        "zmin": 1,
        "zwindows": 10
    }
}
```



Here are some options for the different parameters and a brief description

"acquisition"
```
"DAPI_channel": "ch00",
"RNA_channel": "ch01",
"fileNameRegExp": "scan_(?P<runNumber>[0-9]+)_(?P<cycle>[\\w|-]+)_(?P<roi>[0-9]+)_ROI_converted_decon_(?P<channel>[\\>. This regular expression encodes our filename format
"barcode_channel": "ch01",
"fiducialBarcode_channel": "ch00",
"fiducialDAPI_channel": "ch02",
"label": "DAPI", *Options: DAPI, fiducial, barcode, RNA
"pixelSizeXY": 0.1, *lateral pixel size in nm*
"pixelSizeZ": 0.25 *axial pixel size in nm*
"positionROIinformation": 3 *position for the ROI in the filename: will be removed in future versions!*
```

"zProject"

```
"folder": "zProject",  *Description:* output folder
"operation": "overwrite",  *Options:* overwrite | skip
"mode": "full",  *Options:* full | manual | automatic | laplacian
"display": True,
"blockSize": 128,
"saveImage": True,
"zmin": 1,
"zmax": 59,
"zwindows": 10,
"windowSecurity": 2,
"zProjectOption": "sum",  *Options:* **sum** | **MIP**
```

"alignImages"
```
"folder": "alignImages",  *Description:* output folder
"operation": "overwrite",  *Options:* overwrite | skip
"outputFile": "alignImages",
"referenceFiducial": "RT18"
"alignByBlock": True, # alignByBlock True will perform block alignment
"tolerance": 0.1, #Used in blockAlignment to determine the % of error tolerated
"lower_threshold": 0.999, # lower threshold to adjust image intensity levels before xcorrelation
"higher_threshold": 0.9999999, # higher threshold to adjust image intensity levels before xcorrelation
"background_sigma": 3.0,  # used to remove inhom background
"localShiftTolerance": 1,
"bezel": 20,
```

"projectsBarcodes"
```
"folder": "projectsBarcodes",  *Description:* output folder
"operation": "overwrite",  *Options:* overwrite | skip
"outputFile": "projectsBarcodes",
```

"buildsPWDmatrix"

```
"folder": "buildsPWDmatrix",  # output folder
"flux_min": 1000, *Description:* minimum flux per spot. If flux is smaller, localization will be discarded
"toleranceDrift":1, *Description*: tolerance used for block drift correction, in px
```

"segmentedObjects"

```
"folder": "segmentedObjects",  *Description:* output folder
"operation": "overwrite",  *Options:* overwrite | skip
"outputFile": "segmentedObjects",
"background_method": "inhomogeneous",  *Options:* **flat** |**inhomogeneous** | **stardist** (AI)
"stardist_network": "stardist_nc14_nrays:64_epochs:20_grid:2", *Description*: name of network
"stardist_basename": "/mnt/grey/DATA/users/marcnol/models", *Description*: location of AI models
"background_sigma": 3.0,  *Description:* used to remove inhomogenous background
"threshold_over_std": 1.0,  *Description:* threshold used to detect sources
"fwhm": 3.0,  *Description:* source size in pixels
"brightest": 1100,  *Description:* max number of objects segmented per FOV (only for barcodes!)
"intensity_min": 0,  *Description:* min intensity to keep object
"intensity_max": 59,  *Description:* max intensity to keep object
"area_min": 50,  *Description:* min area to keep object
"area_max": 500,  *Description:* max area to keep object
"residual_max": 2.5, *Description:*  maximum difference between axial spot intensity and gaussian fit.
"3Dmethod":"zASTROPY", # options: zASTROPY, zProfile
"sigma_max": 5,*Description:* maximum gaussian fit sigma allowed (axial spot intensity)
"centroidDifference_max": 5,  *Description:* max difference between z centroid position determined by moment and by gaussian fitting       
"3DGaussianfitWindow": 3,*Description:* size of window in xy to extract 3D subVolume, in px. 3 means subvolume will be 7x7.
"3dAP_window": 5, # constructs a YZ image by summing from xPlane-window:xPlane+window
"3dAP_flux_min": 2, # # threshold to keep a source detected in YZ
"3dAP_brightest": 100, # number of sources sought in each YZ plane
"3dAP_distTolerance": 1, # px dist to attribute a source localized in YZ to one localized in XY
```



#### MakeProjections

This function will take 3D stacks and project them into 2D.

There are many choices of how to do this:

- ```manual```: indicate the planes in zmin and zmax and set <mode> to <u>manual</u>.

- ```automatic```:  the function estimates focal plane using the maximum of the std deviation from plane to plane, then projects around ```zwindows``` of the focal plane. Set <mode> to <u>automatic</u>.

- ```full```: projects all planes into a 2D image.  Set <mode> to <u>full</u>.

  There are some additional options that can be provided to indicate how projections are made:

- ```laplacian```: breaks the image into blocks of size ```blockSize```. Then calculates the laplacian variance in each block, and estimates the focal position per block by maximizing the laplacian variance. The overall focal plane for the image will be outputed to the terminal and to the block image (see title in image below). The 2D image is reconstructed block by block by using the optimal focal plane for each block. If the parameter ```zwindows``` is set to zero, then only the image at the focal point will be used. Otherwise we will do an MIP in the subvolume: ``` focalPlane-zwindows/2:focalPlane+zwindows/2```.Set <mode> to <u>laplacian</u>.

  

  There are some additional options that can be provided to indicate how projections are made:

- ```windowSecurity```: during automatic focal plane search, it will discard maxima located this number of planes away from the border.

- ```zProjectOption```: how it converts a 3D stack into a 2D projection:
  - sum: sums all planes
  - MIP: maximum intensity projection

"zProject"

```
"folder": "zProject",  *Description:* output folder
"operation": "overwrite",  *Options:* overwrite | skip
"mode": "full",  *Options:* full | manual | automatic | laplacian
"display": True,
"blockSize": 128,
"saveImage": True,
"zmin": 1,
"zmax": 59,
"zwindows": 10,
"windowSecurity": 2,
"zProjectOption": "sum",  *Options:* **sum** | **MIP**
```





#### Drift Correction



There are several ways of correcting for drift within pyHiM:

1. **Global drift correction by cross-correlation.** This option just runs a x-correlation between the 2D projected images for the reference and cycle <i>  fiducials. It is the fastest, but will ignore local deformations in the sample and, sometimes, can get fooled by bright dirt in the image that will drive the x-correlation to the wrong place. If your sample is clean and does not show much deformation, this is the way to go. The method will output overlap images that should be used whether the method worked as expected, and to what extent a local correction is needed.
2. **Block drift correlation.** This option will also use the 2D projection  images of reference and cycle <i> fiducials, but it will first break them up into blocks and will perform a block-by-block optimization of XY drift. This means that this method is very robust and is not easily fooled by dirt in the sample. However, the method will find a consensus global drift that will be applied and therefore local drift issues are not solved. An additional advantage to method 1 is that it can estimate how much local drift is present in each block and will use this to discard blocks where the local drift is higher than a user-provided tolerance (see below). After you run this method, you will get the uncorrected and corrected images so you can evaluate whether it worked properly and whether local drift correction methods need to be applied.
3. **2D Local drift correction.** This method will be applied after methods 1 and 2. It will iterate over the DAPI masks detected in the segmentation function (see below), extract a 2D region around each mask, and x-correlate the reference and cycle <i> fiducials in this 2D sub-region. Thus, this method is slower than methods 1 and 2, but provides for local corrections that account for deformations of the sample. The method will output images with the uncorrected and corrected overlaps for each DAPI mask sub-region so you can evaluate its performance. 
4. **3D local drift correction.** None of the methods above takes into account the drift of the sample in the z-plane. While this is typically very small given the good performance of autofocus, it could be an issue in some instances. This method will first apply the 2D drift obtained using methods 1 or 2 to the 3D stack of cycle <i>. Then it will background-substract and level-normalize the reference and cycle <i> fiducial images and will break them into 3D blocks (somewhat similar to method 2, which was breaking images into 2D blocks). Next, it will x-correlate every single 3D block in the reference image to the corresponding, pre-aligned block in the cycle <i> image to obtain a local 3D drift correction. The results are outputted as 3 matrices that indicate the correction applied to each block in z, x and y. In addition, a reassembled image made of XY, XZ and YZ projections is outputted to evaluate performance. Needless to say, this is the slowest but most performant method in the stack. 

##### 1- Global drift correction by cross-correlation

Global drift correction will be run by default if ```"alignByBlock": false```.

The method first finds the optimal x-y translation from fiducial images and and applies them to DAPI, barcodes and RNA images. The reference fiducial to be used needs to be indicated in ```infoList_fiducial.json``` (see above)

This can be run by typing 
```sh
runAlignImages -F .
```
in the working directory. Otherwise provide full path.

In most cases, this works very well. A good example is:

<img src="Running_pyHiM.assets/image-20200929135403059.png" alt="image-20200929135403059" style="zoom:150%;" />



##### Applying drift corrections



You can apply registrations by running:

```sh
runAppliesRegistrations.py -F .

usage: runAppliesRegistrations.py [-h] [-F ROOTFOLDER] [--parallel]
                                  [--localAlignment] [--refit]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --parallel            Runs in parallel mode
  --localAlignment      Runs localAlignment function
  --refit               Refits barcode spots using a Gaussian axial fitting
                        function.

```



##### Problems with drift correction

Exceptions are:

1. there is a bit of crap in one cycle that distorts the cross correlation. This is addressed with blockAlignment.
2. The cross correlation is biased towards a bright spot in the image that is shifted with respect to the others. This can be addressed by adapting thresholds.
3. There is a deformation in the sample that cannot be corrected by a global translation. See solution in LocalDriftCorrection (below)

##### Adapting thresholds
The function will first take the pixels whose intensities reside within the ```lower_threshold``` and the ```higher_threshold``` highest pixel intensities. This reduces problems with intense backgrounds skewing the cross-correlation.
By default, the values are ```lower_threshold=0.999``` and the ```higher_threshold=0.9999999```.

In the example below a bright spot is affecting the alignment. 

In this case this can be solved by using lower thresholds, for instance ```lower_threshold=0.99``` and the ```higher_threshold=0.999```.



##### 2- Block drift correlation

Block drift correction will be activated if ```"alignByBlock": true```.

Inspired by the use of blocks to restore large images, I implemented a new registration routine that:

1. breaks the image into non-overlapping blocks of 256x256

2. calculates the optimal translational shift between fiducial and reference in each block

3. estimates the root mean square error (RMS) between the reference and the shifted image for each block (using the optimal shift for that block). This is performed for the entire image, not the block.

4. finds the blocks with RMS within ```tolerance``` of the minimum RMS. For instance, if ```tolerance=0.1``` then a block with RMS 10% higher than that of the minimum is tolerated. 

5. Mean and standard deviation of the XY shifts are calculated from the blocks selected in step 4. Mean shifts are used for shifting the image and getting the final RMS error (reported now in the output alignment Table).

6. We re-calculate the global shift by using the entire image, and the associated RMS error. We estimate the shift errors from the polled blocks.

   1.  If the global RMS error is lower than the polled RMS error (see 4) 
   2. or if the number of blocks within the ```tolerance``` is lower than ```minNumberPollsters``` (this is not yet provided as a parameter in infoList.json, but should in future)
   3. or if the shift errors are higher than 5 pixels (this could be added as a parameter in infoList.json if needed)

   then we fall back to global alignment.

To turn the routing on, just set ```alignByBlock``` to True.

**Important note**: When you use blockAlignment, pyHiM will produce 8x8 output matrices with the abs(shift-global_shift) maps. These will be stored in ```alignImages``` folder with names such as ```scan_006_DAPI_003_ROI_converted_decon_ch02_errorAlignmentBlockMap.npy```.

When alignBarcodesMasks() runs, it will search for this files. If it finds them it will filter out barcode localizations that displayed an absolute drift larger than the ```toleranceDrift``` parameter in the ```segmentedObjects``` segment of the ```infoList_DAPI.json``` file.

So if you see a large drop in the barcodes that are used (this can be seen by matrices with empty rows/columns) it may be that your inaccurate barcode localizations are being dropped. Check for this the ```Assigned: 266 | discarded: 285``` outputted by alignBarcodesMasks(), and change ```toleranceDrift``` if needed.



###### Examples

<u>Nice ROI</u>

The right plot shows the RMS between the whole reference fiducial image and the fiducial of cycle <i> for each block. In each block, the optimal drift for that block was used for the calculation. Regions with extended blue mean that the optimal drift encountered was equally good for the whole image. Regions with high RMS indicate local block drifts not optimal for the whole image. The program will automatically select the regions with lowest RMS and use them to derive a consensus global drift.

The left plots shows the relative change between the consensus global drift and the local drift. Values are in px. White region indicates blocks where the global consensus drift is very similar to the optimal local block drift.

![scan_001_RT29_001_ROI_converted_decon_ch00_block_alignments](Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00_block_alignments.png)

This image shows the reference fiducial and the drift-corrected cycle <i> fiducial images overlapping. Regions in yellow mean good overlap.

![scan_001_RT29_001_ROI_converted_decon_ch00_overlay_corrected](Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00_overlay_corrected.png)

Left image shows the reference minus the <u>uncorrected</u> cycle <i>  fiducial images (red and blue). Blue and red regions represent places where the two images do not overlap. White regions represent regions with the same pixel values in both images. The right panel shows the reference and <u>corrected</u> cycle <i>  fiducial images (red and blue). Same colormap. So, the right image should be mostly white with almost no spot (in either blue or red) if the correction worked well. 

<img src="Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00_referenceDifference.png" alt="scan_001_RT29_001_ROI_converted_decon_ch00_referenceDifference" style="zoom:200%;" />



<u>Challenging ROI</u>

These are typical results in a challenging ROI with a bit of 'dirt' on the edge that typically skews the global cross correlation:

Example: left image is a 2D projection of a fiducial without contamination. The right panel a fiducial with a contamination in the top right.

![image-20200928145733668](Running_pyHiM.assets/image-20200928145733668.png)



Alignment using global alignment:

![scan_001_RT31_001_ROI_converted_decon_ch00_overlay_corrected](Running_pyHiM.assets/scan_001_RT31_001_ROI_converted_decon_ch00_overlay_corrected-1601972166646.png)

![scan_001_RT31_001_ROI_converted_decon_ch00_referenceDifference](Running_pyHiM.assets/scan_001_RT31_001_ROI_converted_decon_ch00_referenceDifference-1601972172822.png)

Alignment using blockAlignment:



![scan_006_DAPI_001_ROI_converted_decon_ch02_block_alignments](Running_pyHiM.assets/scan_006_DAPI_001_ROI_converted_decon_ch02_block_alignments.png)

![scan_001_RT31_001_ROI_converted_decon_ch00_overlay_corrected](Running_pyHiM.assets/scan_001_RT31_001_ROI_converted_decon_ch00_overlay_corrected.png)

![scan_001_RT31_001_ROI_converted_decon_ch00_referenceDifference](Running_pyHiM.assets/scan_001_RT31_001_ROI_converted_decon_ch00_referenceDifference.png)







#####  3- 2D LocalDriftCorrection

2D Local drift correction will be run after you run a global drift correction method either using methods 1 (global) or 2 (block alignment).  To select between these, use the ```alignByBlock``` flag.

To properly run this method, use the ```--localAlignment``` flag when you call pyHiM.py. Otherwise, run ```runLocalAlignment.py -F .```  directly from the command line, either in the working directory or providing a full path.

Deformation of samples means a simple translation will not be enough to correct drift. Typical example where most fiducial spots are corrected apart from one on the top right, very likely due to the embryo getting deformed in this cycle:

![image-20200928150007248](Running_pyHiM.assets/image-20200928150007248.png)



The ```localDriftCorrection``` function will iterate over the DAPI masks, retrieve a bounding box that is ```bezel``` pixels larger than the mask for both the reference fiducial and the fiducial of each cycle. It will then apply the same cross-correlation algorithm to find an additional local shift.  If this shift is larger than ```localShiftTolerance``` in any direction, then it will not apply it.

#####  4- 3D Local Drift Correction 

<u>How to invoke: need to find how to invoke other that runAlignImages3D.py !</u>

<u>Parameters: need to program in!</u> 

**Output of method**

The first diagnostic image shows 2D projections of the 3D reference and  cycle <i> fiducial images. Uncorrected in the top row and 3D background-subtracted and level-renormalized images. These latter will be used for the x-correlations.

![scan_001_RT29_001_ROI_converted_decon_ch00.tif_bkgSubstracted](Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00.tif_bkgSubstracted.png)



Matrices indicating the correction applied to each block in z, x and y. Values are in pixel units. You should look for roughly homogeneous corrections within embryos. Typically small local corrections are found in X and Y after if the global correction was successful. If the autofocus run fine, then the Z shifts should be ~ 1 px.

![scan_001_RT29_001_ROI_converted_decon_ch00.tif_shiftMatrices](Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00.tif_shiftMatrices-1615108552490.png)

Reassembled image made of XY, XZ and YZ projections is outputted to evaluate performance. Reference fiducial is in red, cycle <i> fiducial is in green.

![scan_001_RT29_001_ROI_converted_decon_ch00.tif_3Dalignments](Running_pyHiM.assets/scan_001_RT29_001_ROI_converted_decon_ch00.tif_3Dalignments.png)



#### Segmenting masks

To manually segment masks, run

```sh
runSegmentMasks.py -F .

usage: runSegmentMasks.py [-h] [-F ROOTFOLDER] [--parallel] [--localAlignment]
                          [--refit]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --parallel            Runs in parallel mode
  --localAlignment      Runs localAlignment function
  --refit               Refits barcode spots using a Gaussian axial fitting
                        function.

```





An example of a segmentation of a nice barcode (color indicates flux, with red being high:2000 and blue low:0, *jet colormap*):

![image-20201009113923988](Running_pyHiM.assets/image-20201009113923988.png)



An example of a barcode where localization signals are far from optimal: note that crosses are now all blue (low fluxes)

![image-20201009114117140](Running_pyHiM.assets/image-20201009114117140.png)



##### 3D fits of barcode positions using zProfiling

```"3Dmethod":"zProfile"```

Barcode 3D positions are now calculated as follows.

First ASTROPY calculates the 2D positions using 2D projected images.

Then the ```fittingSession.refitFolders``` class function draws a subVolume around each localized barcode (default: 7x7 window), estimates the axial intensity profile, substracts background, calculates the z position by calculating the weighted moment, then uses it as a seed to Gaussian fit the axial intensity distribution and find a better estimate of the z-position based on Gaussian fitting.
Finally, the localizations are filtered by following various criteria

- the residual of the fit has to be small (see parameter in infoList)
- the difference between Moment and Gaussian positions has to me small (see parameter in infoList)
- the sigma of the Gaussian fit has to be smaller than a threshold (see parameter in infoList)

The results for any given ROI and barcode appear as a figure with two subplots where these values are shown. The dots in black represent the spots that will be kept for further analysis.

![segmentedObjects_3Drefit_ROI:1_barcode:29](Running_pyHiM.assets/segmentedObjects_3Drefit_ROI1_barcode29.png)



##### 3D fits of barcode positions using ASTROPY

```"3Dmethod":"zASTROPY"```

In this case, the 2D XY source positions calculated using projected images are loaded from disk.

Then the ```fittingSession.refitFolders``` will 

- reinterpolate the 3D stack to correct for XY drift
- Make a YZ plane projection at a specific xPlane position. For this it will sum the YZ images from ```xPlane-3dAP_window:xPlane+3dAP_window``` . From this YZ projection, it will call ASTROPY and localize new sources using the parameters: 
  - ```3dAP_flux_min``` it will discard sources that have a flux smaller than this value
  - ```3dAP_brightest```: it will recover this number of sources per YZ plane
- Then, it will iterate over the sources in that YZ plane and match them to a source in the original 2D XY source list. For this, it will find the sources within a Y-distance of  ```3dAP_distTolerance```, and match to the closest source found. If no source is found closer than this threshold, it will not match it to anything. *There is an experimental setting ```addSources``` that can be used to add the sources as new sources if no match was found in the original 2D XY list of sources. By default this is set to ```False```.* 
- If the flux of the 2D XY source is smaller than that of the ZY source, the flux will be updated accordingly.
- This process is repeated for all xPlanes using ```range(3dAP_window, imageSizeX,3dAP_window)``` (i.e. the x-range of the image will be split into ```imageSizeX/3dAP_window``` equally-spaced planes). Typically ```3dAP_window=5``` works fine.

The results for any given ROI and barcode appear with three figures representing:

- a typical YZ profile in the middle of the x-range. Sources will be represented as crosses.
- the 2D-projected image with sources represented as crosses. Color-code represents z-position. Colormap is from ```0``` to the max number of z-planes in the image.
- Finally, a scatter plot displaying the z-position and the flux of every source detected is displayed. 

Important: Sources are not flux-filtered at this stage. The original or updated flux values will be outputted in the Table, and used in the ```alignBarcodesMasks``` module to filter the localization using the value of ```flux_min``` provided in the ```infoList_DAPI.json``` file.

Output examples:

![image-20201228111034499](Running_pyHiM.assets/image-20201228111034499.png)



![image-20201228111115388](Running_pyHiM.assets/image-20201228111115388.png)



![image-20201228111138959](Running_pyHiM.assets/image-20201228111138959.png)



![image-20201228111209001](Running_pyHiM.assets/image-20201228111209001.png)

#### Align DAPI masks and barcodes

This last function will align DAPI masks and barcodes and construct the single cell contact matrix.

```sh
runAlignBarcodesMasks.py -F .

usage: runAlignBarcodesMasks.py [-h] [-F ROOTFOLDER] [--parallel]
                                [--localAlignment] [--refit]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --parallel            Runs in parallel mode
  --localAlignment      Runs localAlignment function
  --refit               Refits barcode spots using a Gaussian axial fitting
                        function.

```


<u>Example outputs:</u>

Mean pairwise distance matrix. By default means are calculated using Kernel Density Estimators of the PWD distributions.

<img src="Running_pyHiM.assets/buildsPWDmatrix_HiMmatrix.png" alt="buildsPWDmatrix_HiMmatrix" style="zoom:50%;" />

In addition, the function outputs the distribution of distances for each combination of barcodes:

<img src="Running_pyHiM.assets/buildsPWDmatrix_PWDhistograms.png" alt="buildsPWDmatrix_PWDhistograms" style="zoom: 25%;" />


##### Filtering barcode localizations

There are several filters:
1. Properties of 2D localization algorithm (e.g. brightness)

2. Accuracy of 3D localization: sigma of fit, correlation between z-position from weighted moment and from gaussian fit, etc

3. Accuracy of drift correction in the region where the barcode was localized. 

   This is only applied if LocalDrift correction was **not** run. 

   

*Examples.*

| Filtering | Matrix |
| --- |  ---- |
| Unfiltered matrix. Total barcode localizations: 18700 | <img src="Running_pyHiM.assets/buildsPWDmatrix.png" alt="buildsPWDmatrix" style="zoom:25%;" /> |
|```toleranceDrift = 1px```. Barcode localizations kept: 12377 of a total: 18700.| <img src="Running_pyHiM.assets/buildsPWDmatrixFilterBlockDrift.png" alt="buildsPWDmatrixFilterBlockDrift" style="zoom:25%;" />|
| ```toleranceDrift = 1px```  ```Flux = 100```. Barcode localizations kept: 5562 of a total: 18700. | <img src="Running_pyHiM.assets/buildsPWDmatrix_filterFlux100.png" alt="buildsPWDmatrix_filterFlux100" style="zoom:25%;" /> |
|```toleranceDrift = 1px```  ```Flux = 200```. Barcode localizations kept: 4528 of a total: 18700. | <img src="Running_pyHiM.assets/buildsPWDmatrix_filterFlux.png" alt="buildsPWDmatrix_filterFlux" style="zoom:25%;" />|
|```toleranceDrift = 1px``` ```Flux = 1000```. Barcode localizations kept: 1923 of a total: 18700.| <img src="Running_pyHiM.assets/buildsPWDmatrix_filterFlux1000.png" alt="buildsPWDmatrix_filterFlux1000" style="zoom:25%;" />|



##### Outputs

In addition to the PWD matrix, we now also have a map of the alignment accuracy  and scatter plots showing the flux of each barcode, its sharpness, magnitude and roundness. These are used in order to validate the segmentation process and help with the selection of the ```flux``` threshold used in this filtering step.

*Alignment accuracy*

This provides a map of all barcode localizations in an ROI, colorcoded by the accuracy of localization. Colorbar scale is in pixels.

<img src="Running_pyHiM.assets/BarcodeAlignmentAccuracy_ROI1_2D2.png" alt="BarcodeAlignmentAccuracy_ROI:1_2D2" style="zoom: 50%;" />

*Barcode localization statistics*

This provides the localization statistics from ASTROPY. The main use of these plots is to determine if the threshold ```flux``` used is correct. Default is *200*.

<img src="Running_pyHiM.assets/BarcodeStats_ROI1_2D.png" alt="BarcodeStats_ROI:1_2D" style="zoom: 67%;" />



### Process second channel (i.e RNA, segments, etc)

```pyHiM.py``` will project all TIFFS, and align them together using the fiducial. This will include the second channel of DAPI containing RNA intensities. Now, we need to mask these files so that we can tell which cell was expressing or not a specific RNA. For this, you will run ```processSNDchannel.py```

Go to the ```destination_directory``` and run  ```processSNDchannel.py --addMask sna``` for manually segmenting all the ROIs in the destination_directory and label them with the ```sna``` tag. You can repeat this process as many times as you want using different tags. For instance, if you also want to label ```doc``` then you can run   ```processSNDchannel.py --addMask doc```. This second command will not overwrite what you did with the first command but rather accumulate different tags.

After you run this command for all the tags you want to identify, you now need to assign these tags to the cells that were previously identified during the run of ```pyHiM.py```.

For this, just run ```processSNDchannel.py``` on the command line.

#### Folder

If you don't want run ```processSNDchannel.py``` in the current directory, just choose another folder by using the ```--rootFolder``` argument.

#### erasing segmentations

If you want to start all over again and erase all manual segmentations, run with the ```--cleanAllMasks``` argument.

#### output

The output of ```processSNDchannel.py``` will be stored in ```./segmentedObjects/SNDassignedCells.ecsv``` Astrpy Table. Example:

```bash
# %ECSV 0.9
# ---
# datatype:
# - {name: 'ROI #', datatype: int64}
# - {name: 'CellID #', datatype: int64}
# - {name: 'MaskID #', datatype: string}
# schema: astropy-2.0
"ROI #" "CellID #" "MaskID #"
1 0 sna
1 12 sna
1 13 sna
1 14 sna
1 15 sna
1 16 sna
1 18 sna
1 20 sna
1 21 sna
...
```

where the first column contains the ROI number, the second the number of the cell mask, the third the tag assigned.

This file can then be loaded within ```replotHiMmatrix.py``` to identify which cells of the matrix have which tag. More on this will be added to the section below LATER.



### Analysis of several samples at once

You can now use a new script to call several samples in one go:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse

nArgs=len(sys.argv)
print("Total arguments passed: {}".format(nArgs))
EmbryoTag='Embryo_'
if nArgs > 2:
    rootDir=sys.argv[1]
    print("parameters> rootFolder: {}".format(rootDir))

    for i in range(2,nArgs):
        print("Processing Embryo #{}".format(sys.argv[i]))
        command2Run1='nice -19 pyHiM.py -F '+rootDir+EmbryoTag+sys.argv[i]
        os.system(command2Run1)
        command2Run2='zipHuMrun.py -F '+rootDir+EmbryoTag+str(i)
        os.system(command2Run2)
        print('Commands: {}\n{}'.format(command2Run1,command2Run2))
else:
    print('not enough arguments.')
```



To run just do:

```bash
pyHiM.py -F rootFolder 0 1 33
```

to run Emrbyo_0, Embryo_1 and Embryo_33 from rootFolder



***I will change this in the future to replace it by either Dask or snakeMake...***



### Parallel Computations

Several routines are now fitted with the possibility of performing parallel computations using the Dask package.

The use of parallel computations can be invoked using the ```--parallel``` flag in the command line or for any individual functions (e.g. ```runMakeProjections```) or for ```pyHiM.py```.

The gain is roughly two fold currently, but it can be up to 10 fold for some functions.

To monitor parralel computations, SSH into the cluster using

```sh
ssh -L 8787:localhost:8787 marcnol@lopevi
```

Invoke using the ```--parallel``` flag in ```pyHiM.py```.

Go to your browser and open ```http://localhost:8787``` to see the progress.



### zipping and erasing run

#### zip and retrieve results

Other utilities have been written to retrieve data from a run to a remote server. For this, go to the directory with the data an run ```zipHiM_run.py```. This will archive all the png, ecsv, and dat files with the results in addition to the MD file with the output of the run. The name of the archive will be HiMrun.tar.gz.

#### clean run

If you want to erase a run, for instance to make sure you can run it again without any leftover, you can run ```cleanHiM_run.py` in the directory with the data.

#### Link folder

If you want to link many files in a folder without having to copy very large  files (e.g. TIFF) then you can run  ```lndir```.

This script  links of files in a directory to a second directory. it selects the files in the first folder by following rules provided by the user where wildcards are possible.

For instance, In the command line, 

```sh
$ lndir.py "/home/marcnol/Repositories/pyHiM/*py" ~/Downloads/test
```
to **link** all files with extension ```py``` from ```/home/marcnol/Repositories/pyHiM/``` to ```~/Downloads/test```.

Make sure that the first argument has quotation marks if you use wildcards!



## 2- Combining results from different experiments

Once you run a bunch of datasets, you will want to combine the PWD matrices together. For this:

1. Retrieve the matrices and barcodes files by scp:

```bash
scp rata@lopevi:/mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/buildsPWDmatrix/*ecsv /home/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix/

scp rata@lopevi:/mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/buildsPWDmatrix/*npy /home/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix/
```

Now, you can run ```processHiMmatrix.py``` locally. You should setup your files in a directory. For instance the directory ```/mnt/disk2/marcnol/data/Experiment_19``` contains three folders:

```bash
006_Embryo  009_Embryo  026_Embryo
```

containing each of the analysis from different embryos of the same experiment. Now, in this directory, you should create a file called ```folders2Load.json``` with the following:

```python
{
    "wt_docTAD": {
        "Folders": [
            "/mnt/disk2/marcnol/data/Experiment_19/026_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_19/009_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_19/006_Embryo/buildsPWDmatrix"
        ],
        "PWD_clim": 1.4,
        "PWD_mode": "median",
        "PWD_cm": "terrain",
        "iPWD_clim": 6,
        "iPWD_mode": "median",
        "iPWD_cm": "terrain",
        "ContactProbability_scale": 12,
        "ContactProbability_cmin": 0.0,
        "ContactProbability_distanceThreshold": 0.35,
		"ContactProbability_cm": "coolwarm",
		"BarcodeColormap": [4, 4, 4, 4, 8, 4, 3, 4, 4, 8, 3, 4, 4, 8, 4, 8, 3],
		"3wayContacts_anchors": [7, 11, 17, 5, 10, 14]
    }
}
```

```

This file contains the directories with the data to be analyzed and some parameters for the analysis (more on this later).

Go to the root folder (```/mnt/disk2/marcnol/data/Experiment_19```) and run

​```bash
 processHiMmatrix.py
```

If you use a parameter file with another name (e.g. ```myparameters.json```) then run:

```bash
 processHiMmatrix.py --parameters myparameters.json
```

This will produce an MD file (e.g. ```processHiMmatrixAnalysis__wt_docTAD_27052020_140943.md```) with the following output:

- PWD matrix for each dataset
- Inverse distance matrices for each dataset
- Contact probability matrices for each dataset
- Combined contact probability matrix

This last matrix will be outputed in the ```scHiMmatrices``` directory as two files. Examples:

- ```CombinedMatrixwt_docTAD.dat```: plain text ensemble HiM contact probability matrix (can be opened in MATLAB). Each row in the matrix is separated by a ```\n``` .
-  ```UniqueBarcodeswt_docTAD.dat```: plain text file with the barcodes used.

### Analyzing labeled datasets

If you run ```processSNDchannel.py``` before, you may want now to look at cells with different labels (ON, OFF, etc). For this, you need to run ```processHiMmatrix.py``` with two more parameters:
- ```--label```: indicates the name used when you run processSNDchannel.py with the option ```--addMask```. Typical names: doc, sna.Running
- ```--action```: three options are available:
  - ```all```: selects all cells for analysis irrespective of whether they are labeled
  - ```labeled``` only runs analysis on labeled cells
  - ```unlabeled```: only runs analysis on unlabeled cells

### options in parameter file

Options:

- ```PWD_clim```: 1.4. Maximum of colormap for PWD matrices.
- ```PWD_mode```: Mode used to calculate value for each mean from many measurements. ```median``` is the median excluding NaNs, ```KDE``` uses kernel density estimator and peaks up the maximum (uses 0.2 as size of kernel as this works for most situations).
- ```iPWD_clim```: 6. Same as for ```PWD_clim.```
- ```iPWD_mode```: Same as for ```PWD_mode```.
- ```ContactProbability_scale```:  normalization factor for the contact probability. Plotted contact probability is calculated as the ratio of the calculated contact probability and ```ContactProbability_scale```.
- ```ContactProbability_cmin```: Minimum of colormap in the the contact probability map.
- ```ContactProbability_distanceThreshold```: distance used for the calculation of the contact probabilities in pixel units.


### Analyzing MATLAB datasets

One can use this same pipeline to analyze data that was segmented using MATLAB routines (i.e merfish_main.m). For this, a few steps need to be taken.

#### Make directory structure

Make a directory structure where you will load your data to be read by processHiMmatrix.py. Follow this example closely:

```bash
.
├── 000_Embryo
│   └── buildsPWDmatrix
│       ├── buildsPWDmatrix_uniqueBarcodes.ecsv
│       └── HiMscMatrix.mat
```

place a MAT file with the name ```HiMscMatrix.mat``` in the ```buildsPWDmatrix``` directory. This file should contain the ```distanceMatrixCumulative``` variable, a 3-D matrix with the format barcode x barcode x nCells. Make sure you remove empty barcode rows by saving ```distanceMatrixCumulative(p.listofRTsusedGenomicallySorted,p.listofRTsusedGenomicallySorted,:)```

Then, you need to have a file called ```buildsPWDmatrix_uniqueBarcodes.ecsv``` with the barcodes in a line separated list, for example:

```bash
3
4
5
6
7
13
23
25
29
37
43
48
60
62
66
67
68
89
90
91
92
```

#### Create folders2Load.json file

Now you create the parameters file. An example follows:

```bash
{
    "wt_Pc_Chr3R": {
        "Folders": [
            "/home/marcnol/data/Experiment_Julian/000_Embryo/buildsPWDmatrix"
        ],
        "PWD_clim": 1.4,
        "PWD_mode": "median",
        "PWD_cm": "terrain",
        "iPWD_clim": 6,
        "iPWD_mode": "median",
        "iPWD_cm": "terrain",
        "ContactProbability_scale": 12,
        "ContactProbability_cmin": 0.0,
        "ContactProbability_distanceThreshold": 0.35,
	"ContactProbability_cm": "coolwarm",
	"BarcodeColormap": [4, 4, 4, 4, 8, 4, 3, 4, 4, 8, 3, 4, 4, 8, 4, 8, 3],
	"3wayContacts_anchors": [7, 11, 17, 5, 10, 14]
    }
}
```

Make sure you edit the name of the dataset (here ```wt_Pc_Chr3R```) and add a folder for each dataset you want to analyze.

#### Run processHiMmatrix.py

You can now run the script. For instance, do

```bash
processHiMmatrix.py  --matlab
```

to run with default options. The important thing is to add the ```--matlab``` flag.

You should be now set.

## 3- Plotting publication-quality figures

The processing with ```processHiMmatrix.py``` produces plots but this routing is primarily concerned with the collection of different datasets, thus it does not have options for customization.

This is performed in a series of ***py*** routines  that plot publication-quality figures with many options.

#### Plotting single HiM matrices

```figureHiMmatrix.py``` will plot a single matrix. It has many options, described when you run it with the --help parameter:

```bash
usage: figureHiMmatrix.py [-h] [-F ROOTFOLDER] [-O OUTPUTFOLDER]
                          [-P PARAMETERS] [-A LABEL] [-W ACTION]
                          [--fontsize FONTSIZE] [--axisLabel] [--axisTicks]
                          [--barcodes] [--scalingParameter SCALINGPARAMETER]
                          [--plottingFileExtension PLOTTINGFILEEXTENSION]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with dataset
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder for outputs
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A LABEL, --label LABEL
                        Add name of label (e.g. doc)
  -W ACTION, --action ACTION
                        Select: [all], [labeled] or [unlabeled] cells plotted
  --fontsize FONTSIZE   Size of fonts to be used in matrix
  --axisLabel           Use if you want a label in x and y
  --axisTicks           Use if you want axes ticks
  --barcodes            Use if you want barcode images to be displayed
  --scalingParameter SCALINGPARAMETER
                        Scaling parameter of colormap
  --plottingFileExtension PLOTTINGFILEEXTENSION
                        By default: svg. Other options: pdf, png
  --shuffle SHUFFLE     Provide shuffle vector: 0,1,2,3... of the same size or
                        smaller than the original matrix. No spaces! comma-
                        separated!    
  --scalogram           Use if you want scalogram image to be displayed

```



Most of these options are self-explanatory. Here is the command to run to get Figure 1F made:

```bash
figureHiMmatrix.py -F "$DATA1" --fontsize 22 --label doc --action labeled --scalingParameter 1 --barcodes --outputFolder "$FIGS"/Figure1  --plottingFileExtension png
```

And here is the output:

![Fig_HiMmatrix_dataset1:wt_docTAD_nc14_label:doc_action:labeled](Running_pyHiM.assets/Fig_HiMmatrix_dataset1wt_docTAD_nc14_labeldoc_actionlabeled.png)

##### Altering order or number of barcodes in displayed matrix

If you want to remove some of the barcodes or change their order, then you need to use the --shuffle option. For this, you need to provide with the order of the new barcodes. For instance if you want to plot only the first 9 barcodes then run:

```bash
figureHiMmatrix.py -F "$DATA1" --fontsize 22 --label doc --action labeled --scalingParameter 1 --barcodes --outputFolder "$FIGS"/Figure1  --plottingFileExtension png --shuffle 0,1,2,3,4,5,6,7,8
```

If you want to change the order, for instance, put the last bin (bin 16) first. The run:

```bash
figureHiMmatrix.py -F "$DATA1" --fontsize 22 --label doc --action labeled --scalingParameter 1 --barcodes --outputFolder "$FIGS"/Figure1  --plottingFileExtension png --shuffle 16,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
```

Be careful with the format of the ```--shuffle``` parameter. It should contain **no spaces**, and the indexes should be **comma-separated integers**. If the number of integers in the vector exceeds the number of dimensions of the matrix, the output plot will appear blue (matrix will contain zeros). If any given index exceeds the matrix dimensions, then it will be ignored but the matrix substitution will be wrong!

##### Scalograms

You can now produce a scalogram from the matrix. For this, use the --scalogram option, which will output a separate output image with the scalogram. As example:

![Fig_HiMmatrix_scalogram_dataset1:wt_docTAD_nc14_label:doc_action:labeled](Running_pyHiM.assets/Fig_HiMmatrix_scalogram_dataset1wt_docTAD_nc14_labeldoc_actionlabeled.png)

#### Plotting 3-way contact matrices

figure3wayInteractions.py``` will plot 6 3-way contact matrices matrix. The anchors used are those provided in the folders2Load.json file of the analysis. This function has many options, described when you run it with the --help parameter:

```bash
usage: figure3wayInteractions.py [-h] [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2]
                                 [-O OUTPUTFOLDER] [-P PARAMETERS]
                                 [-A1 LABEL1] [-W1 ACTION1] [-A2 LABEL2]
                                 [-W2 ACTION2] [--fontsize FONTSIZE]
                                 [--scalingParameter SCALINGPARAMETER]
                                 [--colorbar]
                                 [--plottingFileExtension PLOTTINGFILEEXTENSION]

optional arguments:
  -h, --help            show this help message and exit
  -F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
                        Folder with dataset 1
  -F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
                        Folder with dataset 2
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder for outputs
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A1 LABEL1, --label1 LABEL1
                        Add name of label for dataset 1 (e.g. doc)
  -W1 ACTION1, --action1 ACTION1
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  -A2 LABEL2, --label2 LABEL2
                        Add name of label for dataset 1 (e.g. doc)
  -W2 ACTION2, --action2 ACTION2
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  --fontsize FONTSIZE   Size of fonts to be used in matrix
  --scalingParameter SCALINGPARAMETER
                        Scaling parameter of colormap
  --colorbar            Use if you want a colorbar
  --plottingFileExtension PLOTTINGFILEEXTENSION
                        By default: svg. Other options: pdf, png

```



Run to produce Figure 1G:

```bash
figure3wayInteractions.py -F1 "$DATA1" --label1 doc --action1 labeled --fontsize 12 --scalingParameter 1.0 --outputFolder "$FIGS"/Figure1  --plottingFileExtension png
```



which produces:

![Fig_3wayContacts_dataset1:wt_docTAD_nc14_label1:doc_action1:labeled](Running_pyHiM.assets/Fig_3wayContacts_dataset1wt_docTAD_nc14_label1doc_action1labeled.png)



It is possible also to compare the 3-way matrices of two datasets, for this, run:

```bash
figure3wayInteractions.py -F1 "$DATA2" --label1 doc --action1 labeled -F2 "$DATA2" --label2 NE --action2 labeled --fontsize 22 --scalingParameter 1.0  --outputFolder "$FIGS"/Figure2 --plottingFileExtension png
```

to obtain:

![Fig_3wayContacts_dataset1:wt_docTAD_nc14_label1:doc_action1:labeled_dataset2:wt_docTAD_nc14_label2:NE_action2:labeled](Running_pyHiM.assets/Fig_3wayContacts_dataset1wt_docTAD_nc14_label1doc_action1labeled_dataset2wt_docTAD_nc14_label2NE_action2labeled.png)

The first quadrant (top right) will represent the first dataset.



#### Compare two HiM matrices

```figureCompare2Matrices.py``` will make two panels to compare datasets. One will contain the HiM contact matrices of each datasets in different quadrants. The first dataset will be plotted on top right. The second plot will produce either a difference matrix (default) or the log(ratio) (see --ratio option). The options are described when you run it with the --help parameter:

```bash
usage: figureCompare2Matrices.py [-h] [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2]
                                 [-O OUTPUTFOLDER] [-P PARAMETERS]
                                 [-A1 LABEL1] [-W1 ACTION1] [-A2 LABEL2]
                                 [-W2 ACTION2] [--fontsize FONTSIZE]
                                 [--axisLabel] [--axisTicks] [--ratio]
                                 [--cAxis CAXIS]
                                 [--plottingFileExtension PLOTTINGFILEEXTENSION]

optional arguments:
  -h, --help            show this help message and exit
  -F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
                        Folder with dataset 1
  -F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
                        Folder with dataset 2
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder for outputs
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A1 LABEL1, --label1 LABEL1
                        Add name of label for dataset 1 (e.g. doc)
  -W1 ACTION1, --action1 ACTION1
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  -A2 LABEL2, --label2 LABEL2
                        Add name of label for dataset 1 (e.g. doc)
  -W2 ACTION2, --action2 ACTION2
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  --fontsize FONTSIZE   Size of fonts to be used in matrix
  --axisLabel           Use if you want a label in x and y
  --axisTicks           Use if you want axes ticks
  --ratio               Does ratio between matrices. Default: difference
  --cAxis CAXIS         absolute cAxis value for colormap
  --plottingFileExtension PLOTTINGFILEEXTENSION
                        By default: svg. Other options: pdf, png

```

To produce Figure 2B you run

```bash
figureCompare2Matrices.py -F1 "$DATA2" -F2 "$DATA2" --label1 doc --label2 M --action1 labeled --action2 labeled --cAxis 2 --fontsize 22 --outputFolder "$FIGS"/Figure2 --plottingFileExtension png --ratio
```

and will obtain:

![Fig_mixedHiMmatrices_dataset1:wt_docTAD_nc14_label1:doc_action1:labeled_dataset2:wt_docTAD_nc14_label2:M_action2:labeled](Running_pyHiM.assets/Fig_mixedHiMmatrices_dataset1wt_docTAD_nc14_label1doc_action1labeled_dataset2wt_docTAD_nc14_label2M_action2labeled.png)

and:

<img src="Running_pyHiM.assets/Fig_ratio2HiMmatrices_dataset1wt_docTAD_nc14_label1doc_action1labeled_dataset2wt_docTAD_nc14_label2M_action2labeled.png" alt="Fig_ratio2HiMmatrices_dataset1:wt_docTAD_nc14_label1:doc_action1:labeled_dataset2:wt_docTAD_nc14_label2:M_action2:labeled" style="zoom: 80%;" />



#### Plot 4M profiles

```figure4Mmatrix.py``` will plot 4M profiles using the anchors defined in folders2Load.json parameters file. The options are described when you run it with the --help parameter:

```bash
usage: figure4Mmatrix.py [-h] [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2]
                         [-O OUTPUTFOLDER] [-P PARAMETERS] [-A1 LABEL1]
                         [-W1 ACTION1] [-A2 LABEL2] [-W2 ACTION2]
                         [--fontsize FONTSIZE] [--axisLabel] [--axisTicks]
                         [--splines] [--cAxis CAXIS]
                         [--plottingFileExtension PLOTTINGFILEEXTENSION]

optional arguments:
  -h, --help            show this help message and exit
  -F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
                        Folder with dataset 1
  -F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
                        Folder with dataset 2
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder for outputs
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A1 LABEL1, --label1 LABEL1
                        Add name of label for dataset 1 (e.g. doc)
  -W1 ACTION1, --action1 ACTION1
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  -A2 LABEL2, --label2 LABEL2
                        Add name of label for dataset 1 (e.g. doc)
  -W2 ACTION2, --action2 ACTION2
                        Select: [all], [labeled] or [unlabeled] cells plotted
                        for dataset 1
  --fontsize FONTSIZE   Size of fonts to be used in matrix
  --axisLabel           Use if you want a label in x and y
  --axisTicks           Use if you want axes ticks
  --splines             Use if you want plot data using spline interpolations
  --cAxis CAXIS         absolute cAxis value for colormap
  --plottingFileExtension PLOTTINGFILEEXTENSION
                        By default: svg. Other options: pdf, png
```



To produce Figure S1H, one would run:

```bash
figure4Mmatrix.py  --rootFolder1 "$DATA1" --label1 doc --action1 all --cAxis 1 --outputFolder "$FIGS"/Figure1  --plottingFileExtension png
```

to obtain:

![Fig_4Mcontacts_dataset1:wt_docTAD_nc14_label1:doc_action1:all](Running_pyHiM.assets/Fig_4Mcontacts_dataset1wt_docTAD_nc14_label1doc_action1all.png)

One can also plot multiple datasets, for example:

```bash
figure4Mmatrix.py --rootFolder1 "$DATA2" --rootFolder2 "$DATA2" --label1 doc --label2 M --action1 labeled --action2 labeled --cAxis 1 --outputFolder "$FIGS"/Figure2 --plottingFileExtension png
```

produces

![Fig_4Mcontacts_dataset1:wt_docTAD_nc14_label1:doc_action1:labeled_dataset2:wt_docTAD_nc14_label2:M_action2:labeled](Running_pyHiM.assets/Fig_4Mcontacts_dataset1wt_docTAD_nc14_label1doc_action1labeled_dataset2wt_docTAD_nc14_label2M_action2labeled.png)

#### Plotting different matrices together

This function is useful when you want to plot side by side the matrices from different datasets. One typical example is the matrices for different segments of an embryo.

```sh
figureN_HiMmatrices.py

usage: figureN_HiMmatrices.py [-h] [-F ROOTFOLDER] [-O OUTPUTFOLDER]
                              [-P PARAMETERS] [-A LABEL] [-W ACTION]
                              [--fontsize FONTSIZE] [--axisLabel]
                              [--axisTicks] [--barcodes]
                              [--scalingParameter SCALINGPARAMETER]
                              [--plottingFileExtension PLOTTINGFILEEXTENSION]
                              [--shuffle SHUFFLE] [--scalogram] [--type TYPE]
                              [--pixelSize PIXELSIZE] [--cAxis CAXIS]
                              [--ratio] [--normalizeMatrix]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with dataset
  -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder for outputs
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A LABEL, --label LABEL
                        Add name of label (e.g. doc)
  -W ACTION, --action ACTION
                        Select: [all], [labeled] or [unlabeled] cells plotted
  --fontsize FONTSIZE   Size of fonts to be used in matrix
  --axisLabel           Use if you want a label in x and y
  --axisTicks           Use if you want axes ticks
  --barcodes            Use if you want barcode images to be displayed
  --scalingParameter SCALINGPARAMETER
                        Scaling parameter of colormap
  --plottingFileExtension PLOTTINGFILEEXTENSION
                        By default: svg. Other options: pdf, png
  --shuffle SHUFFLE     Provide shuffle vector: 0,1,2,3... of the same size or
                        smaller than the original matrix. No spaces! comma-
                        separated!
  --scalogram           Use if you want scalogram image to be displayed
  --type TYPE           Provide one of the following: PWD, contact, iPWD
  --pixelSize PIXELSIZE
                        Provide pixelSize in um
  --cAxis CAXIS         absolute cAxis value for colormap
  --ratio               Does ratio between matrices. Default: difference
  --normalizeMatrix     Normalizes matrices by maximum. Default: True

```
The two outputs are:

1. the entire matrices 

![image-20200928162116627](Running_pyHiM.assets/image-20200928162116627.png)



2. A subset of bins, for instance to highlight how a TAD changes (normalized by the first matrix):

![image-20200928162310354](Running_pyHiM.assets/image-20200928162310354.png)

3. A line scan over all the matrices normalized by a segment

![image-20200928162417318](Running_pyHiM.assets/image-20200928162417318.png)



#### Displaying single cells

A thorough description of the work in displaying single cells will be maintained in the following MD file within the */doc* directory of pyHiM:

[Single Cell MD tutorial](SingleCellsAnalysis.md)

