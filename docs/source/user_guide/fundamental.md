# pyHiM fundamentals

## Pipeline overview
### Default pyHiM flow

To run default pipeline, pyHiM need two kind of data:
- A dictionary of initialization parameters as a JSON file, name `infoList.json`
- 3D images with TIFF format, 4 types of images are accepted and they are processed in this order:
	1. Fiducial
	2. Barcode
	3. Mask (like DAPI)
	4. RNA

These types of images are called **labels**.

Steps of pipeline are executed sequentially in this order:

1. **makeProjections**: Projects all 3D images in 2D 
2. **alignImages**: Compute the best shift to align all 2D fiducials.
3. **appliesRegistrations**: Shifts 2D barcodes, masks and RNA with result of alignImages step.
4. **alignImages3D**: Takes 2D aligned fiducials and find the best shift on Z-axis. This shift will be apply on the 3D segmented barecodes at buildHiMmatrix step.
5. **segmentMasks**: Segments 2D aligned barcodes and masks with two different ways for each.
6. **segmentSources3D**: Applies 2D shift, computed at alignImages step, to 3D barcodes. Then, segments them in 3D.
7. **buildHiMmatrix**: Filters the segmentation results, associates barcode coordinates with the good mask and makes the PWD matrix for each mask.

But there are specific routines for specific labels, so this is the real running order:

|Command|Fiducial|Barcode|DAPI|RNA|
|:-:|:-:|:-:|:-:|:-:|
|makeProjections|1|4|8|12|
|alignImages|2||||
|appliesRegistrations||5|9|13|
|alignImages3D|3||||
|segmentMasks||6|10||
|segmentSources3D||7|||
|buildsPWDmatrix|||11||

### Input / Output data

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

### Processing flowchart for different data sources

Different data sources exist and need to processed differently by pyHiM:
- fiducials
- DNA-FISH barcodes
- masks (e.g. DAPI, cell, genomic region)
- RNA-FISH images

Each data source is handled by a different workflow in pyHiM. In each workflow, you will see different symbols for `features`, `input data` and `I/O data`, as follows:

```{mermaid}
flowchart TD
	subgraph graph legend
		Z1((feature))
		Z2([initial data])
		Z3[I/O data]
	end 
```

#### fiducial flow
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

#### barcode flow
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
#### masks flow
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
#### RNA-FISH image flow
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

## Main features
### makeProjections
*Projects 3D images in 2D*
#### Invoke
To run this function exclusively, run *pyHiM* using the ``` -C makeProjections ``` argument. This routine take all 3D images and project its in 2D. Depending on the chosen *mode*, this feature start to find the good set of Z-plans, where there is the least noise. This step give a range centered on a focal plan, named *zRange*. After, projection is done on this range either by sum or by maximum intensity projection.

#### Relevant options
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

### alignImages
*Registers fiducials using a barcode as reference*

#### Invoke

To run this function exclusively, run *pyHiM* using the ``` -C alignImages ``` argument. 
In the set of *fiducial* images, one is chosen by initialization parameters to be the reference. 
The algorithm takes images one by one and align with the reference.
There are two ways to compute the shift:
- Global alignement make simple cross-correlation with tow images
- Split image in block and make cross-correlation block by block. Then we have one shift by block and to align the global image an average of those shifts are made. This method is more robust against a bright noise spot.

#### Relevant options
Parameters for this script will be read from the  ```alignImages``` field of ```infoList.json```

|Name|Option|Description|
|:-:|:-:|:-:|
|referenceFiducial| |Selects reference barcode image|

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

### AppliesRegistrations
*Applies registration to DAPI and barcodes*

#### Invoke
To run this function exclusively, run *pyHiM* using the ``` -C appliesRegistrations ``` argument. It loads masks, RNA, and barcodes 2D projected images, and applies registrations to them. The resulting images are saved as npy arrays in the ```alignImages``` folder. 

```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([2D_DAPI])
			A2([2D_RNA])
			A3([2D_Barcodes])
			A4([dictShifts])
		end
		subgraph C[OUTPUT]
			C1([aligned_2D_DAPI])
			C2([aligned_2D_RNA])
			C3([aligned_2D_Barcodes])
		end
		A --> B1
		B1[["appliesRegistrations2fileName()"]] --> C
	
```

### alignImages3D
*Aligns fiducials in 3D*

#### Invoke
To run this function exclusively, run *pyHiM* using the ``` -C alignImages3D ``` argument.


#### Relevant options
Parameters for this script will be read from the  ```alignImages``` field of ```infoList.json```. 
To run, the value for ```localAlignment``` key should be ```block3D```. 

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


### segmentMasks
*Segments DAPI and sources in 2D*

#### Invoke
To run this function exclusively, run *pyHiM* using the ``` -C segmentMasks ``` argument.

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


### segmentSources3D
*Segments sources in 3D*

#### Invoke
To run this function exclusively, run *pyHiM* using the ``` -C segmentSources3D ``` argument.


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

## Building chromatin traces

### Build traces: new method

The new method requires executing several modules:

- `filter_localizations`
- `register_localizations`
- `build_traces`
- `build_matrices`



#### filter_localizations



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



#### register_localizations



**Invoke**

To run this function exclusively, run *pyHiM* using the ```-C register_localizations``` argument. This function will find and process all the localization files in the `segmentedObjects` folder. To avoid overwriting data, existing files will be renamed with the extension `_version_n` where `n`will be incremented from run to run. The output of `register_localizations` will be saved with the original localizations filename. A comment in the header will be added to indicate that a *registration* operation was run on this file. `register_localizations` **will not be run on files that were previously registered.**



**Relevant options**

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

```
"toleranceDrift": 1 # tolerance drift in pixels. Above this value localizations will not be locally registered
```



outputs images:

- `_registered_barcode_stats.png`
- `_registered_barcode_localizations_ROI*.png`

| statistics of registration                                   | localization map                                             |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20220210221852444](../_static/user_guide/image-20220210221852444.png) | ![image-20220210221942291](../_static/user_guide/image-20220210221942291.png) |
| ![image-20220210222028835](../_static/user_guide/image-20220210222028835.png) | ![image-20220210222006297](../_static/user_guide/image-20220210222006297.png) |



#### build_traces



**Invoke**

To run this function exclusively, run *pyHiM* using the ```-C build_traces``` argument. This function will find and process all the localization files in the `segmentedObjects` folder. The output of `register_localizations` will be saved in the `buildsPWDmatrix` folder with the name starting with `Trace_`. The reminder of the name will contain the kind of operation run (mask/KDtree) the identity of the mask (e.g. mask0, DAPI), and whether localizations used were from a *2D* or a *3D* analysis. 



**Relevant options**

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

```
"tracing_method": ["masking","clustering"], # list of methods it will use
"mask_expansion": 8,# number of pixels masks will be expanded to assign localizations
"masks2process":{"nuclei":"DAPI","mask1":"mask0"}, # masks identities to process
"KDtree_distance_threshold_mum": 1,# threshold distance used for KDtree clustering
```



Output images:

- `_XYZ_ROI*.png`

|  | full image | zoomed images |
| --- |   ---- | --- |
| 3D **mask** | ![image-20220210221402082](../_static/user_guide/image-20220210221402082.png) |![image-20220210221430543](../_static/user_guide/image-20220210221430543.png)|
| 3D **mask** | ![image-20220210222233148](../_static/user_guide/image-20220210222233148.png) |![image-20220210222354093](../_static/user_guide/image-20220210222354093.png)|
| 3D **KDtree** |  ||



#### trace_selector

The objective of this script is to attribute labels (e.g. gene expression, morphological, etc) to single traces. 

**Invoke**

To run this function, run ```trace_selector.py``` at the command line. This function will find and process all the `Trace` files within the `buildsPWDmatrix` folder with the name starting with `Trace_`. To provide a specific list of `Trace` files, please see section below (Piping file lists using the `--pipe`  option). Each trace file will be edited so that the `labels` colums now contains a comma-separated list of the labels where it was found. It will then search and find the numpy (NPY) files within `segmentedObjects` with the keyword `SNDmask` in the filename. These files can be generated by running `processSNDchannel.py`.


**Relevant options**

```
usage: trace_selector.py [-h] [-F ROOTFOLDER] [--pixel_size PIXEL_SIZE] [--pipe]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --pixel_size PIXEL_SIZE
                        Lateral pixel size un microns. Default = 0.1
  --pipe                inputs Trace file list from stdin (pipe)
```

As coordinates in the `Trace` file are in micrometers, `trace_selector.py` needs to know how to convert them to pixel units to match `SNDchannel` masks.

**Piping file lists using the `--pipe`  option** 

When this argument is not invoked, `trace_selector` will read all `Trace_` files within `BuildsPWDmatrix`.

If this argument is provided, you can send a file list of `Trace_` files to process using piping. For instance:

```sh
ls buildsPWDmatrix/Trace_*3D*2.ecsv | trace_selector.py --pipe
```

will list a selection of Trace files and send it to `trace_selector` to process.



#### trace_combinator

The objective of trace_combinator is to combine traces from different ROIs or different experiments.

**Invoke**

`trace_combinator` can be run from the command line in two manners: 

- Provide input parameters in a `folders2Load.json` dictionary (see example below) at execution folder or at a `rootFolder` provided as argument. The script will read ALL the `Trace` files in the paths  provided in `folders2Load.json` and will combine them into a unique  `Trace_` files. 
  The output of `trace_combinator` will be generated at the `rootFolder` in the directory `combined_traces` that will contain a new `Trace` file with the combined traces.

  <u>*warning*</u>: `trace_combinator` will combine **ALL** traces found in the `buildPWDmatrix` folders within `folders2load.json`. If you want to filter which *Trace* files are used, we recommend you to use the second method described below. Alternatively, you can use the `--method` command-line argument to filter which `Trace` files are used. In that case, make sure to verify the terminal output so ensure the right *Trace* files have been used.

- Al alternative way of using `trace_combinator` is by using pipes. In this case, a list of `Trace` files to process is piped into `trace_combinator` and the `--pipe` argument is invoked. Piping can be either done using a command or by sending a list in a file. See examples below.

  

  Example 1: using `ls`

  ```sh
  # use ls to select which files you want to combine
  $ ls folder1/Trace*3D*6.ecsv
  folder1/Trace_3D_barcode_KDtree_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:DAPI_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:mask0_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ ls folder1/Trace*3D*6.ecsv | trace_combinator.py --pipe
  ```

  this will process the three `Trace` files listed using `ls`.

  

  Example 2: using cat

  ```sh
  # first make list of files to process and write it in a file
  $ cat files_to_combine 
  folder1/Trace_3D_barcode_KDtree_ROI:15.ecsv
  folder2/Trace_3D_barcode_KDtree_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ cat files_to_combine | trace_combinator.py --pipe
  ```

  this will process the files within `files_to_combine`.



**Use `trace_combinator` to select traces with specific *labels***

`trace_combinator` can also be used to combine traces with specific, user-provided  `labels`, which can be attributed using `trace_selector` (see above). These can be indicated by using the `--label` command-line argument. In addition, you need to indicate whether you want to use only the traces that contain the label (`labeled`), the traces that <u>do not</u> contain the label (`unlabeled`), or all the traces (`all`) irrespective of whether they contain or not the label. It is good practice to verify the `label` column of your output `Trace_` file to check that you selected the expected traces.



**Relevant options**

```
usage: trace_combinator.py [-h] [-F ROOTFOLDER] [-P PARAMETERS] [-A LABEL]
                           [-W ACTION] [--saveMatrix] [--ndims NDIMS]
                           [--method METHOD] [--pipe]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  -P PARAMETERS, --parameters PARAMETERS
                        Provide name of parameter files. folders2Load.json
                        assumed as default
  -A LABEL, --label LABEL
                        Add name of label (e.g. doc)
  -W ACTION, --action ACTION
                        Select: [all], [labeled] or [unlabeled] cells plotted
  --saveMatrix          Use to load matlab formatted data
  --ndims NDIMS         Dimensions of trace
  --method METHOD       Method or mask ID used for tracing: KDtree, mask, DAPI
  --pipe                inputs Trace file list from stdin (pipe)
```



**folders2Load.json template**

```sh
{
    "dataset_name": {
        "Folders": [
            "/mnt/disk2/marcnol/data/Experiment_4/0_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_4/1_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_4/2_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_4/4_Embryo/buildsPWDmatrix",
            "/mnt/disk2/marcnol/data/Experiment_4/5_Embryo/buildsPWDmatrix"
        ],
        "PWD_clim": 1.4,
        "PWD_mode": "median",
        "iPWD_clim": 6,
        "iPWD_mode": "median",
        "ContactProbability_scale": 15,
        "ContactProbability_cmin": 0.0,
        "ContactProbability_distanceThreshold": 0.25
    }
}
```



#### trace_filter

This script is intended to filter the traces with a low number of spots. It acts on a single or several trace files and takes as input the minimum number of spots a trace needs to have to be kept.

How to use:

- Example 1: Invoke from `rootFolder`

  ```sh
  # Call trace_filter from rootFolder to process all Trace files in `buildPWDmatrix`
  $ trace_filter.py -- N_barcodes 2
  ```

  this will process all  `Trace` files in `buildPWDmatrix`

  

- Example 2: using cat or ls to provide a specific list of Trace files to process

  ```sh
  # either make list of files to process and write it in a file
  $ cat files_to_combine 
  folder1/Trace_3D_barcode_KDtree_ROI:15.ecsv
  folder2/Trace_3D_barcode_KDtree_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ cat files_to_combine | trace_combinator.py --N_barcodes 2 --pipe
  
  # OR use `ls` to select which files you want to combine
  $ ls folder1/Trace*3D*6.ecsv
  folder1/Trace_3D_barcode_KDtree_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:DAPI_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:mask0_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ ls folder1/Trace*3D*6.ecsv | trace_combinator.py --N_barcodes 2 --pipe
  ```




**Relevant options**

```sh
usage: trace_filter [-h] [-F ROOTFOLDER] [--N_barcodes N_BARCODES] [--pipe]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --N_barcodes N_BARCODES
                        minimum_number_barcodes. Default = 2
  --pipe                inputs Trace file list from stdin (pipe)
```



#### build_matrices

This script will build single-cell pair-wise distance (PWD) matrices, proximity frequency maps, and N-matrices from each `Trace_` file in the `buildsPWDmatrix` folder. 

**Invoke**

To run this function exclusively, run *pyHiM* using the ```-C build_matrix``` argument. This function will find and process all the `Trace_` files in the `buildsPWDmatrix` folder. The outputs of `build_matrix` will be saved in the `buildsPWDmatrix` folder. Output files will be created with the root filename of `Trace_`files. They will contain Numpy arrays with single cell PWD matrices  (`_PWDscMatrix.npy`) and N-matrices (`_Nmatrix.npy`), and an `.ecsv` list of barcode names (`_unique_barcodes.ecsv`).



**Relevant options**

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

```
"colormaps":{"PWD_KDE":"terrain","PWD_median":"terrain","contact":"coolwarm","Nmatrix":"Blues"},    
```



Output images:

- `_PWDhistograms.png`
- `_Nmatrix.png`
- `HiMmatrix.png`
- `_PWDmatrixMedian.png`
- `_PWDmatrixKDE.png`

*example uses fiducial mask*

| method    | contact matrices                                             | **PWD matrix**                                               |
| --------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 2D - mask | ![image-20220212093032574](../_static/user_guide/image-20220212093032574.png) | ![image-20220212093119700](../_static/user_guide/image-20220212093119700.png) |
| 3D - mask | ![image-20220212093245315](../_static/user_guide/image-20220212093245315.png) | ![image-20220212093210913](../_static/user_guide/image-20220212093210913.png) |
| KDtree 3D | ![image-20220213120843091](../_static/user_guide/image-20220213120843091.png) | ![image-20220213120807698](../_static/user_guide/image-20220213120807698.png) |
| Nmatrices | Masking![image-20220212093324905](../_static/user_guide/image-20220212093324905.png) | KDTREE![image-20220213120921749](../_static/user_guide/image-20220213120921749.png) |



### build traces: old method

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

Parameters to run this script will be read from the ```buildsPWDmatrix``` field of ```infoList.json```.

If you want to run this function exclusively, run *pyHiM* using the ```-C buildHiMmatrix``` argument.

```sh
usage: pyHiM.py [-h] [-F ROOTFOLDER] [-C CMD] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  -C CMD, --cmd CMD     Comma-separated list of routines to run (order matters!): 
  						makeProjections, appliesRegistrations,
                        alignImages,alignImages3D, segmentMasks,
                        segmentSources3D,refitBarcodes3D,
                        localDriftCorrection,projectBarcodes,buildHiMmatrix
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

<img src="../_static/user_guide/buildsPWDmatrix_HiMmatrix.png" alt="buildsPWDmatrix_HiMmatrix" style="zoom:50%;" />

In addition, the function outputs the distribution of distances for each combination of barcodes:

<img src="../_static/user_guide/buildsPWDmatrix_PWDhistograms.png" alt="buildsPWDmatrix_PWDhistograms" style="zoom: 25%;" />


#### Filtering barcode localizations

There are several filters:
1. Properties of 2D localization algorithm (e.g. brightness)

2. Accuracy of 3D localization: sigma of fit, correlation between z-position from weighted moment and from gaussian fit, etc

3. Accuracy of drift correction in the region where the barcode was localized. 

   This is only applied if LocalDrift correction was **not** run. 


*Examples.*

| Filtering | Matrix |
| --- |  ---- |
| Unfiltered matrix. Total barcode localizations: 18700 | <img src="../_static/user_guide/buildsPWDmatrix.png" alt="buildsPWDmatrix" style="zoom:25%;" /> |
|```toleranceDrift = 1px```. Barcode localizations kept: 12377 of a total: 18700.| <img src="../_static/user_guide/buildsPWDmatrixFilterBlockDrift.png" alt="buildsPWDmatrixFilterBlockDrift" style="zoom:25%;" />|
| ```toleranceDrift = 1px```  ```Flux = 100```. Barcode localizations kept: 5562 of a total: 18700. | <img src="../_static/user_guide/buildsPWDmatrix_filterFlux100.png" alt="buildsPWDmatrix_filterFlux100" style="zoom:25%;" /> |
|```toleranceDrift = 1px```  ```Flux = 200```. Barcode localizations kept: 4528 of a total: 18700. | <img src="../_static/user_guide/buildsPWDmatrix_filterFlux.png" alt="buildsPWDmatrix_filterFlux" style="zoom:25%;" />|
|```toleranceDrift = 1px``` ```Flux = 1000```. Barcode localizations kept: 1923 of a total: 18700.| <img src="../_static/user_guide/buildsPWDmatrix_filterFlux1000.png" alt="buildsPWDmatrix_filterFlux1000" style="zoom:25%;" />|



#### Other outputs

In addition to the PWD matrix, we now also have a map of the alignment accuracy  and scatter plots showing the flux of each barcode, its sharpness, magnitude and roundness. These are used in order to validate the segmentation process and help with the selection of the ```flux``` threshold used in this filtering step.

*Alignment accuracy*

This provides a map of all barcode localizations in an ROI, colorcoded by the accuracy of localization. Colorbar scale is in pixels.

<img src="../_static/user_guide/BarcodeAlignmentAccuracy_ROI1_2D2.png" alt="BarcodeAlignmentAccuracy_ROI:1_2D2" style="zoom: 50%;" />

*Barcode localization statistics*

This provides the localization statistics from ASTROPY. The main use of these plots is to determine if the threshold ```flux``` used is correct. Default is *200*.

<img src="../_static/user_guide/BarcodeStats_ROI1_2D.png" alt="BarcodeStats_ROI:1_2D" style="zoom: 67%;" />


## Post-processing scripts

### npy_to_tiff

This script will convert Numpy array files into imageJ-readable TIFs. Images will be rescaled to (0, 2^14) range and will be histogram normalized using `skimage.exposure.equalize_adapthist()`.

You can invoke by two ways:

- Use `find` and send the list of files as arguments:

```sh
npy_to_tiff $(find -name "*ch0*_2d_registered.npy")
```

- Otherwise you can pipe the results as follows:
```sh
ls *ch0*_2d_registered.npy | npy_to_tiff
```


