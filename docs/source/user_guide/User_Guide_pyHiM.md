# User guide of pyHiM software TODELETE

*This guide provides the key concepts of pyHiM.*

## Summary

## 1. What is pyHiM ?
### 1.1. Main goals

In the context of structural biochemistry, pyHiM software provides the **processing** of **multiplexed FISH data** produced with HiM protocol and **visualization** tools.

Basic concept of this software is to **determine 3D position** of each fluorescent spot in FISH images. These coordinates will be associated with their corresponding cell by using masks like DAPI.

We have as input **HiM experience images** and as output **position matrices** of spots.

The main use of pyHiM is to detect a set of **DNA loci**, localized with _barcodes_. In order to produce a **distance matrix** between each _barcode_. This _barcodes_ must often be detected on different images, hence the interest of _fiducials_ (see below).

FISH techniques can also be used to observe RNA and proteins. This software provides a basic way to detect a **type of RNA** but it doesn't include processing of data FISH for proteins.

### 1.2. Advantages

Each image taken during acquisition is associated with a spatial reference, the _fiducial_, common to all images. So pyHiM can :

-   **Rectifies slight shifts** from one image to another due to experimental conditions.
    
-   **Associates spots** from different images in a same area delimited by a mask.
    

Acquisition of HiM protocol experiment often generates huge volume data which makes the analysis very time consuming. So, you can run pyHiM with **parallelization** on computing servers.

### 1.3. Global structure
The pyHiM software follow a pipeline pattern. A default pipeline is defined with main features but they can be used independently if you have the right input data.

Features in the default pipeline can be classified into five categories :
1. **Preprocessing** manipulate or delete a part of data before it is used in order to ensure or enhance performance.
2. **Identification** segment areas for masks and detect coordinates for spots using starDist IA models.
3. **Matching** allow to associate spots or mask with their mask.
4. **Postprocessing** make output data easier to manipulate.
5. **Visualization** communicate progress and results to the user.

Each step can be optimized with **parallel computations** using the Dask package.

The use cases of pyHiM can be summarized in this diagram below :

```{mermaid}
flowchart LR
subgraph uc[pyHiM]
	f1([Preprocessing])
	f12([Identification])
	f7([Matching])
	f17([Postprocessing])
	f9([Visualization])
end
Biologist --- f1 & f12 & f7 & f17 & f9
uc -.- a2[Parallelization servers]
f12 --- a1[IA models]
```

## 2. Installation
To install pyhiM follow the instructions in [Installing_pyHiM file].
## 3. Quickstart
### 3.1. Basic run
To run pyHiM in the simplest way, follow this steps :
1. Identify the folder containing the TIFF images you want to process, make sure there is only this data to process in this folder and no other TIFF files or folders. This folder will be called `input_directory`.
2. Copy or create a `infoList.json` file to your `input_directory`. This file contains all input parameters. You can find a template in [modelParameterFiles_JSON folder]. For a description of `infoList.json`, see section below.
3. Change the fiducial RT by running `changeRT_infoList.py` at the command line in the `input_directory`. The input arguments are the RT currently present in the infoList files and the RT that you want to change it for. For instance : `changeRT_infoList.py RT33 RT95`changes RT33 to RT95 in all the infoList files.
4. Still in the `input_directory`, run pyHiM by the following command at the command line :
	```bash
	pyHiM.py
	```
5. When computing is done, you can find the results inside `input_directory`.

### 3.2. Separate input and output data
This last step produce output data in the same folder than your input data. If you want to separate the both, run pyhiM from an output directory with this command :
```bash
pyHiM.py -F <input_directory_path>
```
With `<input_directory_path>` the relative or absolute path to the folder containing your input data.

### 3.3. Optional arguments

If you have any doubt, you can use the help option, `pyhiM.py -h`, from anywhere to see every optional arguments and their description like this :
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



```-F ``` or ```--rootFolder``` indicates the rootFolder where pyHiM expects to find the dataset.

```--threads``` argument will make it run in parallel mode. Be ready to open your browser in ```http://localhost:8787``` and make sure you connect by ```ssh -L 8787:localhost:8787 username@servername```. Think to change your username and server name.

```-C or --cmd``` is an optional argument that can be used to run a specific set of functions detailed as a comma separated list. If you don't provide this argument, the full list of functions will be run and the mode of action will be determined from the ```infoList.json``` file (see below for details).

## 4. pyHiM fundamentals

### 4.1. Pipeline overview
#### 4.1.1. Default pyHiM flow

To run default pipeline, pyHiM need two kind of data :
- A dictionary of initialization parameters as a JSON file, name `infoList.json`
- 3D images with TIFF format, 4 types of images are accepted and they are processed in this order :
	1. Fiducial
	2. Barcode
	3. Mask (like DAPI)
	4. RNA

These types of images are called **labels**.

Steps of pipeline are executed sequentially in this order :

1. **makeProjections** : Projects all 3D images in 2D 
2. **alignImages** : Compute the best shift to align all 2D fiducials.
3. **appliesRegistrations** : Shifts 2D barcodes, masks and RNA with result of alignImages step.
4. **alignImages3D** : Takes 2D aligned fiducials and find the best shift on Z-axis. This shift will be apply on the 3D segmented barecodes at buildHiMmatrix step.
5. **segmentMasks** : Segments 2D aligned barcodes and masks with two different ways for each.
6. **segmentSources3D** : Applies 2D shift, computed at alignImages step, to 3D barcodes. Then, segments them in 3D.
7. **buildHiMmatrix** : Filters the segmentation results, associates barcode coordinates with the good mask and makes the PWD matrix for each mask.

But there are specific routines for specific labels, so this is the real running order :

|Command|Fiducial|Barcode|DAPI|RNA|
|:-:|:-:|:-:|:-:|:-:|
|makeProjections|1|4|8|12|
|alignImages|2||||
|appliesRegistrations||5|9|13|
|alignImages3D|3||||
|segmentMasks||6|10||
|segmentSources3D||7|||
|buildsPWDmatrix|||11||

#### 4.1.2. Input / Output data

Here is a table summarizing the type of input and output data for each routine :

|Routine|Input|Output|
|:-:|---|---|
|**makeProjections**|3D_raw.tif|2D_raw.npy|
|**alignImages**|2D_raw.npy + 2D_reference.npy|alignImages.ecsv|
|**AppliesRegistrations**|2D_raw.npy + alignImages.ecsv|2D_registered.npy|
|**alignImages3D**|3D_raw.tif + 3D_reference.tif + alignImages.ecsv|alignImages_block3D.ecsv|
|**segmentMasks**|2D_registered.npy|segmented_barecode.ecsv + segmented_mask.npy|
|**segmentSources3D**|3D_raw.tif + alignImages.ecsv|3D_segmented_barcode.ecsv|
|**buildHiMmatrix**|segmented_barcode.ecsv + segmented_mask.npy + alignImages_block3D.ecsv + 3D_segmented_barcode.ecsv|PWDMatrix.ecsv + 3D_PWDMatrix.ecsv|

#### 4.1.3. Flowchart of labels
```{mermaid}
flowchart TD
	subgraph graph legend
		Z1((feature))
		Z2([initial data])
		Z3[I/O data]
	end 
```

##### 4.1.3.1. Fiducial flow

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

##### 4.1.3.2. Barcode flow
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
##### 4.1.3.3. Mask flow
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
##### 4.1.3.4. RNA flow
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

### 4.2. Initialization parameters
The file of initialization parameters is called *infoList.json*. It's composed by 2 main keys :

- **labels** : fiducial / DAPI / barcode / RNA.

 Label order is important for the process.

- **common**, parameters common to the modules that pyHiM can execute :
  [label] (routine)
  1. **acquisition** [fiducial / DAPI / barcode / RNA] (all)
  1. **alignImages** [fiducial / DAPI / barcode / RNA] (alignImages / alignImages3D / appliesRegistrations)
  4. **segmentedObjects** [DAPI / barcode] (segmentMasks / segmentSources3D)
  5. **zProject** [fiducial / DAPI / barcode / RNA] (makeProjections)
  6. **buildsPWDMatrix** [DAPI] (buildHiMmatrix)

Label specific parameters can be added by adding a key in the dictionary with the module name.

Key descriptions of *common* :

- **acquisition**

  This key contains all the information related to the acquisition during the experiment.

  There are 2 types of acquisition (called *cycles*) :
  1. **DAPI (+ RNA) + fiducial** : Required but possible without RNA *channel* if not done in this experiment
  2. **(barecode + fiducial)**

	Description :
	  - DAPI_channel" : code of the DAPI for the DAPI/RNA cycle
	  - RNA_channel" : code of the RNA for the DAPI/RNA cycle
	  - "fiducialDAPI_channel" : fiducial code for the DAPI/RNA cycle
	  - "barcode_channel" : code of the barecode for the barecode cycle
	  - "fiducialBarcode_channel" : fiducial code for the barecode cycle
	  - "fileNameRegExp" : Nomenclature of the name of the data files translated in regex
	  - "pixelSizeXY" : pixel-nanometer relation for X and Y axis
	  - "pixelSizeZ" : pixel-nanometer relation for the Z axis
	  - "zBinning" : When loading the image in array to process the data, it is possible to take only one plane on two on the Z axis. For this, the value must be "2" otherwise it will be "1" to have all the planes.
	  - positionROIinformation" : To be deleted in the next version ... Allows to know the position of the ROI number in the file name
	  - parallelizePlanes" : Boolean allowing to parallelize by plans rather than by files if the parallelization is activated.


- **alignImages**

  This key contains all information related to the alignment.

  - folder": name of the folder where the aligned images are written
  - operation": unused option ?
    - overwrite
    - skip
  - "outputFile": output file name ?
  - referenceFiducial": code of reference fiducial cycle
  - "alignByBlock": If False, apply global cross-correlation, otherwise we do the one by sub-block
  - "tolerance": Used in blockAlignment to determine the % of tolerated error
  - lower_threshold": lower threshold to adjust image intensity levels before X-correlation, removes some background noise
  - higher_threshold": upper threshold to adjust the intensity levels of the image before x-correlation, avoids being biased by very bright noise
  - "background_sigma": used to remove "inhom" background?
  - "localShiftTolerance": used for 2D local alignment, if the shift exceeds this value then nothing is done
  - "bezel": pixel margin around the DAPI mask ?
  - "blockSize": size in (X,Y) used for 3D local alignment
  - "localAlignment":
    - If not used == "None".
    - if local 2D alignment == "mask2D" 
    - if 3D local alignment == "block3D"
  - "3D_lower_threshold":
  - "3D_higher_threshold":


- **segmentedObjects** TODO

- **zProject**

  Parameters concerning the projection on the Z axis.

  - folder": name of the folder for writing projected images

  - operation": unused option ?

    - overwrite
    - skip

  - "mode": 

    - full: projects everything in 2D depending on "zProjectOption"

    - manual : take the planes between the values entered in "zmin" and "zmax" then apply the "zProjectOption

    - automatic : estimates the focal plane using the maximum standard deviation plane by plane. Then, we project (according to "zProjectOption") around this focal plane more or less the value entered in "zwindows".

    - laplacian : 

      Divides the image into blocks of the size given by "blockSize".
      Then we look for the maximum of the "laplacian variance" (blur estimation) for each block in order to estimate the focal plane.

      Finally, we perform a block by block reconstruction of the 2D image using the optimal focal plane of each block. If "zwindows" != 0, the MIP option of "zProjectOption" is necessarily applied.

  - display": display in terminal or logging markdown ?

  - "blockSize": used for *laplacian* mode

  - "saveImage": boolean to keep the projected image, just logging?

  - "zmin": minimum for *manual mode

  - zmax": maximum for *manual mode*.

  - "zwindows": used for both *automatic* and *laplacian* modes

  - windowSecurity": used for *automatic* mode, removes the lowest and highest X planes on the Z axis

  - zProjectOption": used for all modes ?

    - sum : sum of all concerned planes
    - MIP : *maximum intensity projection*, takes the maximum of the concerned planes, pixel by pixel.

- **buildsPWDMatrix** TODO

### 4.3. Main features
#### 4.3.1. makeProjections
*Projects 3D images in 2D*

Initialization parameters :

 |Name|Option|Description|
 |:-:|:-:|:-:|
 |mode|manual|Assign plans between *zmin* and *zmax* to "zRange"|
 ||automatic|Estimates the focal plane using the maximum standard deviation plane by plane. Use *zwindows* to set "zRange" arround this focal plane.|
 ||full|Assign all plans to "zRange"|
 ||laplacian|Split 3D image into blocks of the size given by *blockSize*. Find Laplacian Variance maximum (blur estimation) for each block in order to estimate the focal plane. Rebuild block-by-block 2D image with optimal focal plane of each block. if *zwindows* option is activated, project each block with MIP option.|
 |windowSecurity||Used for *automatic* mode, removes the lowest and highest Z-plans.|
 |zProjectOption|sum|Sum plans in "zRange"|
 ||MIP|Maximum Intensity Projection of plans in "zRange"|

This routine take all 3D images and project its in 2D.
Depending on the chosen *mode*, this feature start to find the good set of Z-plans, where there is the least noise. This step give a range centered on a focal plan, named *zRange*. After, projection is done on this range either by sum or by maximum intensity projection.

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

#### 4.3.2. alignImages
*Registers fiducials using a barcode as reference*

In the set of *fiducial* images, one is chosen by initialization parameters to be the reference. 
The algorithm takes images one by one and align with the reference.
There are two ways to compute the shift :
- Global alignement make simple cross-correlation with tow images
- Split image in block and make cross-correlation block by block. Then we have one shift by block and to align the global image an average of those shifts are made. This method is more robust against a bright noise spot.
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

#### 4.3.3. AppliesRegistrations
*Applies registration to DAPI and barcodes*
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

#### 4.3.4. alignImages3D
*Aligns fiducials in 3D*
This feature run for *fiducial* images and with "block3D" value for "localAlignment" key in *infoList.json* file.
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
#### 4.3.5. segmentMasks
*Segments DAPI and sources in 2D*
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

#### 4.3.6. segmentSources3D
*Segments sources in 3D*

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

#### 4.3.7. buildHiMmatrix TOCHANGE
*Builds PWD matrix for all folders with images*
Assigns barcode localizations to DAPI masks and constructs single cell cummulative PWD matrix.
```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([.])
		end
		subgraph C[OUTPUT]
			C1([.])
		end
		A --> C
		
	
```

### 4.4. Secondary features
#### 4.4.1. fileProcessing
- changeRT_infoList.py
- cleanHiM_run.py
- Indir.py
- zipHiM_run.py

#### 4.4.2. plots
- figure3wayInteractions.py
- figure4Mmatrix.py
- figureCompare2Matrices_fig3a.py
- figureCompare2Matrices_figS1n_v3.py
- figureCompare2Matrices.py
- figureHiMmatrix.py
- figureN_HiMmatrices.py
- figurePlotImageProfile.py
- figureSingleCell.py

#### 4.4.3. postProcessing
- processHiMmatrix.py
- processingMultipleDatasets.py
- processSNDchannel.py

#### 4.4.4. IA
- trainStarDist.pye''' 