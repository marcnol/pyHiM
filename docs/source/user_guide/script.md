# WIP - Secondary scripts


## fileProcessing Scripts

### changeRT_infoList.py

This script is used to modify the reference fiducial signal to be used for shift correction in all parameters files (.json files)

```
usage: cahngeRT_infolist.py old_RT new_RT
```
Example: 

changeRT_infolist.py RT33 RT1

### cleanHiM_run.py

Cleans the directories created by pyHiM in the analysis folder. 

```
Usage: cleanHiM_run.py [-F ROOTFOLDER] [-P PARAMETERS] [-A ALL]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the analysis has been performed
  -P PARAMETERS, --fileParameters PARAMETERS
                        Parameters file. Default: infolist.json
  -A ALL, --all ALL
                        Delete all folders and all created files
```
### lndir.py

Creates links for files in a directory into a second directory.

```
Usage: lndir.py "/home/marcnol/Repositories/pyHiM/\*py" ~/Downloads/test
```

Use quotation marks in the first argument if using wildcards.

### zipHiM_run.py

Zip all output files from a pyHiM run. It excludes .npy and .tif files. 

```
Usage: zipHiM_run.py [-F ROOTFOLDER] [-P PARAMETERS] [-R RECURSIVE]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the analysis has been performed
  -P PARAMETERS, --fileParameters PARAMETERS
                        Parameters file. Default: infolist.json
  -R RECURSIVE, --recursive RECURSIVE
                        Zip files inside folders of current directory
```

### unzipHiM_run.py

Unzips HiM_run.tar.gz recursively

```
Usage: unzipHiM_run.py [-F ROOTFOLDER] [-R RECURSIVE]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the HiM_run.tar.gz is located
  -R RECURSIVE, --recursive RECURSIVE
                        Unzip files inside folders of current directory
```

### runHiM_cluster.py

Launches pyHiM on a cluster using slurm srun job. 

```
Usage: runHiM_cluster.py
```


### fileManagement.py --> ?? We dont use it directly.
Contains classes and functions for file management.

### functionCaller.py --> ?? We dont use it directly
Contains classes, methods and functions that call the functions that perform different steps in the pyHiM pipeline (for example, segmentSources3D).

## plots
### figure3wayInteractions.py

Plots 3-way contact probability matrices for a given anchor (or set of anchors), defined in folders2Load.json file. This can be done for two datasets simultaneously. 
The calculation of 3-way contact probability matrices needs to be previously done using processHiMmatrix.py script. 

```
Usage: figure3wayInteractions.py [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2] 
								 [-O OUTPUTFOLDER] [-P PARAMETERS]
								 [-P2 PARAMETERS2] [-A1 LABEL1] [-A2 LABEL2]
								 [-W1 ACTION1] [-W2 ACTION2] [--fontsize]
								 [--scalingParameter] [--colorbar] 
								 [--plottingFileExtension] [--normalize]
	
	-F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
		 Folder with dataset 1
	-F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
		 Folder with dataset 2
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -P2 PARAMETERS2, --parameters2 PARAMETERS2
		 Name of parameters file for dataset 2. Default: folders2Load.json
	 -A1 LABEL1, --label1 LABEL1
		 Name of label for dataset 1
	 -A2 LABEL2, --label2 LABEL2
		 Name of label for dataset 2
	 -W1 ACTION1, --action1 ACTION1
		 Selects: all, labeled or unlabeled for dataset 1
	 -W2 ACTION2, --action2 ACTION2
		 Selects: all, labeled or unlabeled for dataset 2
	 --fontsize
		 Size of fonts to be used in plots
	 --scalingParameter
		 Scaling parameter of the colormap
	 --colorbar
		 Use if a colorbar is required
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	 --normalize
		 Normalizes matrices by their maximum.

```


### figure4Mmatrix.py
Creates 4M profiles of interaction freaquency for a given list of anchors (similar analysis to a 4C experiment, but using HiM data). Works with up to two datasets.

```
Usage: figure4Mmatrix.py [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2] [-O OUTPUTFOLDER]
						 [-P PARAMETERS] [-A1 LABEL1] [-A2 LABEL2] [-W1 ACTION1]
						 [-W2 ACTION2] [--fontsize] [--axisLabel] [--axisTicks]
						 [--splines] [--cAxis] [--plottingFileExtension]
						 [--legend] [--normalize]
 

	-F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
		 Folder with dataset 1
	-F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
		 Folder with dataset 2
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -A1 LABEL1, --label1 LABEL1
		 Name of label for dataset 1
	 -A2 LABEL2, --label2 LABEL2
		 Name of label for dataset 2
	 -W1 ACTION1, --action1 ACTION1
		 Selects: all, labeled or unlabeled for dataset 1
	 -W2 ACTION2, --action2 ACTION2
		 Selects: all, labeled or unlabeled for dataset 2
	 --fontsize
		 Size of fonts to be used in plots
	 --axisLabel
		 Select optional label in x and y axis
	 --axisTicks
		 Display axis ticks
	 --splines 
		 Plots data using spline interpolations
	 --cAxis
		 Absolute axis value for colormap
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	 --legend
		 Shows legends for datasets in plot
	 --normalize
		 Matrix normalization factor: maximum, none, single value. Default: none
```

### figureCompare2Matrices_fig3a.py
Plots either the ration between two HiM matrices, or the difference. It also plots both matrices together, with one in the upper triangular part, and the other in the lower triangular part. 

```
Usage: figureCompare2Matrices_fig3a.py [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2] 
									   [-O OUTPUTFOLDER] [-P PARAMETERS]
									   [-A1 LABEL1] [-A2 LABEL2] [-W1 ACTION1]
									   [-W2 ACTION2] [--fontsize] [--axisLabel] 
									   [--axisTicks] [--ratio] [--cAxis] 
									   [--plottingFileExtension] [--shuffle1]
									   [--shuffle2] [--cMinMax]
 

	-F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
		 Folder with dataset 1
	-F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
		 Folder with dataset 2
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -A1 LABEL1, --label1 LABEL1
		 Name of label for dataset 1
	 -A2 LABEL2, --label2 LABEL2
		 Name of label for dataset 2
	 -W1 ACTION1, --action1 ACTION1
		 Selects: all, labeled or unlabeled for dataset 1
	 -W2 ACTION2, --action2 ACTION2
		 Selects: all, labeled or unlabeled for dataset 2
	 --fontsize
		 Size of fonts to be used in plots
	 --axisLabel
		 Select optional label in x and y axis
	 --axisTicks
		 Display axis ticks
	 --ratio 
		 Performs the ratio between matrices. Defaukt: difference
	 --cAxis
		 Absolute axis value for colormap
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	 --shuffle1
		 Provide shuffle vector of the same size or smaller than the original
		 matrix for dataset 1. The vector should be formatted as: 0,1,2,... 
		 with no spaces, and comma-separated
	 --shuffle2
		 Provide shuffle vector of the same size or smaller than the original
		 matrix for dataset 2. The vector should be formatted as: 0,1,2,... 
		 with no spaces, and comma-separated
	 --cMinMax
		 Define min and max value for the colormap. It should be 
		 comma-separated and with no spaces. For example: 0,0.5
```

### figureCompare2Matrices_figS1n_v3.py
Plots a comparison matrix between nuclei displaying at least one RNA-FISH signal, and a subset of 33% nuclei displaying the brightest RNA-FISH signals.

```
Usage: figureCompare2Matrices_figS1n_v3.py [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2] 
										   [-O OUTPUTFOLDER] 
 

	-F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
		 Folder with dataset 1
	-F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
		 Folder with dataset 2
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 
```

### figureCompare2Matrices.py
Plots either the ration between two HiM matrices, or the difference. It also plots both matrices together, with one in the upper triangular part, and the other in the lower triangular part.

```
Usage: figureCompare2Matrices.py [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2] 
								 [-O OUTPUTFOLDER] [-P PARAMETERS]
							     [-A1 LABEL1] [-A2 LABEL2] [-W1 ACTION1]
							     [-W2 ACTION2] [--fontsize] [--axisLabel] 
							     [--axisTicks] [--ratio] [--cAxis] 
							     [--plottingFileExtension] [--normalize]
							     [--inputMatrix] [--pixelSize]
 

	-F1 ROOTFOLDER1, --rootFolder1 ROOTFOLDER1
		 Folder with dataset 1
	-F2 ROOTFOLDER2, --rootFolder2 ROOTFOLDER2
		 Folder with dataset 2
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -A1 LABEL1, --label1 LABEL1
		 Name of label for dataset 1
	 -A2 LABEL2, --label2 LABEL2
		 Name of label for dataset 2
	 -W1 ACTION1, --action1 ACTION1
		 Selects: all, labeled or unlabeled for dataset 1
	 -W2 ACTION2, --action2 ACTION2
		 Selects: all, labeled or unlabeled for dataset 2
	 --fontsize
		 Size of fonts to be used in plots
	 --axisLabel
		 Select optional label in x and y axis
	 --axisTicks
		 Display axis ticks
	 --ratio 
		 Performs the ratio between matrices. Defaukt: difference
	 --cAxis
		 Absolute axis value for colormap
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png		
	 --normalize
		 Matrix normalization factor: maximum, none, single value, 
		 bin pair. Default: none
	 --inputMatrix
		 Source of input matrix: contact (default), PWD matrix, 
		 iPWD matrix
	 --pixelSize
		 Pixel size in microns. Default: 0.1 microns
		 
```

### figureHiMmatrix.py

Produces and plots a HiM matrix for a given dataset.

```
Usage figureHiMmatrix.py [-F ROOTFOLDER] [-O OUTPUTFOLDER] [-P PARAMETERS] 
						 [-A LABEL] [-W ACTION] [--fontsize] [--axisLabel]
						 [--axisTicks] [--barcodes] [--scalingParameter]
						 [--cScale] [--plottingFileExtension] [--shuffle]
						 [--scalogram] [--inputMatrix] [--pixelSize]
						 [--cmap] [--PWDmode]						 
	 -F ROOTFOLDER, --rootFolder ROOTFOLDER
		 Folder with datasets
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -A LABEL, --label LABEL
		 Name of label
	 -W ACTION, --action ACTION
		 Selects: all, labeled or unlabeled for the datasets.
	 --fontsize
		 Size of fonts to be used in plots
	 --axisLabel
		 Select optional label in x and y axis
	 --axisTicks
		 Display axis ticks
	 --barcodes
		 Display barcode images
	 --scalingParameter
		 Scaling parameter of colormap
	 --cScale
		 Colormap absolute scale
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	 --shuffle
		 Provide shuffle vector: 0,1,2,3,.. of the same size or
		 smaller than the original matrix. 
	 --scalogram
		 Display scalogram image
	 --inputMatrix
		 Select plot type among one of the following: PWD, contact, iPWD.
		 Default: contact
	 --pixelSize
		 Pixel size in µm
	 --cmap
		 Select colormap. Default: coolwarm
	 --PWDmode
		 Mode used to calculate the mean distance. 
		 Options are: 'median' or 'KDE'. Default: 'median'
```


### figureN_HiMmatrices.py
Plots N HiM matrices in the same plot, using N datasets specified in folders2Load.json.

It also plots a submatrix containing the difference of contact probability for a subset of barcodes with respect to a particular dataset. The subset of barcodes and the reference dataset are defined in folders2Load.json by the options "barcodes2plot" and "plotSegment_anchor" respectively. 

```
Usage figureN_HiMmatrices.py [-F ROOTFOLDER] [-O OUTPUTFOLDER] [-P PARAMETERS] 
							 [-A LABEL] [-W ACTION] [--fontsize] [--axisLabel]
							 [--axisTicks] [--barcodes] [--scalingParameter]
							 [--plottingFileExtension] [--shuffle] [--scalogram]
							 [--type] [--pixelSize] [--cAxis] [--ratio]
							 [--normalizeMatrix]
	 -F ROOTFOLDER, --rootFolder ROOTFOLDER
		 Folder with datasets
	 -O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	 -P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	 -A LABEL, --label LABEL
		 Name of label
	 -W ACTION, --action ACTION
		 Selects: all, labeled or unlabeled for the datasets.
	 --fontsize
		 Size of fonts to be used in plots
	 --axisLabel
		 Select optional label in x and y axis
	 --axisTicks
		 Display axis ticks
	 --barcodes
		 Display barcode images
	 --scalingParameter
		 Scaling parameter of colormap
	 --plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	 --shuffle
		 Provide shuffle vector: 0,1,2,3,.. of the same size or
		 smaller than the original matrix. 
	 --scalogram
		 Display scalogram image
	 --type
		 Select plot type among one of the following: PWD, contact, iPWD
	 --pixelSize
		 Pixel size in µm
	 --cAxis
		 Absolute axis value for colormap
	 --ratio
		 Calculates ration between matrices for submatrices plots. 
		 Default: difference
	 --normalizeMatrix
		 Normalize matrices by maximum. Default: True
```

### figurePlotImageProfile.py ???
Plots a line profile from an image in npy format. Loads a maximum intensity projection of a barcode image (2D image). 

```
Usage: figurePlotImageProfile.py
```

### figureSingleCell.py ??
Produces movies and structures from single cell PWD matrices. 

```
Usage: figureSingleCell.py [-F ROOTFOLDER] [-O OUTPUTFOLDER] [-P PARAMETERS] 
						   [-A LABEL] [-W ACTION] [--fontsize] [--axisLabel]
						   [--axisTicks] [--barcodes] [--nRows] [--pixelSize]
						   [--maxDistance] [--plottingFileExtension] [--shuffle]
						   [--ensembleMatrix] [--video] [--videoAllcells]
						   [--plotHistogramMatrix] [--minNumberPWD] [--threshold]
	
	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		 Folder with datasets
	-O OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
		 Folder for outputs
	-P PARAMETERS, --parameters PARAMETERS
		 Name of parameters file. Default: folders2Load.json
	-A LABEL, --label LABEL
		 Name of label
	-W ACTION, --action ACTION
		 Selects: all, labeled or unlabeled for the datasets.
	--fontsize
		 Size of fonts to be used in plots
	--axisLabel
		 Select optional label in x and y axis
	--axisTicks
		 Display axis ticks
	--barcodes
		 Display barcode images
	--nRows
		 The number of cells is determined by nRows^2. Default: 10
	--pixelSize
		 Pixel size in microns. Default: 0.1 microns
	--maxDistance
		 Maximum distance for histograms in microns. Default: 4 microns
	--plottingFileExtension
		 Select file extension to save images. Default: svg. 
		 Other options: pdf, png
	--shuffle
		 Provide shuffle vector: 0,1,2,3,.. of the same size or
		 smaller than the original matrix. 
	--ensembleMatrix
		Use if ensemble matrix should be plot alongside single cell
		matrices
	--video
		Use if you want to output video
	--videoAllcells
		Use if you want all nRows^2 single cells to be output in video
	--plotHistogramMatrix
		Use if you want to plot the PWD histograms for all bin combinations
	--minNumberPWD
		Minimum number of PWD to calculate radius of gyration Rg. Default: 6
	--threshold
		Maximum accepted PWD (in pixels) to calculate radius of gyration Rg.
		Default: 8
```

## postProcessing

### processHiMmatrix.py

This script performs the post-processing of one or more previously pyHiM-analysed datasets, defined in folders2Load.json file. 

It performs the following operations:

* Calculates and plots ensemble pairwise distance (PWD) matrix.
* Calculates and plots the inverse of the PWD matrix. 
* Calculates and plots the contact probability matrix for each dataset.
* Calculates and plots the ensemble contact probability matrix.
* Calcualtes and plots the ensemble 3-way contact probability matrix, for the set of anchors defined in folders2Load.json file. 
* Optional: Read MATLAB single-cell PWD matrices, and perform all previous operations.

```
Usage: processHiMmatrix.py [-F ROOTFOLDER] [-P PARAMETERS] [-A LABEL] [-W ACTION]
						   [--matlab] [--saveMatrix] [--getStructure] [--pixelSize]
						   [--HiMnormalization] [--d3]
Optional arguments:

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
			Folder with folders2Load.json file
	-P PARAMETERS, --parameters PARAMETERS
			File with parameters. Default: folders2Load.json
	-A LABEL, --labal LABEL
			Name of label for the dataset
	-W ACTION, --action ACTION
			Selects: all, labeled or unlabeled for the datasets. 
	--matlab
			Loads MATLAB data (e.g. .mat files)
	--saveMatrix
			????
			Saves the combined PWD matrix from all datasets. Default: False
	--getStructure
			????
			Multi-dimensional scaling to get coordinates from PWDs. Default: False
	--pixelSize
			Specify images pixel size. Default: 100 nm.
	--HiMnormalization
			Normalization of contact matrix: nonNANs (default) or nCells.
	--d3
			Loads data segmented in 3D. Default: False
```

### processingMultipleDatasets.py ??

Script to process several datasets at once. 

```
Usage: processingMultipleDatasets.py rootfolder datasetID 

	rootfolder 
		Directory containing the files to analyse. 
	datasetID 
		Number of the datasets to be analysed, e.g. 1 2 3 will analyse dataset_1 
		and dataset_2
```


### processSNDchannel.py

Process secondary masks for RNA label

```
Usage: processSNDchannel.py [-F ROOTFOLDER] [-A ADDMASK] [--cleanAllMasks]

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		Folder with images
	-A ADDMASK, --addMask ADDMASK
		Add manual segmentation
	--cleanAllMasks
		Clear all masks
```

### process_segmentMasks3D.py

Projects the 3D labeled numpy arrays produced by segmentMasks3D, and replaces those produced by segmentMasks

```
Usage: process_segmentMasks3D.py
```

### trace_combinator.py
This script takes the JSON file with folders where datasets are stored. It searches for Trace files calculated with the expected methods, loads them, and combines them into a single table that is outputed to buildPWDmatrix folder. 

Outputs: ChromatinTraceTable() object, and output .ecsv formatted file with assembled trace tables.

```
Usage: trace_combinator.py [-F ROOTFOLDER] [-P PARAMETERS] [-A LABEL] [-W ACTION]
						   [--saveMatrix] [--ndims] [--method]

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		Folder with folders2Load.json file
	-P PARAMETERS, --parameters PARAMETERS
		File with parameters. Default: folders2Load.json
	-A LABEL, --labal LABEL
		Name of label for the dataset
	-W ACTION, --action ACTION
		Selects: all, labeled or unlabeled for the datasets. 
	--saveMatrix
		Saves the combined PWD matrix from all datasets. Default: False
	--ndims
		Dimensions of the trace (2 or 3). Default: 3
	--method
		Method or mask ID used for tracing: KDtree, mask, mask0
```

### trace_selector.py

This scipt loads a trace file and a number of numpy masks, and assings them labels.

```
Usage: trace_selector.py [-F ROOTFOLDER] [--pixel_size]

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		Folder with fimages
	--pixel_size
		Lateral pixel size in microns. Default = 0.1
```

## IA???
- trainStarDist.pye''' 
