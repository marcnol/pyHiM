# General-use scripts


## File Processing, handling HPC runs, etc

### changeRT_infoList.py

This script is used to modify the reference fiducial signal to be used for shift correction in the parameter file

```
usage: cahngeRT_infolist.py old_RT new_RT
```
Example: 

changeRT_infolist.py RT33 RT1

### cleanHiM_run.py

Cleans the directories and log files created by pyHiM from a previous run. 

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

Creates links for files in a directory into a second directory. Useful to analyze data in a new folder without copying the raw data files.

```
Usage: lndir.py "/home/marcnol/Repositories/pyHiM/\*py" ~/Downloads/test
```

Use quotation marks in the first argument if using wildcards.

### zipHiM_run.py

Zip all output files from a pyHiM run. It excludes .npy and .tif files. Useful to retrieve results from a run from an HPC cluster.

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

Unzips HiM_run.tar.gz recursively. Useful to unzip the results from several folders retrieved from a run in an HPC cluster.

```
Usage: unzipHiM_run.py [-F ROOTFOLDER] [-R RECURSIVE]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the HiM_run.tar.gz is located
  -R RECURSIVE, --recursive RECURSIVE
                        Unzip files inside folders of current directory
```

### runHiM_cluster.py

Launches pyHiM on a cluster using slurm `srun`. 

```
Usage: runHiM_cluster.py
```

## Plotting scripts


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


### figure4Mmatrix.py
Creates proximity frequency 4M profiles from a given list of anchors (similar analysis to a 4C experiment, but using HiM data). Works with up to two datasets.

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

### figureCompare2Matrices.py
Comparison of proximity matrices. Plots either the ratio between two HiM matrices, or the difference. It also plots both matrices together, with one in the upper triangle, and the other in the lower triangle.

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

### figure3wayInteractions.py

Plots 3-way proximity probability matrices for a given anchor (or set of anchors), as defined in the folders2Load.json configuration file. Comparative analysis can be performed for two datasets simultaneously. The calculation of 3-way proximity probability matrices needs to be previously performed using the `processHiMmatrix.py` script. 

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



### figureN_HiMmatrices.py
Plots several (`N`) HiM matrices in the same plot, using `N` datasets specified in folders2Load.json.

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

### figureSingleCell.py
This scripts: 
- produces movies and trajectories from single cell PWD matrices. 
- calculates barcode detection efficiencies and number of barcodes per cell.
- plots single cell matrices.
- plots distance histograms and distributions of Rg.


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

## Matrix post-processing scripts

### processHiMmatrix.py

This script performs the post-processing of one or more previously pyHiM-analysed datasets, defined in folders2Load.json file. 

It performs the following operations:
- Merges datasets from different experiments.
- Calculates and plots ensemble pairwise distance (PWD) matrix.
- Calculates and plots the inverse of the PWD matrix. 
- Calculates and plots the contact probability matrix for each dataset.
- Calculates and plots the ensemble contact probability matrix.
- Calcualtes and plots the ensemble 3-way contact probability matrix, for the set of anchors defined in folders2Load.json file. 
- Optional: Read MATLAB single-cell PWD matrices, and perform all previous operations.

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
			Saves the combined PWD matrix from all datasets. Default: False
	--getStructure
			Multi-dimensional scaling to get coordinates from PWDs. Default: False
	--pixelSize
			Specify images pixel size. Default: 100 nm.
	--HiMnormalization
			Normalization of contact matrix: nonNANs (default) or nCells.
	--d3
			Loads data segmented in 3D. Default: False
```

### processSNDchannel.py

This script will:
- allow the user to manually draw ROI based on secondary labels, such as RNA-FISH images.
- use the ROIs defined by the user to attribute labels to traces.


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

Script used to 2D-project the 3D-labeled numpy array images generated by segmentMasks3D. By running this script, you will replace the output files produced by segmentMasks (2D labeled masks).

```
Usage: process_segmentMasks3D.py
```

### trace_combinator.py
This script combines trace tables from different experiments/ROIs together into a single trace table.
The folders containing the trace tables of the experiments to be combined are provided in as a JSON file. The possibility to selecting only a subset of trace tables within the folders provided is allowed using the `methods` parameter. Merged trace table is outputed to the buildPWDmatrix folder. 

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

This scipt loads a trace file and a number of numpy masks, and assings them the labels produced by `processSNDchannel.py`.

```
Usage: trace_selector.py [-F ROOTFOLDER] [--pixel_size]

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		Folder with fimages
	--pixel_size
		Lateral pixel size in microns. Default = 0.1
```

