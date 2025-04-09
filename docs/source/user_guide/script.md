# General-use scripts

## File Processing, handling HPC runs, etc

### cleanHiM_run.py

*Cleans the directories and log files created by pyHiM in previous runs.*

```
Usage: clean_him_run [-F ROOTFOLDER] [-P PARAMETERS] [-A ALL]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the analysis has been performed
  -P PARAMETERS, --fileParameters PARAMETERS
                        Parameters file. Default: parameters.json
  -A ALL, --all ALL
                        Delete all folders and all created files
```
### lndir.py

*Creates link for files in a second directory (useful to analyze data in a new folder without copying raw data files).*

```
Usage: lndir "/user_home/Repositories/pyHiM/*py" ~/Downloads/test
```

Use quotation marks in the first argument if using wildcards.

### zipHiM_run.py

Zip all output files from a *pyHiM* run. It excludes .npy and .tif files. Useful to retrieve results from a run from an HPC cluster.

```
Usage: zip_him_run [-F ROOTFOLDER] [-P PARAMETERS] [-R RECURSIVE]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the analysis has been performed
  -P PARAMETERS, --fileParameters PARAMETERS
                        Parameters file. Default: parameters.json
  -R RECURSIVE, --recursive RECURSIVE
                        Zip files inside folders of current directory
```

### unzipHiM_run.py

Unzips HiM_run.tar.gz recursively. Useful to unzip the results from several folders retrieved from a run in an HPC cluster.

```
Usage: unzip_him_run [-F ROOTFOLDER] [-R RECURSIVE]

optional arguments:
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder where the HiM_run.tar.gz is located
  -R RECURSIVE, --recursive RECURSIVE
                        Unzip files inside folders of current directory
```

## Plotting scripts

### figureCompare2Matrices.py
Comparison of proximity matrices. Plots either the ratio or the difference between two HiM matrices. It also plots both matrices together, with one in the upper triangle, and the other in the lower triangle.

```
Usage: figure_compare_2_matrices [-F1 ROOTFOLDER1] [-F2 ROOTFOLDER2]
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
		 Performs the ratio between matrices. Default: difference
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

## Post-processing scripts

### processSNDchannel.py

This script will:
- allow the user to manually draw ROI based on secondary labels, such as RNA-FISH images.
- use the ROIs defined by the user to attribute labels to traces.


```
Usage: process_snd_channel [-F ROOTFOLDER] [-A ADDMASK] [--cleanAllMasks]

	-F ROOTFOLDER, --rootFolder ROOTFOLDER
		Folder with images
	-A ADDMASK, --addMask ADDMASK
		Add manual segmentation
	--cleanAllMasks
		Clear all masks
```

### npy_to_tiff

This script will convert Numpy array files into imageJ-readable TIFs. Images will be rescaled to (0, 2^14) range and will be histogram normalized using `skimage.exposure.equalize_adapthist()`.

You can invoke this in two ways:

- Use `find` and send the list of files as arguments:

```sh
npy_to_tiff $(find -name "*ch0*_2d_registered.npy")
```

- Otherwise you can pipe the results as follows:
```sh
ls *ch0*_2d_registered.npy | npy_to_tiff
```
