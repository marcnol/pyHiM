# Classic usage

## Default pipeline
The default pyhiM pipeline loads your images, and follow this steps :
- Project
- Align
- Detect spots
- Segment masks
- Match spots in masks
- Build matrix of pairwise distances

## Classic run

To run pyHiM in the simplest way, follow this steps :
1. Identify the folder containing the TIFF images you want to process. This folder will be called `input_directory`.
2. Copy or create a `infoList.json` file into your `input_directory`. This file contains all input parameters. You can find an example in [pyHiM/modelParameterFiles_JSON folder](https://github.com/marcnol/pyHiM/blob/master/modelParameterFiles_JSON/infoList.json).
3. Change the fiducial RT by running `changeRT_infoList.py` at the command line in the `input_directory`. The input arguments are the RT currently present in the infoList files and the RT that you want to change it for. For instance : `changeRT_infoList.py RT33 RT95`changes RT33 to RT95 in all the infoList files.
4. Still in the `input_directory`, activate your conda environment and run pyHiM by the following command at the command line :
	```bash
    conda activate pyHiM
	pyHiM.py
	```
5. When computing is done, you can find the results inside `input_directory`.