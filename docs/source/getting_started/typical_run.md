# Running pyHiM



## Default pipeline

The default *pyHiM* pipeline will:
- Project images in 3D

- Calculate drift based on fiducial images

- Apply drift to nuclei and barcode images (i.e. registration)

- Segment and localize barcode spots

- Segment nuclei and masks

- Match spots to masks to build chromatin traces

- Build single-trace and ensemble pairwise distance maps

  

## Typical run

To run *pyHiM* in the simplest way, follow these steps:
1. Copy the TIFF images you want to process into a single folder. This folder will be referred as the `input_directory`.
2. Copy or create a model configuration file (called `infoList.json`) into your `input_directory`. This file contains all the input parameters required to run `pyHiM`. You can find an example in [pyHiM/modelParameterFiles_JSON folder](https://github.com/marcnol/pyHiM/blob/master/modelParameterFiles_JSON/infoList.json).
3. Modify the `infoList.json` file to indicate the reference cycle used for drift correction. This can be done by manually editing the `infoList.json` file or by running the graphical user interface provided in the script: `function_parameters.py`.
4. In the `input_directory` activate your conda environment and run pyhim by typing the following commands in a terminal (Linux) or in the anaconda prompt terminal (Windows):
	```bash
   conda activate pyHiM
	pyhim
	```
5. Each module will store its results in a single dedicated folder in the `input_directory`. A summary of all outputs will also be provided as a markdown file `HiM_Analysis<DDMMYYYY_HHMMSS>.md` saved in the `input_directory`.

## Tutorials

Tutorials are provided as Jupyter labs and can be found [here](https://github.com/marcnol/pyHiM/tree/notebooks/notebooks).