# Running pyHiM



## Default pipeline

The default pyHiM pipeline will:
- Project images in 3D

- Calculate drift based on fiducial images

- Apply drift to nuclei and barcode images (i.e. registration)

- Segment and localize barcode spots

- Segment nuclei and masks

- Match spots to masks to build chromatin traces

- Build single-trace and ensemble pairwise distance maps

  

## Typical run

To run pyHiM in the simplest way, follow this steps:
1. Copy the TIFF images you want to process into a single folder. This folder will be called `input_directory`.
2. Copy or create a model configuration file (called `infoList.json`) into your `input_directory`. This file contains all input parameters required to ryb `pyHiM`. You can find an example in [pyHiM/modelParameterFiles_JSON folder](https://github.com/marcnol/pyHiM/blob/master/modelParameterFiles_JSON/infoList.json).
3. Modify the `infoList.json` file to indicate the cycle used for drift correction. This can be used by manually editing the `infoList.json` file or by running the graphical user interface provided in the script: `function_parameters.py`.
4. In the `input_directory` activate your conda environment and run pyHiM by the following command at the command line:
	```bash
   conda activate pyHiM
	pyHiM.py
	```
5. When computing is done, the results will be stored within the folders created at the `input_directory`. Each module or feature will store its results in a single folder. A summary of all the results and outputs is provided in a markdown file that will be created at the `input_directory`.