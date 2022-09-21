# Running pyHiM

## Prepare run

1. Copy the TIFF images you want to process into a single folder. This folder will be referred as the `input_directory`.
2. Copy or create a model configuration file (called `infoList.json`) into your `input_directory`. This file contains all the input parameters required to run `pyHiM`. You can find an example in [pyHiM/modelParameterFiles_JSON folder](https://github.com/marcnol/pyHiM/blob/master/modelParameterFiles_JSON/infoList.json).
3. Modify the `infoList.json` file to indicate the reference cycle used for drift correction. This can be done by manually editing the `infoList.json` file or by running the graphical user interface provided in the script: `function_parameters.py` ([tutorial here](/insert/link/tuto/param)).
4. In the `input_directory` activate your conda environment by typing the following command in a terminal (Linux) or in the anaconda prompt terminal (Windows):
	```bash
   conda activate pyHiM
	```

## Basic run

4. To run without any argument, type the following command in the `input_directory`:
	```bash
	pyhim
	```

5. Each module will store its results in a single dedicated folder in the `input_directory`. A summary of all outputs will also be provided as a markdown file `HiM_Analysis<DDMMYYYY_HHMMSS>.md` saved in the `input_directory`.

## Optional arguments

If you require help, you can call `pyHiM` with the help option as follows: 
```bash
pyhim -h
```

The output should look as follows:

```sh
usage: pyhim [-h] [-F ROOTFOLDER] [-C CMD] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  -C CMD, --cmd CMD     Comma-separated list of routines to run (order matters!):
                            makeProjections alignImages appliesRegistrations
                            alignImages3D segmentMasks segmentMasks3D
                            segmentSources3D buildHiMmatrix
                        optional:
                            filter_localizations register_localizations
                            build_traces build_matrix
  --threads THREADS     Number of threads to run in parallel mode. If none,
                        then it will run with one thread.
```



- ```-F ``` or ```--rootFolder``` indicates the rootFolder where *pyHiM* expects to find the dataset.

- ```-C or --cmd``` is an optional argument that can be used to run a specific set of functions detailed as a comma separated list. If you don't provide this argument, the full list of functions will be run and the mode of action will be determined from the ```infoList.json``` configuration file.

- ```--threads``` will ask *pyHiM* to run in parallel using multiple threads in your computer or computer cluster. To visualize the progress of your run,  open your browser in ```http://localhost:8787``` and make sure you connect by ```ssh -L 8787:localhost:8787 username@servername``` if you are not running `pyHiM` locally.

## Running in multiple folders
To run `pyHiM` in multiple folders we recommend that you create a BASH script and provide the location of each `input_directory` that needs to be processed. For this, run `pyHiM`:
```bash
pyhim -F <input_directory>
```
Where `<input_directory>` is the relative or absolute path to the folder containing your input data.


## Default pipeline

The basic *pyHiM* pipeline will:

- Project images in 3D

- Calculate drift based on fiducial images

- Apply drift to nuclei and barcode images (i.e. registration)

- Segment and localize barcode spots

- Segment nuclei and masks

- Match spots to masks to build chromatin traces

- Build single-trace and ensemble pairwise distance maps