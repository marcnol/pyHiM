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


## Post-processing scripts

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
