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



### fileManagement.py

### functionCaller.py

### runHiM_cluster.py




## plots
- figure3wayInteractions.py
- figure4Mmatrix.py
- figureCompare2Matrices_fig3a.py
- figureCompare2Matrices_figS1n_v3.py
- figureCompare2Matrices.py
- figureHiMmatrix.py
- figureN_HiMmatrices.py
- figurePlotImageProfile.py
- figureSingleCell.py

## postProcessing
- processHiMmatrix.py
- processingMultipleDatasets.py
- processSNDchannel.py

## IA
- trainStarDist.pye''' 
