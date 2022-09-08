# Personalized runtime 

## Running pyHiM in multiple folders
To run `pyHiM` in multiple folders we recommend that you create a BASH script and provide the location of each `input_directory` that needs to be processed. For this, run `pyHiM`:
```bash
pyHiM.py -F <input_directory>
```
Where `<input_directory>` is the relative or absolute path to the folder containing your input data.

## Optional arguments

If you require help, you can call `pyHiM` with the help option as follows: `pyhiM.py -h`. The output should look as follows:
```sh
usage: pyHiM.py [-h] [-F ROOTFOLDER] [-C CMD] [--threads THREADS]

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



```-F ``` or ```--rootFolder``` indicates the rootFolder where *pyHiM* expects to find the dataset.

```--threads``` will ask *pyHiM* to run in parallel using multiple threads in your computer or computer cluster. To visualize the progress of your run,  open your browser in ```http://localhost:8787``` and make sure you connect by ```ssh -L 8787:localhost:8787 username@servername``` if you are not running `pyHiM` locally.

```-C or --cmd``` is an optional argument that can be used to run a specific set of functions detailed as a comma separated list. If you don't provide this argument, the full list of functions will be run and the mode of action will be determined from the ```infoList.json``` configuration file.