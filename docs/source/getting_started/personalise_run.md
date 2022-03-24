# WIP - How-to personalise a pyHiM runtime ?

## Separate input and output data
This last step produce output data in the same folder than your input data. If you want to separate the both, run pyhiM from an output directory with this command:
```bash
pyHiM.py -F <input_directory_path>
```
Where `<input_directory_path>` is the relative or absolute path to the folder containing your input data.

## Optional arguments

If you have any doubt, you can use the help option, `pyhiM.py -h`, from anywhere to see every optional arguments and their description like this:
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



```-F ``` or ```--rootFolder``` indicates the rootFolder where pyHiM expects to find the dataset.

```--threads``` argument will make it run in parallel mode. Be ready to open your browser in ```http://localhost:8787``` and make sure you connect by ```ssh -L 8787:localhost:8787 username@servername```. Think to change your username and server name.

```-C or --cmd``` is an optional argument that can be used to run a specific set of functions detailed as a comma separated list. If you don't provide this argument, the full list of functions will be run and the mode of action will be determined from the ```infoList.json``` file (see below for details).

**TODO:**
- *a example of pyHiM run with -C option and with different parameter options.*