# Installing pyHiM



## Install conda

in the Downloads directory, run:

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
```

or download the package from [conda installation script](https://www.anaconda.com/products/individual)



Now, run the installation by

```
bash Anaconda3-2020.02-Linux-x86_64.sh 

```

and accept all the questions and default installation folder. Then update anaconda by

```bash
bash
conda update anaconda
```

You are set.



## Install and configure pyHiM

Clone the repository. Standard location to do it is: ```$HOME/Repositories/pyHiM```

Open your ~/.bashrc using nano 

```bash
nano ~/.bashrc
```

and add the following line to the end

```sh
export PATH="$PATH:/home/rata/Repositories/pyHiM/"
```

make sure you use a different directory name if this is not where you put pyHiM !

Now, install the necessary packages using conda, by running:

```sh
conda install pandas scikit-image numpy matplotlib tqdm astropy
conda install photutils -c astropy

```

You should be set!



## Test run

If you want to go for a test run, do the following:

Select a folder with data to analyse, for instance: ```/mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000```



then copy the configuration files to the directory where you want to run pyHiM

```bash
cp /home/rata/Repositories/pyHiM/infoListModels/* /mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/
```

specify the RT that should be used as referenceFiducial

```bash
bash infoList_changeRT.sh
```

run by:

```bash
processingPipeline.py -F /mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/
```

## zipping and erasing run

#### zip and retrieve results

Other utilities have been written to retrieve data from a run to a remote server. For this, go to the directory with the data an run ```zipHiM_run.py```. This will archive all the png, ecsv, and dat files with the results in addition to the MD file with the output of the run. The name of the archive will be HiMrun.tar.gz.

#### clean run

If you want to erase a run, for instance to make sure you can run it again without any leftover, you can run ```cleanHiM_run.py` in the directory with the data. 

 

## Combine results from runs

Once you run a bunch of datasets, you will want to combine the PWD matrices together. For this:

1. Retrieve the matrices and barcodes files by for instance:

```bash
scp rata@lopevi:/mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/buildsPWDmatrix/*ecsv /home/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix/

scp rata@lopevi:/mnt/tronador/Sergio/RAMM_experiments/Experiment_3/deconvolved_DAPI/Embryo_000/buildsPWDmatrix/*npy /home/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix/
```

Now, you can run ```replotHiMmatrix.py``` locally and just add the data to the list as follows:

```python
ListRootFolders=[\
                 #'/mnt/disk2/marcnol/data/Experiment_3/019_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/007_Embryo/buildsPWDmatrix',\
                 '/mnt/disk2/marcnol/data/Experiment_3/016_Embryo/buildsPWDmatrix'\
                 '/mnt/disk2/marcnol/data/Experiment_3/000_Embryo/buildsPWDmatrix'\
                 ]

```

You should be able to run the following sections that will load all the datasets together and produce:

- PWD matrix for each dataset
- HistogramMatrix for each dataset
- Inverse distance matrices for each dataset
- Contact probability matrices for each dataset
- Combined contact probability matrix

This last matrix is save as a file ```'CombinedMatrix.dat'```in the ```outputFolder``` directory defined in the first block of code. It is saved as plain text so that it can be opened in MATLAB. Each row in the matrix is separated by a ```\n``` . 

For info, the barcodes used are always stored in the ```buildsPWDmatrix_uniqueBarcodes.ecsv```

file within the ```buildsPWDmatrix``` directory.









