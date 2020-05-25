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
yes
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











