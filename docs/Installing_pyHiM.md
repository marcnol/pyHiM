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
export PATH="$PATH:/home/rata/Repositories/pyHiM/:home/rata/Repositories/pyHiM/fileProcessing"
export PYTHONPATH="/home/rata/Repositories/pyHiM"
```

make sure you use a different directory name if this is not where you put pyHiM !



To install the necessary packages using conda, run:

```sh
conda install scikit-image numpy matplotlib astropy

conda install photutils -c astropy

pip install tqdm roipoly opencv-python stardist csbdeep

pip install --upgrade tensorflow

```

You should be set!



## Test run

If you want to go for a test run, do the following:

Select a folder with data to analyse, for instance: ```/mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug```

run by:

```bash
pyHiM.py -F /mnt/grey/DATA/users/marcnol/test_HiM/merfish_2019_Experiment_18_Embryo0/debug
```



If you want to run it in your data directory, then copy the configuration files to the directory where you want to run pyHiM

```bash
cp /home/rata/Repositories/pyHiM/infoListModels/infoList*json path-to-your-directory
```

