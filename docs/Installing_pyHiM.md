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

Open your $HOME/.bashrc using nano

```bash
nano $HOME/.bashrc
```

and add the following line to the end

```sh
export PATH="$PATH:$HOME/Repositories/pyHiM/src:$HOME/Repositories/pyHiM/src/fileProcessing"
export PYTHONPATH="$HOME/Repositories/pyHiM/src"

export MPLBACKEND=agg

```

make sure you use a different directory name if this is not where you put pyHiM !

To install the necessary packages using conda, run:

```sh
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly opencv-python tqdm stardist csbdeep
pip install --upgrade tensorflow
```

Remember to activate the environment before running pyHiM:

```sh
conda activate pyHiM
```



### Installing bigfish

```bash
cd $HOME/Repositories
git clone https://github.com/fish-quant/big-fish.git
cd big-fish && git checkout develop
ln -s $HOME/Repositories/big-fish/bigfish ~/anaconda3/lib/python3.7/bigfish
```

If you are running pyHiM in a conda environment, you can link bigfish as follows:

```sh
ln -s $HOME/Repositories/big-fish/bigfish $HOME/Repositories/pyHiM/src/bigfish

```





### Upgrade scikit-image to development version

Uninstall any existing installations:

```
pip uninstall scikit-image
```

or, on conda-based systems:

```
conda uninstall scikit-image
```

Now, clone scikit-image on your local computer, and install:

```
git clone https://github.com/scikit-image/scikit-image.git
cd scikit-image
pip install -e .
```

To update the installation:

```
git pull  # Grab latest source
pip install -e .  # Reinstall
```



You should be set!



### Install in Meso-LR super-computer

To access the private repository of pyHiM, please first create an SSH key and put it in your keyring. Follow the steps described [here](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

Then run the following script:

```sh
#!/bin/bash

# load conda
module load  python/Anaconda/3-5.1.0

# create environment and install packages
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly opencv-python tqdm stardist csbdeep pympler
pip install --upgrade tensorflow

# big-fish
cd $HOME/Repositories
git clone https://github.com/fish-quant/big-fish.git
cd big-fish && git checkout develop
ln -s $HOME/Repositories/big-fish/bigfish $HOME/Repositories/pyHiM/src/bigfish

# clone pyHiM
cd $HOME/Repositories
git clone git@github.com:marcnol/pyHiM.git
git checkout development

# settings
ln -s $HOME/Repositories/pyHiM/src/fileProcessing/cleanHiM_run.py $HOME/bin/cleanHiM

```



## Test run

There are two ways to do a test run in this version of pyHiM

### pytest

The first is to use the **pytest** module. For this, you need to first configure the location of the test datasets by


```sh
cd tests/standardTests
python create_testDataJSON.py -F location_of_test_dataset
```

For instance if your debug dataset is in ```/mnt/grey/DATA/users/marcnol/test_HiM/testDataset``` then running

```sh
python create_testDataJSON.py -F /mnt/grey/DATA/users/marcnol/test_HiM/testDataset
```

Now go to the pyHiM root directory (e.g. ```cd /home/marcnol/Repositories/pyHiM```) and run

```
pytest
```

### Do a mock analysis on the test dataset

If you want to go for a test run, do the following:

Select a folder with data to analyse, for instance: ```/mnt/grey/DATA/users/marcnol/test_HiM/testDataset```

run by:

```bash
pyHiM.py -F .
```

If you want to run it in your data directory, then copy the configuration files to the directory where you want to run pyHiM

```bash
cp /home/rata/Repositories/pyHiM/modelParameterFiles_JSON/infoList*json path-to-your-directory
```
