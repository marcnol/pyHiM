# pyHiM installation


## Install and configure pyHiM

### Clone pyHiM repository

1. Create a folder where you want to install pyHiM and go inside to clone the repository. Standard location to do it is: ```$HOME/Repositories/pyHiM```

```bash
mkdir $HOME/Repositories
cd $HOME/Repositories
```

2. Choose your clone method between HTTPS or SSH key:
    - HTTPS
      ```bash
      git clone https://github.com/marcnol/pyHiM.git
      ```
    - SSH
      ```bash
      git clone git@github.com:marcnol/pyHiM.git
      ```

3. Open your $HOME/.bashrc using nano

```bash
nano $HOME/.bashrc
```

4. Add the following line to the end

```sh
export PATH="$PATH:$HOME/Repositories/pyHiM/src"
export PATH="$PATH:$HOME/Repositories/pyHiM/src/fileProcessing"
export PATH="$PATH:$HOME/Repositories/pyHiM/src/postProcessing"

export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/pyHiM/src"

export MPLBACKEND=agg
```

```{note}
Make sure you change ```.../Repositories/...``` with your directory name (step 1.) if this is not where you put pyHiM !
```

### Set up enviroment using conda

#### Install conda

1. In the Downloads directory, run follow command (or download the package from [conda installation script](https://www.anaconda.com/products/individual)):
```
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
```
2. Now, run the installation:
```
bash Anaconda3-2020.02-Linux-x86_64.sh
```

3. Accept all the questions and default installation folder.

4. Then update anaconda:
```bash
bash
conda update anaconda
```

#### Setup conda enviroment

##### Automatic 

Run this command in your terminal within the root 

```sh
conda env create -f environment.yml
```

```{note}
If you get this error:
`ImportError: Dask\'s distributed scheduler is not installed.`

You solve by running `pip install dask[complete] distributed --upgrade`.
```

##### Manual 

To manually install the necessary packages using conda, run:

```sh
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly opencv-python tqdm stardist csbdeep pympler
pip install --upgrade tensorflow
```

Remember to activate the environment before running pyHiM:

```sh
conda activate pyHiM
```

#### Upgrade scikit-image to development version

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

### Install apifish module

1. Navigate where you want install apifish
```bash
cd $HOME/Repositories
```

2. Choose your clone method between HTTPS or SSH key:
    - HTTPS
      ```bash
      git clone https://github.com/apiFISH/apiFISH.git
      ```
    - SSH
      ```bash
      git clone git@github.com:apiFISH/apiFISH.git
      ```
3. Switch on `development` branch
```bash
cd apiFISH && git checkout development
```

4. Update `PYTHONPATH` env variable by adding the following line to your local ~/.bashrc

```sh
export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/apiFISH"
```

### Additional installation to generate documentation

```sh
conda install sphinx
conda install -c conda-forge myst-parser
conda install -c conda-forge sphinxcontrib-mermaid
conda install -c conda-forge sphinx-panels
conda install -c conda-forge sphinx_rtd_theme
```
Update `PYTHONPATH` env variable, for fileProcessing scripts documentation, by adding the following line to your local ~/.bashrc

```sh
export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/pyHiM/src/fileProcessing"
```


### Script installation for super-computer centers (e.g. Meso-LR)

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

# api-fish
cd $HOME/Repositories
git clone git@github.com:apiFISH/apiFISH.git
cd apifish && git checkout development
echo 'export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/apiFISH"'  >> ~/.bashrc

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



### Test the installation

If you want to go for a test run, do the following:

- Download a test dataset from: [test dataset](https://zenodo.org/record/6351755)

- Open a terminal and cd into the folder with this dataset, for instance:

  ```sh
  cd /mnt/grey/DATA/users/marcnol/test_HiM/testDataset
  ```

- Copy a model `infoList.json` file into the data folder.
- Run `pyHiM` by executing the following command:

```bash
pyHiM.py -F .
```

If you want to run it in your data directory, then copy the configuration files to the directory where you want to run pyHiM

```bash
cp /home/rata/Repositories/pyHiM/modelParameterFiles_JSON/infoList*json path-to-your-directory
```
