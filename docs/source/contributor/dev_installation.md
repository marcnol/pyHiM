# Developer installation

## Clone *pyHiM* repository

1. Create a folder where you want to install *pyHiM* and go inside to clone the repository. Standard location to do it is: ```$HOME/Repositories/pyHiM```

```bash
mkdir $HOME/Repositories
cd $HOME/Repositories
```

2. Choose your clone method between HTTPS or SSH key:
    - HTTPS (For latest version user)
      ```bash
      git clone https://github.com/pyHi-M/pyHiM.git
      ```
    - SSH (ONLY for developer)
      ```bash
      git clone git@github.com:pyHi-M/pyHiM.git
      ```

3. Open your $HOME/.bashrc using nano

```bash
nano $HOME/.bashrc
```

4. Add the following line to the end

```sh
export PATH="$PATH:$HOME/Repositories/pyHiM/src"
export PATH="$PATH:$HOME/Repositories/pyHiM/src/toolbox/file_handling"
export PATH="$PATH:$HOME/Repositories/pyHiM/src/postProcessing"
export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/pyHiM/src"
export MPLBACKEND=agg
```

```{note}
Make sure you change ```.../Repositories/...``` with your directory name (step 1.) if this is not where you put *pyHiM* !
```


## Installation using uv

This is much faster than conda or mamba and is now the preferred installation method.

To install using `uv` just go to the Repository folder and run the bash script:

```sh
cd $HOME/Repositories/pyHiM
bash install_pyhim_uv.sh
```

Make sure you added the environmental variables in `.bashrc'.

Everytime you run pyHiM you need to activate the environment doing:

```sh
source $HOME/Repositories/pyHiM/.venv/bin/activate
```

## Install using conda/mamba

Follow the Miniconda instructions:
[Installing miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)

Or follow these instructions to install mamba:

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/$(uname)-$(uname -m)/latest | tar -xvj bin/micromamba
```

## Automatically configure pyHiM

Run this command in your terminal:

```sh
cd $HOME/Repositories/pyHiM
conda env create -f environment.yml
```

or

```sh
cd $HOME/Repositories/pyHiM
mamba env create -f environment.yml
```

### Verify GPU recognition

After activating the conda environment, you can verify that TensorFlow detects your GPU by starting Python and running:

```python
import tensorflow as tf
print(tf.__version__)
print("Built with CUDA:", tf.test.is_built_with_cuda())
print("Built with GPU support:", tf.test.is_built_with_gpu_support())
print("Visible GPUs:", tf.config.list_physical_devices("GPU"))
```

### Selecting a GPU device

If your machine has multiple GPUs, you can select which one pyHiM uses by setting `CUDA_VISIBLE_DEVICES` before launching the application. For example, to use device `2` on a machine with three GPUs, run the following in your shell prior to starting pyHiM:

```bash
export CUDA_VISIBLE_DEVICES=2
```

This environment variable limits visibility to the specified device index, letting you choose a different GPU than the default.

```{note}
If you get this error:
`ImportError: Dask\'s distributed scheduler is not installed.`

You solve by running `pip install dask[complete] distributed --upgrade`.
```

## Install **traceratops**

If you installed `pyHiM` using `uv` you can ignore this section as the bash script already installs traceratops.

Otherwise run one of the following:

### Latest master version [stable]

```sh
conda activate pyhim39
cd $HOME/Repositories
git clone https://github.com/pyHi-M/traceratops.git
cd $HOME/Repositories/traceratops
pip install -e .
```

### Latest dev version [stable]

```sh
conda activate pyhim39
cd $HOME/Repositories
git clone git@github.com:pyHi-M/traceratops.git
cd $HOME/Repositories/traceratops
pip install -e ".[dev]"
```

## Install **apifish**

If you installed `pyHiM` using `uv` you can ignore this section as the bash script already installs traceratops.

Otherwise run one of the following:

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

---

## Test pyHiM

- The tests use the `pytest` module.
- The test resources are inside `pyhim-small-dataset`. It's a sub-module of pyHiM, so to get the dataset you need to run:
  * `git submodule update --init --recursive`
  OR
  * `git clone --recurse-submodules <HTTPS/SSH>`

- To run the tests:

### Using uv

If you used `uv` to install `pyHiM` run:

```bash
cd $HOME/Repositories/pyHiM
source .venv/bin/activate
pytest tests/ -vv
```

### Using conda

If you used `conda` to install `pyHiM` run:

```bash
cd ~Repositories/pyHiM/
conda activate pyhim39
pytest tests/ -vv
```

### Using mamba

If you used `mamba` to install `pyHiM` run:

```bash
cd ~Repositories/pyHiM/
mamba activate pyhim39
pytest tests/ -vv
```

---

## Additional installation to generate documentation

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

## Build documentation locally
Install in your conda env:
```bash
pip install nbsphinx ipython sphinx-book-theme
conda install pandoc
```
Generate documentation:
```bash
cd docs/
make html
```
A `build/html/` folder has been created with a `index.html` file inside, open it with your favorite browser.

## Script installation for super-computer centers (e.g. Meso-LR)

To access the private repository of *pyHiM*, please first create an SSH key and put it in your keyring. Follow the steps described [here](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

Then run the following automatic script:

```sh
#!/bin/bash

# load conda
module load  python/Anaconda/3-5.1.0

# create environment and install packages
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly tqdm stardist csbdeep pympler
pip install --upgrade tensorflow

# api-fish
cd $HOME/Repositories
git clone git@github.com:apiFISH/apiFISH.git
cd apifish && git checkout development
echo 'export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/apiFISH"'  >> ~/.bashrc

# clone pyHiM
cd $HOME/Repositories
git clone git@github.com:pyHi-M/pyHiM.git
git checkout development

# settings
ln -s $HOME/Repositories/pyHiM/src/toolbox/file_handling/cleanHiM_run.py $HOME/bin/cleanHiM

```

## Step to setup pre-commit in local
- environment installation
  `pip install pre-commit`

- check if it's well installed
  `pre-commit --version`

- install command of the file ".pre-commit-config.yaml" inside ".git/hooks/pre-commit"
  `pre-commit install`

- fix strange issue or warning
  `pre-commit autoupdate --repo https://github.com/pre-commit/pre-commit-hooks`

- test pre-commit without any commit
  `pre-commit run --all-files`

- Update pre-commit file
  - `pre-commit clean`
  - `pre-commit autoupdate`
  - `pre-commit install`
