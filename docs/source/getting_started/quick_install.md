# Quick installation

*For the complete description of installation options, please refer to [pyHiM installation section](../user_guide/pyhim_installation.md).*

## Install using conda and pip

### Install conda and create enviroment

Download anaconda following the steps in the official [Anaconda website](https://www.anaconda.com/products/distribution). Briefly, for linux this requires you to download the installation script by running the following command on a terminal

```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```

### Setup conda enviroment

Create a conda environment by running 
```bash
conda create -n pyHiM
```

### Install pyHiM

Now activate the environment and install pyHiM:

```bash
conda activate pyHiM
pip install pyhim
```

### Using jupyter labs

If you want to use pyHiM from jupyter labs, we recommend you also run the following commands:

```
conda install ipykernel matplotlib
ipython kernel install --user --name=pyHiM-kernel
```

Once you spin up a jupyter lab from the `base` environment, select the `pyHiM-kernel` to be able to run pyHiM functions.




