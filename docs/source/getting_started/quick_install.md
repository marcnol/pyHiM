# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|Beta|No| 

## Install conda and create environment

Download anaconda following the steps in the official [Anaconda website](https://www.anaconda.com/products/distribution).

## Setup conda enviroment

Create a conda environment and activate it:
```bash
conda create -n pyHiM
conda activate pyHiM
```

## Install pyHiM

```bash
pip install pyhim
```

## Using jupyter labs

If you want to use *pyHiM* from jupyter labs, we recommend you also run the following commands:

```
conda install ipykernel matplotlib
ipython kernel install --user --name=pyHiM-kernel
```

Once you spin up a jupyter lab from the `base` environment, select the `pyHiM-kernel` to be able to run *pyHiM* functions.




