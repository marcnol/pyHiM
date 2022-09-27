# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|No|No| 

## Install conda

Download anaconda following the steps in the official [Anaconda website](https://www.anaconda.com/products/distribution) (**or** the light version [miniconda](https://docs.conda.io/en/latest/miniconda.html)).

```{note}
To be sure that we don't have any version problems for the software dependencies with other applications, we use conda environment.
```

## Create conda enviroment

Create a conda environment and activate it:
```
conda create -n pyHiM python=3.9
conda activate pyHiM
```

## Install pyHiM

```bash
pip install pyhim
```

```{note}
To check if pyHiM is well installed, run:
`pyhim --help`
```

## For users of jupyter labs

If you want to use *pyHiM* from jupyter labs, we recommend you also run the following commands:

```
conda install ipykernel matplotlib
ipython kernel install --user --name=pyHiM-kernel
```

Once you spin up a jupyter lab from the `base` environment, select the `pyHiM-kernel` to be able to run *pyHiM* functions.




