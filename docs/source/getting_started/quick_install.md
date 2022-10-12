# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|Yes|No| 

## Install conda

Download anaconda following the steps in the official [Anaconda website](https://www.anaconda.com/products/distribution) (**or** the light version [miniconda](https://docs.conda.io/en/latest/miniconda.html)).

```{note}
It's not mandatory. But we use conda environment to be sure that we don't have any version problems for the software dependencies with other applications.
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



