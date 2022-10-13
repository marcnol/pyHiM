# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|Yes|No| 

## Install conda

*We use conda environment to be sure that we don't have any version problem for the software dependencies with other applications.*

We recommand to download the light version `miniconda` if you just want to use conda environment. But if you are developer, you can download `anaconda` (it's the full application including [spyder IDE](https://www.spyder-ide.org/)).

- [Installing conda on Linux](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)
- [Installing conda on Windows](https://conda.io/projects/conda/en/latest/user-guide/install/windows.html)

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



