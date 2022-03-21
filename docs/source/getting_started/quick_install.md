# Quick installation

*For the complete version of installation options, go to [pyHiM installation section](../user_guide/pyhim_installation.md).*

## Clone pyHiM repository

Create a folder where you want to install pyHiM and go inside to clone the repository. Standard location to do it is: ```$HOME/Repositories/pyHiM```

```bash
mkdir $HOME/Repositories
cd $HOME/Repositories
git clone https://github.com/marcnol/pyHiM.git
```

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

make sure you change ```.../Repositories/...``` with your directory name if this is not where you put pyHiM !

## Set up enviroment using conda

### Install conda

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

### Setup conda enviroment

Run this command in your terminal within the root 

```sh
conda env create -f environment.yml
```

## Install apifish

```bash
cd $HOME/Repositories
git clone https://github.com/apiFISH/apiFISH.git
cd apiFISH && git checkout development
```

Update `PYTHONPATH` env variable by adding the following line to your local ~/.bashrc

```sh
export PYTHONPATH="$PYTHONPATH:$HOME/Repositories/apiFISH"
```

**TODO**
- Test and validate instructions