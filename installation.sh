#!/usr/bin/bash3
# Created on Tue May  4 09:23:56 2021
# @author: marcnol
# pyHiM installation script
#!/bin/bash


# for use in muse-LR:
# load conda otherwise install before running script
# module load  python/Anaconda/3-5.1.0

# create environment and install packages
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly opencv-python tqdm stardist csbdeep
pip install --upgrade tensorflow

# big-fish
cd $HOME/Repositories
git clone https://github.com/fish-quant/big-fish.git
cd big-fish && git checkout develop
ln -s $HOME/Repositories/big-fish/bigfish $HOME/Repositories/pyHiM/src/bigfish

# clone pyHiM
# make sure you copied SSH keys to GITHUB: https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
cd $HOME/Repositories
git clone git@github.com:marcnol/pyHiM.git
git checkout development

# settings
ln -s $HOME/Repositories/pyHiM/src/fileProcessing/cleanHiM_run.py $HOME/bin/cleanHiM


