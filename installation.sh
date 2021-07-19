#!/usr/bin/bash3
# Created on Tue May  4 09:23:56 2021
# @author: marcnol
# pyHiM installation script
#!/bin/bash

# load conda otherwise install before running script
# module load  python/Anaconda/3-5.1.0

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
# make sure you copied SSH keys to GITHUB: https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
cd $HOME/Repositories
git clone git@github.com:marcnol/pyHiM.git
git checkout development

# settings
ln -s $HOME/Repositories/pyHiM/src/fileProcessing/cleanHiM_run.py $HOME/bin/cleanHiM

######################################################
#################### to package ######################
######################################################
# python3 setup.py sdist bdist_wheel
#

# Pack using docker:
# sudo docker build -t py_him .
#
# sudo docker save py_him >dist/docker_pyHiM.tar
# gzip dist/docker_pyHiM.tar

# Deploy docker container
# copy docker_pyHiM.tar.gz into server, then run
# docker import dist/docker_pyHiM.tar.gz new_py_him:latest
# then run by :
# docker run py_him

######################################################
################### to install pip version############
######################################################
# pip install pyHiM-0.3.0.tar.gz

######################################################
#################### to test #########################
######################################################
# pytest

######################################################
############# conventional installation ##############
######################################################
# conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
# conda activate pyHiM
# conda install photutils -c astropy
# pip install mrc roipoly opencv-python tqdm stardist csbdeep
# pip install --upgrade tensorflow

# +++ bigfish +++
# cd ~/Repositories
# git clone https://github.com/fish-quant/big-fish.git
# cd big-fish && git checkout develop
# ln -s $HOME/Repositories/big-fish/bigfish $HOME/Repositories/pyHiM/src/bigfish
#
# for an installation without environment:
#
# ln -s $HOME/Repositories/big-fish $HOME/anaconda3/lib/python3.7/bigfish

######################################################
############# upgrades                  ##############
######################################################
# conda update -c astropy astropy
# conda update photoutils -c astropy
# pip install update stardist
# pip install --upgrade pip
# pip install --upgrade dask
# conda install spyder=4.2.0

######################################################
# deprecated
# conda install -c sherpa sherpa numba
# conda install pandas
