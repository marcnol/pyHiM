#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:18:17 2020

@author: marcnol
"""

from setuptools import setup, find_packages
from datetime import datetime

with open("README.md", "r") as fh:
    long_description = fh.read()

version = "0.5.0_" + datetime.now().strftime("%Y/%m/%d %H:%M:%S")

setup(
    name='pyHiM',
    version=version,
    description='pipeline and functions to analyze Hi-M adata',
    license='MIT',
    packages=find_packages(),
    author='Marcelo Nollmann',
    author_email='marcelo.nollmann@cbs.cnrs.fr',
    keywords=[
              'astropy',
              'csbdeep',
              'dask',
              'matplotlib',
              'numpy',
              'opencv-python',
              'pandas',
              'photutils',
              'roipoly',
              'scikit-learn',
              'scikit-image',
              'stardist',
              'tensorflow',
              'tqdm',
              'mrc'
              ],
    python_requires='>=3.7.2',
    install_requires=[''],
    url='https://github.com/marcnol/pyHiM'
)


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
#################### to install ######################
######################################################
# pip install pyHiM-0.3.0.tar.gz

######################################################
#################### to test #########################
######################################################
# pytest

######################################################
# to package, install in a virtual env and run tests #
######################################################
# tox
# configuration lives in tox.ini

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
