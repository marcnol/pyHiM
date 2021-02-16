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

version = "0.4.0_" + datetime.now().strftime("%Y/%m/%d %H:%M:%S")

setup(
    name='pyHiM',
    version=version,
    description='pipeline and functions to analyze Hi-M adata',
    license='MIT',
    packages=find_packages(),
    author='Marcelo Nollmann',
    author_email='marcelo.nollmann@cbs.cnrs.fr',
    keywords=['roipoly',
              'scaleogram',
              'opencv-python',
              'progress',
              'astropy',
              'photutils',
              'tqdm',
              'numpy',
              'matplotlib',
              'scikit-image',
              'stardist',
              'csbdeep',
              'tensorflow',
              'dask',
              'scikit-learn',
              'numba'
              ],
    python_requires='>=3.6.0',
    install_requires=[''],
    url='https://github.com/marcnol/pyHiM'
)


######################################################
#################### to package ######################
######################################################
# python3 setup.py sdist bdist_wheel
#
# using docker
# sudo docker build -t py_him .
#
# sudo docker save py_him >dist/docker_pyHiM.tar
# sudo docker load < dist/docker_pyHiM.tar.gz
# 

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
# conda install pandas numpy matplotlib astropy
# conda install photutils -c astropy scikit-learn

# pip install roipoly opencv-python tqdm stardist csbdeep numba dask
# pip install --upgrade tensorflow

# pip install --upgrade scikit-image

######################################################
############# upgrades                  ##############
######################################################
# conda update -c astropy astropy
# conda update photoutils -c astropy
# pip install update stardist
# pip install --upgrade pip
# pip install --upgrade dask
# conda install spyder=4.2.0

