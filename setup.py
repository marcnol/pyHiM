#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:18:17 2020

@author: marcnol
"""

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='pyHiM',
    version='0.1.0',
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
              'tensorflow'
              ],
    python_requires='>=3.6.0',
    install_requires=[''],
    url='https://github.com/marcnol/pyHiM'
)


# to package
# python3 setup.py sdist bdist_wheel

# to install
# pip install pyHiM-0.1.0.tar.gz


# list of packages needed

# conda install pandas scikit-image numpy matplotlib astropy graphviz
# conda install photutils -c astropy

# pip install roipoly opencv-python tqdm stardist csbdeep snakemake 
# pip install --upgrade tensorflow

