#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:18:17 2020

@author: marcnol
"""

from setuptools import setup, find_packages

setup(
    name='pyHiM',
    version='0.1.0',
    description='pipeline and functions to analyze Hi-M adata',
    license='MIT',
    packages=find_packages(),
    author='Marcelo Nollmann',
    author_email='marcelo.nollmann@cbs.cnrs.fr',
    keywords=['example'],
    python_requires='>=3.6.0',
    url='https://github.com/marcnol/pyHiM'
)

#python setup.py sdist
