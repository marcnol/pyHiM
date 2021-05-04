#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:07:12 2020

@author: marcnol
"""

import toml

with open("test.toml", "r") as myfile:
    dataString = myfile.readlines()

parsed_toml = toml.loads("".join(dataString))
