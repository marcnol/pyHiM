#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 11:09:04 2021

@author: marcnol
"""

import logging

formatter1 = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
formatter2 = logging.Formatter("%(message)s")

logFile = "/home/marcnol/data/Embryo_debug_dataset/test_dataset/test.log"

logger = logging.getLogger()  # root logger - Good to get it only once.
logger.handlers = []
for hdlr in logger.handlers[:]:  # remove the existing file handlers
    if isinstance(hdlr,logging.FileHandler):
        logger.removeHandler(hdlr)

filehandler = logging.FileHandler(logFile, 'w')
ch = logging.StreamHandler()

filehandler.setLevel(logging.INFO)
# ch.setLevel(logging.WARNING)
logger.setLevel(logging.INFO)

logger.addHandler(ch)
logger.addHandler(filehandler)

filehandler.setFormatter(formatter1)
ch.setFormatter(formatter2)
