#!/usr/bin/bash
# -*- coding: utf-8 -*-
#Created on Thu Jun  4 14:18:18 2020
#
#@author: marcnol

############## Folders
FIGS=/home/marcnol/Repositories/Paper_doc_HiM/PDFs
#DATA1=/home/marcnol/data/Experiment_18
DATA1=/home/marcnol/data/fullDatasets/docwt_hiRes_nc14
DATA2=/home/marcnol/data/Experiment_18
DATA4=/run/media/marcnol/eeed4b81-e9df-4e64-ae90-3298d2650005/data_doc/Experiment_5

############## Figure 1
# Figure 1C
#-----------
rm "$DATA1"/scHiMmatrices/*svg
figureHiMmatrix.py -F "$DATA1" --axisLabel --fontsize 22 --scalingParameter 1

# Figure 1D
#-----------
figure3wayInteractions.py -F "$DATA1" --label doc --action all --fontsize 12 --scalingParameter 1.0 
cp "$DATA1"/scHiMmatrices/*svg "$FIGS"/Figure1/


############## Figure 2
# Figure 2C
#-----------
rm "$DATA2"/scHiMmatrices/*svg
figureHiMmatrix.py -F "$DATA2" --label doc --action labeled --axisLabel --fontsize 15 --scalingParameter 1
figureHiMmatrix.py -F "$DATA2" --label doc --action unlabeled --axisLabel --fontsize 15 --scalingParameter 1

# Figure 2F
#-----------
figure4Mmatrix.py  --rootFolder1 "$DATA2" --rootFolder2 "$DATA2" --label1 doc --label2 doc --action1 labeled --action2 unlabeled --cAxis 1

cp "$DATA2"/scHiMmatrices/*svg "$FIGS"/Figure2/

# Figure 2G
#-----------
figure3wayInteractions.py -F "$DATA2" --label doc --action labeled --fontsize 12 --scalingParameter 1.0 
figure3wayInteractions.py -F "$DATA2" --label doc --action unlabeled --fontsize 12 --scalingParameter 1.0 
cp "$DATA2"/scHiMmatrices/*svg "$FIGS"/Figure2/

############## Figure 4
# Figure 4A
#-----------
rm "$DATA4"/scHiMmatrices/*svg
figureHiMmatrix.py -F "$DATA4" --axisLabel --fontsize 12 --scalingParameter 1


# Figure 4B
#-----------
figureCompare2Matrices.py -F1 "$DATA1" -F2 "$DATA4" --action1 all --action2 all --cAxis 1 --outputFolder "$FIGS"/Figure4

# Figure 4C
#-----------
figure3wayInteractions.py -F "$DATA4" --label doc --action all --fontsize 12 --scalingParameter 1.0 

# Figure S4A
#-----------
figure4Mmatrix.py --rootFolder1 "$DATA1" --rootFolder2 "$DATA4" --label1 doc --label2 doc --action1 all --action2 all --cAxis 1 --outputFolder "$FIGS"/Figure4

cp "$DATA4"/scHiMmatrices/*svg "$FIGS"/Figure4/
