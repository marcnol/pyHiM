# Introduction

![A *pyHiM* output example](../_static/welcome_illustration.png)

In the context of structural biochemistry, *pyHiM* software provides the tools for **processing multiplexed DNA-FISH data** produced **by H-iM experiments** as well as the **visualization** tools to explore the results.

```{note}
For more information on Hi-M, please see the following resources:

- [Hi-M protocol](https://github.com/NollmannLab/HiM_protocol)
- [Hi-M method](https://www.cell.com/molecular-cell/fulltext/S1097-2765(19)30011-5)
- [A recent Hi-M application of Hi-M](https://www.nature.com/articles/s41588-021-00816-z)
```

Basic concept of this software is to **calculate the 3D positions** of a set of DNA-loci, referred to as barcodes. Data acquisition is performed sequentially, in a series of **cycles** combining the acquisition of a single barcode with a **fiducial marker** common to all cycles and later used for drift correction. Each set of barcode coordinates will be associated with its corresponding cell by using masks (segmented nuclei stained with DAPI for instance). At the end of the analysis, *pyHiM* will output several files describing the composition (barcodes), localization and 3D-coordinates of all individual traces detected.

*Optionally, pyHiM also provides a basic way to detect a **single type of RNA** but it doesn’t include processing for multiple RNA species or data for proteins-FISH.*


**A minimal HiM set of data contains:**

- images of nuclei for a given Region Of Interest (ROI)
- images of DNA-FISH barcodes acquired in different hybridization cycles for the same ROI as for the nuclei.
- fiducial images for both nuclei and barcodes acquired simultaneously with the images above (typically in a different color).



***pyHiM* has specific modules that allow the user to:**

- register nuclei and barcode images to remove or minimize drift.
- segment nuclei
- segment and localize barcodes
- construct chromatin traces
- build single-trace and ensemble pair-wise distance maps
