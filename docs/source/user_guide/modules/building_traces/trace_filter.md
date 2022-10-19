# trace_filter

This script is intended to filter the traces with a low number of spots. It acts on a single or several trace files and takes as input the minimum number of spots a trace needs to be kept.

How to use it:

- Example 1: Invoke from `rootFolder`

  ```sh
  # Call trace_filter from rootFolder to process all Trace files in `buildPWDmatrix`
  $ trace_filter.py -- N_barcodes 2
  ```

  this will process all  `Trace` files in `buildPWDmatrix`

  

- Example 2: using cat or ls to provide a specific list of Trace files to process

  ```sh
  # either make list of files to process and write it in a file
  $ cat files_to_combine 
  folder1/Trace_3D_barcode_KDtree_ROI:15.ecsv
  folder2/Trace_3D_barcode_KDtree_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ cat files_to_combine | trace_combinator.py --N_barcodes 2 --pipe
  
  # OR use `ls` to select which files you want to combine
  $ ls folder1/Trace*3D*6.ecsv
  folder1/Trace_3D_barcode_KDtree_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:DAPI_ROI:6.ecsv  folder1/Trace_3D_barcode_mask:mask0_ROI:6.ecsv
  
  # then pipe these files into trace_combinator
  $ ls folder1/Trace*3D*6.ecsv | trace_combinator.py --N_barcodes 2 --pipe
  ```




**Relevant options**

```sh
usage: trace_filter [-h] [-F ROOTFOLDER] [--N_barcodes N_BARCODES] [--pipe]

optional arguments:
  -h, --help            show this help message and exit
  -F ROOTFOLDER, --rootFolder ROOTFOLDER
                        Folder with images
  --N_barcodes N_BARCODES
                        minimum_number_barcodes. Default = 2
  --pipe                inputs Trace file list from stdin (pipe)
```
