help_dic = {
    "DAPI_channel": "Indicate the acquisition channel of the DAPI cycle ch00, ch01 or ch02\n\nDefault Value = ch00",
    "fiducialDAPI_channel": "Indicate the acquisition channel of the fiducial DAPI/RNA cycle ch00, ch01 or \
    ch02\n\nDefault Value = ch01",
    "barcode_channel": "Indicate the acquisition channel of the barcode cycle ch00 or ch01\n\nDefault Value = ch01",
    "fiducialBarcode_channel": "Indicate the acquisition channel of the fiducial barcode cycle ch00 or \
    ch01\n\nDefault Value = ch00",
    "mask_channel": "Indicate the acquisition channel of the mask cycle ch00, ch01\n\nDefault Value = ch01",
    "fiducialMask_channel": "Indicate the acquisition channel of the fiducial mask cycle ch00 or \
    ch01\n\nDefault Value = ch00",
    "RNA_channel": "Indicate the acquisition channel of the RNA cycle ch00, ch01 or ch02\n\nDefault Value = ch01",
    "fiducialRNA_channel": "Indicate the acquisition channel of the fiducial DAPI/RNA cycle ch00, ch01 or ch02",
    "pixelSizeXY": "Indicate the XY pixel size in um\n\nDefault Value = 0.1",
    "pixelSizeZ": "Indicate the Z pixel size in um\n\nDefault Value = 0.25",
    "blockSize": "Define the block size for 3D local alignment, this number indicate how the original should be \
    partitioned for the local drift correction\n\nDefault Value = 256",
    "referenceFiducial": "Define the name of the reference fiducial cycle which corresponds to the first cycle \
    acquired during the experience\n\nDefault Value = RT1",
    "flux_min": "Set minimum flux per spot for 2D. If the flux is smaller than the threshold the localization will be \
    discarded\n\nDefault Value = 10",
    "flux_min_3D": "Set minimum flux per spot for 3D. If the flux is smaller than the threshold the localization will \
    be discarded\n\nDefault Value = 4000",
    "toleranceDrift": "Set tolerance used for block drift correction (in pixels), if the drift correction if higher \
    than this threshold the localization will be discarded\n\nDefault Value = 1",
    "mask_expansion": "Set the expansion in pixel that should be used for to dilate the masks\n\nDefault Value = 8",
    "folder": "Give a name of the output folder\n\nDefault Value = buildsPWDmatrix",
    "masks2process": "Set the list of masks that need to be segmented\n\nDefault Value = nuclei:DAPI, mask1:mask1",
    "tracing_method": "Set list of methods it will use : masking,clustering\n\nDefault Value = masking,clustering",
    "KDtree_distance_threshold_mum": "Set distance threshold used to build KDtree\n\nDefault Value = 1",
    "Stardist_basename": "Set path name of AI models for 2D segmentation with StarDist",
    "brightest": "Set max number of objects segmented per FOV (only for barcodes!); “Number of brightest objects to \
    keep after sorting the full object list. If brightest is set to None, all objects will be selected” \
    (photutils.DAOStarFinder)\n\nDefault Value = 1100",
    "segmentObject_Labels_dapi_aeraMax": "Maximum size of the segmented DAPI mask\n\nDefault Value = 3000",
    "segmentObject_Labels_dapi_aeraMin": "Minimum size of the segmented DAPI mask\n\nDefault Value = 150",
    "zProject_dapi_zmax": "Define the zmax plane to use for the zprojection of the DAPI\n\nDefault Value = 59",
    "zProject_dapi_zmin": "Define the zmin plane to use for the zprojection of the DAPI\n\nDefault Value = 1",
    "zProject_Bcd_zmax": "Define the zmax plane to use for the zprojection of the barcodes\n\nDefault Value = 59",
    "zProject_Bcd_zmin": "Define the zmin plane to use for the zprojection of the barcodes\n\nDefault Value = 1",
    "zProject_Mask_zmax": "Define the zmax plane to use for the zprojection of the mask\n\nDefault Value = 59",
    "zProject_Mask_zmin": "Define the zmax plane to use for the zprojection of the mask\n\nDefault Value = 1",
}
