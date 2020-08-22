import glob,os
rootFolder='rawImages'
tag='.tif'
TIFS=[x.split(tag)[0] for x in glob.glob(rootFolder+os.sep+"*.tif") if tag in x]

rule all:
	input: ["rawImages/zProject/{}_2d.npy".format(os.path.basename(x).split(".")[0]) for x in TIFS]
	    	
rule zProject:
    input:
    	tifs=["{}.tif".format(x.split(".")[0]) for x in TIFS]
    output:
        ["rawImages/zProject/{}_2d.npy".format(os.path.basename(x).split(".")[0]) for x in TIFS]
    shell:
        "runMakeProjections.py -F rawImages --parallel -x {input.tifs}"

rule clean:
    input:
    	rootFolder="{}".format(rootFolder)
    shell:
        "cleanHiM_run.py -F {input.rootFolder}"

