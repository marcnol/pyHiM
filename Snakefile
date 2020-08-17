import glob,os
rootFolder='.'
tag='_ch01.tif'
SAMPLES=[x.split(tag)[0] for x in glob.glob(rootFolder+os.sep+"*.tif") if tag in x]

rule zProject:
    input:
    	["{}_ch01.tif".format(x.split(",")[0]) for x in SAMPLES] 
    output:
        ["zProject/{}_ch01_2d.npy".format(x.split(",")[0]) for x in SAMPLES]
    shell:
        "runMakeProjections.py -F . -x {input}"

