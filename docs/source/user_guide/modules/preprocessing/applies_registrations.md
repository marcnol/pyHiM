# AppliesRegistrations
*Applies registration to DAPI and barcodes*

## Invoke
Inside the folder with your input data, run:
```shell
pyhim -C appliesRegistrations
```

## Inputs

|Name shape|Quantity|Mandatory|Description|
|---|---|---|---|
|infoList.json|1|Yes|Parameter file.|
|alignImages.ecsv|1|Yes|XY alignment resulting from the XY alignment produced while running `alignImages`.|

## Outputs
|Name shape|Quantity|Description|
|---|---|---|
|2D_registered.npy|2..n|2D registered images|

## Relevant options

## Description
To run this function exclusively, run *pyHiM* using the ``` -C appliesRegistrations ``` argument. It loads masks, RNA, and barcodes 2D projected images, and applies registrations to them. The resulting images are saved as npy arrays in the ```alignImages``` folder.  


## Graph

```{mermaid}
flowchart TD

		subgraph A[INPUT]
			A1([2D_DAPI])
			A2([2D_RNA])
			A3([2D_Barcodes])
			A4([dictShifts])
		end
		subgraph C[OUTPUT]
			C1([aligned_2D_DAPI])
			C2([aligned_2D_RNA])
			C3([aligned_2D_Barcodes])
		end
		A --> B1
		B1[["appliesRegistrations2fileName()"]] --> C
	
```