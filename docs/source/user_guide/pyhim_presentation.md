# WIP - What is pyHiM ?

## Main goals

In the context of structural biochemistry, pyHiM software provides the **processing** of **multiplexed FISH data** produced with HiM protocol and **visualization** tools.

Basic concept of this software is to **determine 3D position** of each fluorescent spot in FISH images. These coordinates will be associated with their corresponding cell by using masks like DAPI.

We have as input **HiM experience images** and as output **position matrices** of spots.

The main use of pyHiM is to detect a set of **DNA loci**, localized with _barcodes_. In order to produce a **distance matrix** between each _barcode_. This _barcodes_ must often be detected on different images, hence the interest of _fiducials_ (see below).

FISH techniques can also be used to observe RNA and proteins. This software provides a basic way to detect a **type of RNA** but it doesn't include processing of data FISH for proteins.

## Advantages

Each image taken during acquisition is associated with a spatial reference, the _fiducial_, common to all images. So pyHiM can :

-   **Rectifies slight shifts** from one image to another due to experimental conditions.
    
-   **Associates spots** from different images in a same area delimited by a mask.
    

Acquisition of HiM protocol experiment often generates huge volume data which makes the analysis very time consuming. So, you can run pyHiM with **parallelization** on computing servers.

## Global structure
The pyHiM software follow a pipeline pattern. A default pipeline is defined with main features but they can be used independently if you have the right input data.

Features in the default pipeline can be classified into five categories :
1. **Preprocessing** manipulate or delete a part of data before it is used in order to ensure or enhance performance.
2. **Identification** segment areas for masks and detect coordinates for spots using starDist IA models.
3. **Matching** allow to associate spots or mask with their mask.
4. **Postprocessing** make output data easier to manipulate.
5. **Visualization** communicate progress and results to the user.

Each step can be optimized with **parallel computations** using the Dask package.

The use cases of pyHiM can be summarized in this diagram below :

```{mermaid}
flowchart LR
subgraph uc[pyHiM]
	f1([Preprocessing])
	f12([Identification])
	f7([Matching])
	f17([Postprocessing])
	f9([Visualization])
end
Biologist --- f1 & f12 & f7 & f17 & f9
uc -.- a2[Parallelization servers]
f12 --- a1[IA models]
```