

![](README.assets/website_illustration.png)

# pyHiM

pyHiM is a software package developed by the [Nollmann Lab](http://www.nollmannlab.org) at the [Center of Structural Biology](http://www.cbs.cnrs.fr), a department of the [CNRS](http://www.cnrs.fr) and the [INSERM](http://www.inserm.fr). 
pyHiM implements the analysis of multiplexed DNA-FISH data, as described in our [protocols paper](https://github.com/NollmannLab/HiM_protocol).

pyHiM was entirely written in python and makes extensive use of the following packages:

- [Astropy](https://www.astropy.org/)
- [scikit-image](https://scikit-image.org/)
- [starDist](https://github.com/stardist/stardist)

pyHiM also uses functions from [Big-FISH](https://github.com/fish-quant/big-fish) to perform Gaussian 3D fits.

## Documentation

To install, please follow the tutorial [here](docs/Installing_pyHiM.md).

After you installed pyHiM, you may want to consult a [guide](docs/Running_pyHiM.md) on how to use it. More on single cell analysis [here](docs/SingleCells.md).

If you are a developer, follow packaging instructions [here](docs/packaging_pyHiM.md).


## Publications

For more information on Hi-M, please see the following resources:
- [Hi-M protocol](https://github.com/NollmannLab/HiM_protocol)
- [Hi-M method](https://www.cell.com/molecular-cell/fulltext/S1097-2765(19)30011-5)
- [A recent Hi-M application using pyHiM](https://www.nature.com/articles/s41588-021-00816-z)

## Support

If you have any question relative to the repository, please open an issue. You can also contact Marcelo Nollmann.

## License

pyHiM is licensed under GPLv3 (see LICENSE.txt).
Packages used by pyHiM are licensed under the revised 3-clause BSD style license.

Check COPYRIGHT.txt for a list of authors and the git history for their individual contributions.


