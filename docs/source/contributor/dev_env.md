# Development environment

## Source code

Code and issues are hosted on GitHub in [marcnol/pyHiM repository](https://github.com/marcnol/pyHiM).
If you want to contribute with code, checkout: [development process page](./dev_process.md).

## Dependencies

pyHiM relies on [apiFISH](https://github.com/apiFISH/apiFISH) for the 3D localization of barcodes. In future, more `pyHiM` common functions will be progressively moved to apiFISH to provide a commnon API for FISH-based analysis. If you develop pyHiM, you will need to install apiFISH using the following guide: [how work with apiFISH](./work_with_apifish.md).

Main pyHiM dependencies:
- [Astropy](https://www.astropy.org/)
- [scikit-image](https://scikit-image.org/)
- [starDist](https://github.com/stardist/stardist)

Dependencies needed to locally build the `readthedocs` documentation website:
- [Sphinx](https://www.sphinx-doc.org/en/master/)
- [Read the Docs](https://readthedocs.org/)

## Development editor

The development editor choice is up to you, we list here two editors suitable for pyHiM development in Python:
- [Spyder](https://www.spyder-ide.org/) is popular in scientist community
- [VSCode](https://code.visualstudio.com/) is popular in developer community

## Analysis tools

- [`black`](https://pypi.org/project/black/): Linter to apply an auto format for [PEP8](https://www.python.org/dev/peps/pep-0008/).
- [`pylint`](https://pypi.org/project/pylint/): A Python static code analysis tool.
