#!/bin/bash

echo ">>> Initialising Conda..."
source ~/miniconda3/etc/profile.d/conda.sh

echo ">>> Creating fresh environment..."
conda create -n test_build_doc python=3.9 --yes
conda activate test_build_doc

echo ">>> Installing pandoc (simulating RTD docker image)..."
conda install -c conda-forge pandoc -y

echo ">>> Installing the package without dependency..."
pip install --no-deps .

echo ">>> Installing doc dependencies..."
pip install -r docs/requirements.txt

echo ">>> Cleaning previous build output..."
rm -r docs/build/html

echo ">>> Building documentation with Sphinx..."
sphinx-build -b html docs/source docs/build/html

echo ">>> Opening documentation in browser..."
firefox docs/build/html/index.html

echo ">>> Clean environment..."
conda deactivate
conda env remove --name test_build_doc
