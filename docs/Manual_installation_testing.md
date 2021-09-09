# Installing and testing pyHiM

## Install pyHiM using pip version

```sh
pip install pyHiM-0.3.0.tar.gz
```

## Test pyHiM 

```sh
pytest
```

## Manual installation of pyHiM

```sh
conda create --name pyHiM python=3.7.2 dask numpy matplotlib astropy scikit-learn pandas
conda activate pyHiM
conda install photutils -c astropy
pip install mrc roipoly opencv-python tqdm stardist csbdeep
pip install --upgrade tensorflow
```

Install bigfish

```sh
cd ~/Repositories
git clone https://github.com/fish-quant/big-fish.git
cd big-fish && git checkout develop
ln -s $HOME/Repositories/big-fish/bigfish $HOME/Repositories/pyHiM/src/bigfish
```


for an installation without environment:

```
ln -s $HOME/Repositories/big-fish $HOME/anaconda3/lib/python3.7/bigfish
```

sometimes, you will need to upgrade packages to their last default version. For this, run:

```sh
conda update -c astropy astropy
conda update photoutils -c astropy
pip install update stardist
pip install --upgrade pip
pip install --upgrade dask
conda install spyder=4.2.0
```

