# Tutorial for pyHiM notebook

*A jupyter notebook is an interactive file where you can find both markdown text and executable code with its outputs dispayed.*

To run by yourself the [pyHiM tutorial](notebooks/full_pyHiM_run.ipynb), follow these steps:

## Install and configurate JupyterLab

1. Activate your [conda environment](../quick_install.md#create-conda-enviroment) for pyHiM:
```sh
conda activate pyHiM
```

```{note}
Replace `pyHiM` by your `environment name` if you called it something else.
```

2. Install a tool to manage jupyter notebook like [JupyterLab](https://jupyter.org/install#jupyterlab):
```sh
conda install jupyterlab
```

3. We recommend to create a specific kernel to run pyHiM on JupyterLab with the good environment:

```
conda install ipykernel matplotlib
ipython kernel install --user --name=pyHiM-kernel
```

## Open tutorial with JupyterLab

1. To download and unzip the pyHiM notebook with its python file, [click here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/marcnol/pyHiM/tree/development/docs/source/getting_started/tutorials/notebooks).

2. Open a terminal inside your downloaded folder and activate your [conda environment](../quick_install.md#create-conda-enviroment) for pyHiM
```sh
conda activate pyHiM
```

3. Open pyHiM tutorial with JupiterLab (or jupyter notebook):
```sh
jupyter-lab full_pyHiM_run.ipynb
```

4. Once you spin up a jupyter lab from the `base` environment, select the `pyHiM-kernel` (click on panel Kernel > Change Kernel...) to be able to run *pyHiM* functions.

![select_kernel_screenshot](../../_static/select_kernel.png)

5. Now you can follow the tutorial by running each cell with the `run` icon (or `Shift+Enter` on keyboard):

![run_notebook_screenshot](../../_static/run_notebook.png)