**Cell2location**:

[Installation]: 

The pipeline implemented here, is a direct application of the original  pipeline provided by the cell2location paper (https://github.com/BayraktarLab/cell2location). Thus it requires nearly the same procedures to run. 

A seperate conda environment for installing cell2location is suggested.

First, a conda environment is created, and cell2location package is installed.

```{python}
conda create -y -n cell2loc_env python=3.9

conda activate cell2loc_env
pip install git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location[tutorials]
```

Then add jupyter kernel for this environment in the notebook to use it

```{python}
conda activate cell2loc_env
python -m ipykernel install --user --name=cell2loc_env --display-name='Environment (cell2loc_env)'
```

Note that before installing cell2location and it's dependencies, it could be necessary to make sure that you are creating a fully isolated conda environment by telling python to **NOT** use user site for installing packages by running this line before creating conda environment and every time before activating conda environment in a new terminal session.

```{python}
export PYTHONNOUSERSITE="literallyanyletters"
```

[Run]:
The .ipynb file contains the codes necessary to run the example on seqFISH_10000. The input files and output files may require changes, if you wish to run your own example. 


