**SPICEMIX**:

[Installation]: 

Please download from https://github.com/ma-compbio/SpiceMix.

The environment of the project requires the following dependencies: `Python`, `scipy`, `gurobi=8.1.1`, `pytorch=1.4.0`, `numpy`, `scikit-learn`, `h5py`. The authors of the original SPICEMIX paper provides a environment file for the installation in the Anaconda environment. We can simply run, if we are using anaconda,

```{python}
conda env create -f SpiceMix.yml
```

Note that the implementation is based on an older version of pytorch, some APIs may not be compatiable with the newest version of pytorch. Gurobi is a commercial and outstanding package. It is used to solve quadratic programmings in SpiceMix. To use it, you will need to have a license, which can be requested on https://www.gurobi.com/.

The computing resource needed here includes a CPU, and 8GB RAM. A GPU is optional, but if you plan to use one, please install GPU versions of pytorch.

[Run]:
The detailed Running instructions on running and all parameters can be found at https://github.com/ma-compbio/SpiceMix. However, in our case, to run the seqFISH_10000 experiments, we simply need to run

```{python}
python main.py -K=6 --path2dataset="../datasets/seqFISH+_10000" --repli_list="[1]" --use_spatial="[True]" --result_filename="SPICEMIX"
```

Note that the inputs and outputs filenames here may require minor changes 
if you want to run your own examples.