## Environment Installation
Detailed environment installation can be found in the Jupyter Notebook Tutorial
```
!pip install --quiet scvi-colab
from scvi_colab import install
install()
!pip install --quiet git+https://github.com/yoseflab/destvi_utils.git@main
import pyreadr
!pip install scikit-misc --force
```

## Data Preparation
```
├── datasets
│   ├── data4DestVI
│       ├── scRNA.h5ad
│       └── ST_10000.h5ad

```

## Run on SeqFISH+ Data
see Jupyter Notebook


## Source:
[DestVI Official Tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/DestVI_tutorial.html)