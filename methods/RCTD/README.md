# RCTD (spacexr)
## Environment Setup
Create an R 3.6 conda environment and install `devtools`
```
conda create -n RCTD r-base=3.6 r-devtools
```
Use `devtools` to install spacexr
```
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```
## Run Example
```
Rscript RCTD.R
```
This will run RCTD (spacexr) on the seqFISH_10000 dataset. Output will be stored in `seqFISH_10000_Result/seqFISH_10000.RCTD.norm.csv`.