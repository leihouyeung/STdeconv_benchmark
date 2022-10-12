# SpatialDecon
## Environment Setup
Create an R 4.1 conda environment
```
conda create -n RCTD r-base=4.1
```
Install Bioconductor and use Bioconductor to install SpatialDecon
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(version="release")

BiocManager::install("SpatialDecon")
```
## Run Example
```
Rscript SpatialDecon.R
```
This will run SpatialDecon on the seqFISH_10000 dataset. Output will be stored in `seqFISH_10000_Result/seqFISH_10000.SpatialDecon.norm.csv`.