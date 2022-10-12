**STdeconvolve**:

[Installation]: 

To install STdeconvolve, using `remotes` is recommended.

```{R}
require(remotes)
remotes::install_github('JEFworks-Lab/STdeconvolve')
```

STdeconvolve is also available through Bioconductor. Note that through Bioconductor (release 3.15), the R version must be >=4.2.

```{R}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("STdeconvolve")
```

[Run]:
The .Rmd file contains the codes necessary to run the example on seqFISH_10000. The input files and output files may require changes, if you wish to run your own example.