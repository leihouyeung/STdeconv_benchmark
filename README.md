# The benchmarking of methods in spatial transcriptomics deconvolution
## Description
We benchmarked 16 existing methods resolving a cellular deconvolution task with 44 real-world and simulated datasets by evaluating the accuracy, robustness, and usability of the methods. For accuracy, We used JSD score, RMSE and Pearson correlation coefficient to quantify the performance of all methods on the datasets from six image-based and sequencing-based spatial transcriptomics technologies. For robustness, we designed several experiments with different conditions as follows: (1) 10,000, 6000, or 3000 genes were randomly chosen in the seqFISH+ dataset; (2) three resolutions were simulated with 12 MERFISH datasets using binning sizes of 20, 50, and 100 μm, which were also used to test performance with different numbers of spots; (3) 17 original cells types and 11 integrated cell types were tested in Slide-seqV2 datasets; and (4) the stability of the performance was assessed by repeating the experiments three times with the seqFISH+ dataset using 10,000 genes per spot. For usability, we recorded the running time with three different spot numbers in the MERFISH datasets and scored several aspects of the tutorials of all methods, including document quality, code quality, installation procedure, and example analysis.


![image](./pipeline.png 'pipeline')

---

## implementation

For implementation, we followed the codes from all proposed methods with their default parameters. In the `/methods`, we showed the implementation of all 16 methods on seqFISH+ dataset as an example. Please make sure the installation of all methods following these versions: [CARD (v1.0.0)](https://github.com/YingMa0107/CARD),[SPOTlight (v0.99.0)](https://github.com/MarcElosua/SPOTlight), [DSTG](https://github.com/Su-informatics-lab/DSTG), [SpatialDWLS](https://github.com/RubD/Giotto/), [SD<sup>2</sup>](https://github.com/leihouyeung/SD2), [NMFreg](https://github.com/tudaga/NMFreg_tutorial), 
[stereoscope (v.03)](https://github.com/almaan/stereoscope), [Tangram (v1.0.3)](https://github.com/broadinstitute/Tangram), [Cell2location (v0.1)](https://github.com/BayraktarLab/cell2location), [STdeconvolve (1.0.0)](https://github.com/JEFworks-Lab/STdeconvolve), [Berglund (0.2.0)](https://github.com/SpatialTranscriptomicsResearch/std-poisson), [SpiceMix](https://github.com/ma-compbio/SpiceMIx), [RCTD (spacexr 2.0.0)]( https://github.com/dmcable/spacexr), [STRIDE](https://github.com/DongqingSun96/STRIDE), [DestVI (scvi-tools 0.16.0)](https://github.com/scverse/scvi-tools) and [SpatialDecon](https://github.com/Nanostring-Biostats/SpatialDecon.git).

In the `/evaluation & visualization`, we implemented the visualization and quantification of all figures in our manuscript and supplementary materials. 
