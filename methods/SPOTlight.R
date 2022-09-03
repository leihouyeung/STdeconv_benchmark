library(SPOTlight)
library(Seurat)
library(SeuratObject)

setwd('~/STreview/STdeconvolution_datasets/seqFISH+/')
org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)

sc_exp = read.table('raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F)
st_location = read.csv('Out_rect_locations.csv',header = T, row.names = 1)

colnames(sc_exp) = gsub('_','',colnames(sc_exp))
rownames(sc_exp) = gsub('_','',rownames(sc_exp))

pbmc <- CreateSeuratObject(counts = sc_exp)
Seurat::Idents(object = pbmc) <- sc_anno[,1]
pbmc = Seurat::SCTransform(pbmc, verbose = FALSE)
cluster_markers_all <- Seurat::FindAllMarkers(object = pbmc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)
pbmc$subclass = sc_anno[,1]
st_data <- CreateSeuratObject(counts = as.matrix(t(org_st_count)))
st_data = Seurat::SCTransform(st_data, verbose = FALSE)

spotlight_ls_pbmc <- spotlight_deconvolution(
  se_sc = pbmc,
  counts_spatial = st_data@assays$RNA@counts,
  clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 3, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 )# Remove those cells contributing to a spot below a certain threshold 
spotlight_pred = as.matrix(spotlight_ls_pbmc[[2]])
write.csv(spotlight_pred, 'SPOTlight_seqFISH_10000.csv')
