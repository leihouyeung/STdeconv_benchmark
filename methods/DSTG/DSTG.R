library(Matrix)
library(data.table)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggsci)
library(Matrix)
library(rdist)
library(rBeta2009)
library(RANN)

source('ours_utils.R')
current_path = getwd()
setwd('~/STreview/STdeconvolution_datasets/seqFISH+/')
org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)
sc_exp = read.table('raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F)

setwd(current_path)
st_count = t(as.matrix(org_st_count))
cell_type = sc_anno[,1]

process_data(as.matrix(sc_exp),st_count,cell_type,spot_num = 300, 
             dropout_extract = F,combine_feature = F,HVG = T,
             HVG_num = 1000,scale_num = 10000,lower_cellnum = 10,upper_cellnum = 20)

system('python train.py ')
