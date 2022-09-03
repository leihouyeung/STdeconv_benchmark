library(Giotto)

setwd('~/STreview/STdeconvolution_datasets/seqFISH+/')
org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)

sc_exp = read.table('raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F)
st_location = read.csv('Out_rect_locations.csv',header = T, row.names = 1)
cell_type = sc_anno[,1]

my_python_path= "~/anaconda3/envs/st3/bin/python3"
instrs = createGiottoInstructions(python_path = my_python_path)

sc_obj = createGiottoObject(raw_exprs = sc_count, instructions = instrs)
sc_obj = filterGiotto(gobject = sc_obj, expression_threshold = 0.1, gene_det_in_min_cells = 10,
                      min_det_genes_per_cell = 10, expression_values = c('raw'), verbose = T)
sc_obj = normalizeGiotto(gobject = sc_obj, scalefactor = 6000, verbose = T)
sc_obj <- addStatistics(gobject = sc_obj)


# add cell type annotation
anno = data.table::data.table(cell_ID = sc_obj@cell_ID, cell_type = cell_type)
colnames(anno) = c('cell_ID','cell_type')

sc_obj@cell_metadata = data.table::merge.data.table(sc_obj@cell_metadata, anno, by ='cell_ID')
gini_markers = findMarkers_one_vs_all(gobject = sc_obj,
                                      method = 'gini',
                                      expression_values = 'normalized',
                                      cluster_column = 'cell_type',
                                      min_genes = 20,
                                      min_expr_gini_score = 0.5,
                                      min_det_gini_score = 0.5)

sign_markers = unique(gini_markers$genes[which(gini_markers$comb_rank <= 1000)])
topgenes_gini = gini_markers[, head(.SD, 2), by = 'cluster']
average_cell_type_expr = Giotto:::create_average_DT(gobject = sc_obj, 
                                                    meta_data_name = 'cell_type', 
                                                    expression_values = 'normalized')
average_cell_type_expr = average_cell_type_expr[sign_markers,]
colnames(average_cell_type_expr) = gsub('cluster_', '', colnames(average_cell_type_expr) )

locs = data.table::data.table(cell_ID = rownames(st_location), st_location)


instrs = createGiottoInstructions(python_path = my_python_path)
data.table::setnames(locs, new = c('sdimx', 'sdimy', 'cell_ID'))

st_obj = createGiottoObject(raw_exprs = st_count, 
                            spatial_locs = locs,
                            instructions = instrs)
st_obj = filterGiotto(gobject = st_obj, expression_threshold = 1, gene_det_in_min_cells = 5,
                      min_det_genes_per_cell = 5, expression_values = c('raw'), verbose = T)
st_obj = normalizeGiotto(gobject = st_obj, scalefactor = 6000, verbose = T)
st_obj <- addStatistics(gobject = st_obj)


st_obj <- calculateHVG(gobject = st_obj, method = 'cov_loess', 
                       difference_in_cov = 0.1,show_plot = FALSE, save_param = list(save_name = '3_a_HVGplot', base_height = 5, base_width = 5))

gene_metadata = fDataDT(st_obj)
featgenes = gene_metadata[hvg == 'yes']$gene_ID

st_obj <- runPCA(gobject = st_obj, genes_to_use = featgenes, scale_unit = F, center = F)
st_obj <- runUMAP(st_obj, dimensions_to_use = 1:5, n_threads = 20)
st_obj <- createNearestNetwork(gobject = st_obj, dimensions_to_use = 1:5, k = 10)
## Leiden clustering
st_obj <- doLeidenCluster(gobject = st_obj, resolution = 1, n_iterations = 1000)
st_obj <- newrunDWLSDeconv(gobject = st_obj, 
                           cluster_column = "leiden_clus",
                           sign_matrix = average_cell_type_expr)
write.csv(st_obj@spatial_enrichment$DWLS,'SpatialDWLS_seqFISH_10000.csv')