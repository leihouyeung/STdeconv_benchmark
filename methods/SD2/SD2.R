
source(SD2_utiles.R)
current_path = getwd()
setwd('~/STreview/STdeconvolution_datasets/seqFISH+/')
org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)

sc_exp = read.table('raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F)
st_location = read.csv('Out_rect_locations.csv',header = T, row.names = 1)
colnames(st_location) = c('x','y')

cell_type = sc_anno[,1]

### apply the SD2 to PDAC
setwd(current_path)
SD2(as.matrix(sc_exp),
    t(as.matrix(org_st_count)),
    cell_type,
    ST_location = st_location[,c('x','y')],
    spot_num = 1000, 
    lower_cellnum = 10,
    upper_cellnum = 20)

system('python train.py')
