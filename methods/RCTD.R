library(spacexr)
library(Matrix)

load_seqFISH=function(n_genes){
    st_counts_fp=sprintf("../datasets/seqFISH+/Out_gene_expressions_%dgenes.csv",n_genes)
    st_locations_fp="../datasets/seqFISH+/Out_rect_locations.csv"
    sc_counts_fp="../datasets/seqFISH+/raw_somatosensory_sc_exp.txt"
    sc_labels_fp="../datasets/seqFISH+/somatosensory_sc_labels.txt"
    
    st_counts=read.csv(st_counts_fp,sep=",",row.names=1)
    st_counts=t(st_counts)
    st_locations=read.csv(st_locations_fp,sep=",",row.name=1)
    st_locations=st_locations[,c("X","Y")]
    
    sc_counts=read.csv(sc_counts_fp,sep="\t",row.names=1)
    sc_labels=read.csv(sc_labels_fp,header=FALSE)$V1
    names(sc_labels)=colnames(sc_counts)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    return(ret_list)
}

data = load_seqFISH(10000)
out_dir="../seqFISH_10000_Result"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_norm_fp=file.path(out_dir,sprintf("seqFISH_10000.RCTD.norm.csv"))

sc_reference=Reference(
    counts=data$sc_counts,
    cell_types=data$sc_labels
)

st_data=SpatialRNA(
    counts=data$st_counts,
    coords=data$st_location,
    require_int=FALSE
)

start_time <- Sys.time()

myRCTD <- create.RCTD(st_data, sc_reference, max_cores = 1, CELL_MIN_INSTANCE = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

end_time <- Sys.time()

weights=myRCTD@results$weights
norm_weights=normalize_weights(weights)
print(end_time-start_time)

write.csv(as.matrix(norm_weights),out_matrix_norm_fp)