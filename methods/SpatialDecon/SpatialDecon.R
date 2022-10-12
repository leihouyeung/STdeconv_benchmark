library(SpatialDecon)
library(SeuratObject)
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

preprocess=function(data){
    st_counts_norm = sweep(data$st_counts, 2, colSums(data$st_counts), "/") * mean(colSums(data$st_counts))
    st_object=CreateSeuratObject(counts=st_counts_norm,assay="Spatial")
    stopifnot(setequal(colnames(st_object),rownames(data$st_location)))
    st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),1],col.name="x")
    st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),2],col.name="y")
    
    stopifnot(all(colnames(data$sc_counts)==names(data$sc_labels)))
    
    sc_counts_matrix=as.matrix(data$sc_counts)
    sc_counts_matrix=Matrix::Matrix((sc_counts_matrix),sparse=TRUE)
    sc_labels_df=data.frame(cell_barcodes=names(data$sc_labels),sc_labels=data$sc_labels)
    sc_matrix <- create_profile_matrix(
        mtx = sc_counts_matrix,            # cell x gene count matrix
        cellAnnots = sc_labels_df,  # cell annotations with cell type and cell name as columns 
        cellTypeCol = "sc_labels",  # column containing cell type
        cellNameCol = "cell_barcodes",           # column containing cell ID/name
        matrixName = "custom_cell_type_matrix", # name of final profile matrix
        outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
        normalize = TRUE,                # Should data be normalized? 
        minCellNum = 1,
        minGenes = 1
    ) 
    
    return(
        list(
            st_object=st_object,
            sc_matrix=sc_matrix
        )
    )
}

data = load_seqFISH(10000)
out_dir="../seqFISH_10000_Result"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_norm_fp=file.path(out_dir,sprintf("seqFISH_10000.SpatialDecon.norm.csv"))

processed_data=preprocess(data)
start_time <- Sys.time()
res = runspatialdecon(object = processed_data$st_object,
                      bg = 0.01,
                      X = processed_data$sc_matrix,
                      align_genes = TRUE)
end_time <- Sys.time()

weights=t(res$beta)
norm_weights=sweep(weights, 1, rowSums(weights), "/")
write.csv(as.matrix(norm_weights),out_matrix_norm_fp)