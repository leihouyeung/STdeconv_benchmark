
benchmark_performance <- function (test_spots_metadata_mtrx, spot_composition_mtrx) 
{
  if (!is.matrix(test_spots_metadata_mtrx)) 
    stop("ERROR: test_spots_metadata_mtrx must be a matrix object!")
  if (!is.matrix(spot_composition_mtrx)) 
    stop("ERROR: syn_spots_ls must be the list obtained from the function syn_spot_comb_topic_fun().")
  colnames(spot_composition_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
                                          ".", x = colnames(spot_composition_mtrx), perl = TRUE)
  colnames(test_spots_metadata_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
                                             ".", x = colnames(test_spots_metadata_mtrx), perl = TRUE)
  suppressMessages(require(philentropy))
  true_jsd_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), 
                          ncol = 1)
  for (i in seq_len(nrow(test_spots_metadata_mtrx))) {
    x <- rbind(test_spots_metadata_mtrx[i, ], spot_composition_mtrx[i, 
    ])
    if (sum(spot_composition_mtrx[i, ]) > 0) {
      true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, 
                                                  unit = "log2", est.prob = "empirical"))
    }
    else {
      true_jsd_mtrx[i, 1] <- 1
    }
  }
  ##calculate RMSE for each cell type
  RMSE = matrix(0L,nrow = 1,ncol = ncol(test_spots_metadata_mtrx))
  all_rmse = 0
  for (i in 1:ncol(test_spots_metadata_mtrx)){
    mse = sum((test_spots_metadata_mtrx[ ,i] - spot_composition_mtrx[ ,i])^2)
    all_rmse = all_rmse + mse
    RMSE[1,i] = sqrt(mse/nrow(test_spots_metadata_mtrx))
  }
  colnames(RMSE) = colnames(test_spots_metadata_mtrx)
  all_rmse = sqrt(all_rmse / (nrow(test_spots_metadata_mtrx) * ncol(test_spots_metadata_mtrx)))
  
  quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx, 
                                                    na.rm = TRUE), c(0.25, 0.5, 0.75)), 5)
  
  return(list(
    JSD = quants_jsd[[2]],
    RMSE = RMSE,
    Sum_RMSE = all_rmse,
    corr = cor.test(test_spots_metadata_mtrx,spot_composition_mtrx)
  ))
}