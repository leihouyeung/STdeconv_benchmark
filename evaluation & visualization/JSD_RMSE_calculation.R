### calculation of JSD and RMSE for seqFISH+ and MERFISH

library(stringr)
source('metrics.R')


### MERFISH

data = readRDS('~/STreview/STdeconvolution_datasets/MERFISH/simMERFISH_20.RDS')
ground_truth = data$gtSpotTopics
ground_truth = as.matrix(ground_truth)
setwd('~/STreview/Results/MERFISH/20/')
ls_method = list.files()
Results = data.frame(matrix(ncol = ncol(ground_truth) +2, nrow = length(ls_method)))
for(i in 13:length(ls_method)){
  method = read.csv(ls_method[i],row.names = 1)
  colnames(method)[str_detect(colnames(method),'_')] = gsub('_','.',colnames(method)[str_detect(colnames(method),'_')])
  method = method[,colnames(ground_truth)]
  method = as.matrix(method)
  res = benchmark_performance(method,ground_truth)
  Results[i,1] = res$JSD
  Results[i,2] = res$Sum_RMSE
  Results[i,3:ncol(Results)] = res$RMSE
}
colnames(Results) = c('JSD','total_RMSE',colnames(ground_truth))
rownames(Results) = gsub('.csv','',ls_method)
write.csv(Results,'~/STreview/Results/metrics/MERFISH_20.csv')



### seqFISH+

ground_truth = read.csv('~/ST/CellGridDivide3/Out_cell_ratio_1x.csv',row.names=1)
delete_cells_index = which(is.na(ground_truth[,1]))
ground_truth = ground_truth[-delete_cells_index,]
ground_truth = as.matrix(ground_truth)
setwd('~/STreview/Results/seqFISH/10000/')
ls_method = list.files()
Results = data.frame(matrix(ncol = ncol(ground_truth) +2, nrow = length(ls_method)))
for(i in 13:length(ls_method)){
  method = read.csv(ls_method[i],row.names = 1)
  colnames(method)[str_detect(colnames(method),'_')] = gsub('_','.',colnames(method)[str_detect(colnames(method),'_')])
  method = method[,colnames(ground_truth)]
  method = as.matrix(method)
  res = benchmark_performance(method,ground_truth)
  Results[i,1] = res$JSD
  Results[i,2] = res$Sum_RMSE
  Results[i,3:ncol(Results)] = res$RMSE
}
colnames(Results) = c('JSD','total_RMSE',colnames(ground_truth))
rownames(Results) = gsub('.csv','',ls_method)

write.csv(Results,'~/STreview/Results/metrics/seqFISH_10000.csv')
