###calculation of variance among three times seqFISH+ experiments

source('metrics.R')
library(stringr)


ground_truth = read.csv('~/ST/CellGridDivide3/Out_cell_ratio_1x.csv',row.names=1)
delete_cells_index = which(is.na(ground_truth[,1]))
ground_truth = ground_truth[-delete_cells_index,]
ground_truth = as.matrix(ground_truth)
setwd('~/STreview/Results/seqFISH/1/')
ls_method_1 = list.files()
setwd('~/STreview/Results/seqFISH/2/')
ls_method_2 = list.files()
setwd('~/STreview/Results/seqFISH/3/')
ls_method_3 = list.files()
Results = data.frame(matrix(ncol = 2, nrow = length(ls_method_1)))
for(i in 1:length(ls_method_1)){
  setwd('~/STreview/Results/seqFISH/1/')
  method1 = read.csv(ls_method_1[i],row.names = 1)
  colnames(method1)[str_detect(colnames(method1),'_')] = gsub('_','.',colnames(method1)[str_detect(colnames(method1),'_')])
  method1 = method1[,colnames(ground_truth)]
  method1= as.matrix(method1)
  res1 = benchmark_performance(method1,ground_truth)
  
  setwd('~/STreview/Results/seqFISH/2/')
  method2 = read.csv(ls_method_2[i],row.names = 1)
  colnames(method2)[str_detect(colnames(method2),'_')] = gsub('_','.',colnames(method2)[str_detect(colnames(method2),'_')])
  method2 = method2[,colnames(ground_truth)]
  method2= as.matrix(method2)
  res2 = benchmark_performance(method2,ground_truth)
  
  setwd('~/STreview/Results/seqFISH/3/')
  method3 = read.csv(ls_method_3[i],row.names = 1)
  colnames(method3)[str_detect(colnames(method3),'_')] = gsub('_','.',colnames(method3)[str_detect(colnames(method3),'_')])
  method3 = method3[,colnames(ground_truth)]
  method3= as.matrix(method3)
  res3 = benchmark_performance(method3,ground_truth)
  Results[i,1] = var(c(res1$JSD,res2$JSD,res3$JSD))
  Results[i,2] = var(c(res1$Sum_RMSE,res2$Sum_RMSE,res3$Sum_RMSE))
}
colnames(Results) = c('JSD_var','total_RMSE_var')
rownames(Results) = gsub('.csv','',ls_method_1)

write.csv(Results,'~/STreview/Results/seqFISH_var.csv')

