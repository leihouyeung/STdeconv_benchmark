###  seqFISH+ visualization

library(Rmisc)
library(ggplot2)
library(stringr)
source('metrics.R')

ground_truth = read.csv('~/ST/CellGridDivide3/Out_cell_ratio_1x.csv',row.names=1)
delete_cells_index = which(is.na(ground_truth[,1]))
ground_truth = ground_truth[-delete_cells_index,]
location = read.csv('~/ST/CellGridDivide3/Out_rect_locations_1x.csv',header = T,row.names = 1)
location = location[-delete_cells_index,]
for(i in 1:nrow(location)){
  location[i,c('X','Y')] = c(as.numeric(strsplit(location[i,1],split = '_')[[1]][1]),
                             as.numeric(strsplit(location[i,1],split = '_')[[1]][2]))
}

cell_type = 'microglia'
cell_type = 'eNeuron'
cell_type = 'iNeuron'
cell_type = 'endo_mural'
cell_type = 'astrocytes'
cell_type = 'Olig'


setwd('~/STreview/Results/seqFISH/10000/')
ls_method = list.files()
ls_results = list()

df = data.frame(x1 = location[,'X'], y1=location[,'Y'],
                x2 = location[,'X']+1,y2 =location[,'Y']+1,Proportion = ground_truth[,cell_type] )

title = cell_type
if(title == 'endo_mural'){
  title = 'endothelial-mural'
}else if(title == 'Olig'){
  title = 'Oligodendrocytes'
}
title = 'Ground Truth'
pt = ggplot()+geom_rect(data=df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill = Proportion),color = 'black',alpha=0.5)+
  ggtitle(title)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 22),
        axis.line = element_blank(),
        legend.position = 'none')+scale_fill_gradient(low = 'white',high = '#EE4000', limits = c(0,1))
for(i in 1:length(ls_method)){
  predict <- read.csv(ls_method[i],header=T,row.names = 1 )
  title = gsub('.csv','',ls_method[i])
  if(title == 'SD2'){
    title = expression(SD^2)
  }
  colnames(predict)[which(colnames(predict) == 'endo.mural')] = 'endo_mural'
  df = data.frame(x1 = location[,'X'], y1=location[,'Y'],
                  x2 = location[,'X']+1,y2 =location[,'Y']+1,Proportion = predict[,cell_type] )
  assign(paste('p',i,sep = ''),ggplot()+geom_rect(data=df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill = Proportion),color = 'black',alpha=0.5)+
           ggtitle(title)+
           theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 plot.title = element_text(hjust = 0.5,size = 22),
                 axis.line = element_blank(),
                 legend.position = 'none')+scale_fill_gradient(low = 'white',high = '#EE4000', limits = c(0,1)))
}
png(paste('~/STreview/figures/seqFISH_10000_',cell_type,'.png',sep = ''),width = 8000, height = 2000, res = 300)
multiplot(pt,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,cols = 6)
dev.off()
