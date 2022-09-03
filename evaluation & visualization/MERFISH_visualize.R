library(reshape2)
library(easyGgplot2)
library(gridExtra)

source('metrics.R')


data = readRDS('~/STreview/STdeconvolution_datasets/MERFISH/simMERFISH_20.RDS')
colors = c("#C0392B",'#8E44AD',"#2ECC71","#34495E","#F1C40F","#3498DB")
slice_index = data$annotDf[,c('Bregma','patch_id')]
slice_ids = sort(unique(slice_index[,'Bregma']) )
setwd('~/STreview/Results/MERFISH/100/')
ls_method = list.files()
plot_list = list()
for(i in 1:length(slice_ids)){
  selected_spots = slice_index[which(slice_index[,'Bregma'] == slice_ids[i]),'patch_id']
  gt = data$gtSpotTopics
  gt = gt[unique(selected_spots),order(colnames(gt))]
  ct.visualize = colnames(gt)
  location = as.data.frame(data$st_location[unique(selected_spots),])
  location = location/100
  bind_data = cbind(gt,location)
  ct.select = colnames(gt)
  title = 'Ground Truth'
  assign(paste('p',i,'gt',sep = '_'),ggplot() + 
           geom_scatterpie(aes(x=X, y=Y,r = 0.13),data=bind_data,
                           cols=ct.select,color=NA) + coord_fixed(ratio = 1) + 
           scale_fill_manual(values =  colors)+
           theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                 panel.background = element_blank(),
                 plot.background = element_blank(),
                 panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                 axis.text =element_blank(),
                 axis.ticks =element_blank(),
                 axis.title =element_blank(),
                 plot.title = element_text(hjust = 0.5,size = 50),
                 legend.position = 'none'))+guides(fill='none')
  plot_list = append(plot_list,get(paste('p',i,'gt',sep = '_')))
  
  for(j in 1:length(ls_method)){
    predict <- read.csv(ls_method[j],header=T,row.names = 1 )
    title = gsub('.csv','',ls_method[j])
    if(title == 'SD2'){
      title = expression(SD^2)
    }
    print(title)
    #colnames(predict)[which(colnames(predict) == 'endo.mural')] = 'endo_mural'
    predict = predict[unique(selected_spots),order(colnames(predict))]
    ct.visualize = colnames(predict)
    location = as.data.frame(data$st_location[unique(selected_spots),])
    location = location/100
    
    bind_data = cbind(predict,location)
    ct.select = colnames(predict)
    assign(paste('p',i,j,sep = '_'),ggplot() +geom_scatterpie(aes(x=X, y=Y,r = 0.13),data=bind_data,
                                                              cols=ct.select,color=NA) + coord_fixed(ratio = 1) + 
             scale_fill_manual(values =  colors)+
             theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                   panel.background = element_blank(),
                   plot.background = element_blank(),
                   panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                   axis.text =element_blank(),
                   axis.ticks =element_blank(),
                   axis.title =element_blank(),
                   plot.title = element_text(hjust = 0.5,size = 50),
                   legend.position = 'none'))+guides(fill='none')
    plot_list = append(plot_list,get(paste('p',i,j,sep = '_')))
    
    
  }
}

tiff('~/STreview/figures/MERFISH/MERFISH_20.tiff',width = 8000, height = 2000, res = 600)
grid.arrange(p_1_gt,p_1_1,p_1_2,p_1_3,p_1_4,p_1_5,p_1_6,p_1_7,p_1_8,p_1_9,p_1_10,p_1_11,p_1_12,p_1_13,p_1_14,p_1_15,p_1_16, nrow = 1, ncol = 17)
dev.off()

