### dir_current <- "E:/Research paper/IntTemporalSpatial/AnalysisCode/ProFAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/util_funcs.R"))
setwd(paste0(dir_current, "Real_data_results/dataFiles/BC2/") )

##### 3b ###########
load("posList_xeHBC2.rds")
load("umap3_pois_xeHBC2.rds")
indexList <- get_indexList(posList)
pList_umap_RGB <- list()
for(r in 1:2){## each sample
  message("r = ", r)
  ptmp <- plot_RGB(posList[[r]], umap3_pois[indexList[[r]],], pointsize = 0.8) +
    PRECAST:::mytheme_graybox()+  geom_point(alpha=0.5)
  pList_umap_RGB[[r]]<- ptmp 
}


p12 <- cowplot::plot_grid(plotlist =pList_umap_RGB, nrow=1, ncol=2)

ggsave(file="BC2_umap_RGB_iSCMEB_pois.png", plot = p12,
       width = 12.1, height =6, units = "in", dpi = 80)
##### 3c ###########

load("posList_xeHBC2.rds")
color_use <- chooseColors(palettes_name = "Blink 23", n_colors=19)[-16]
simutool::colorbar_adj_transparent(color_use, alpha=1)
names(color_use) <- 1:length(color_use)
load("hyList_iscmeb_profastP_renumberV2_xeHBC2.rds")

pList <- list()
for(r in 1:length(posList)){ # hyList_iscmeb_pois_renumber
  #r <- 1
  p_tmp <- plot_scatter(posList[[r]], meta_data = data.frame(cluster=factor(hyList_iscmeb_pois_renumberV2[[r]], levels=1:18)),
                        label_name = 'cluster', xy_names = c("", ""), no_guides = F,
                        point_size = 1.0, cols =color_use, point_alpha = 0.9) + 
    PRECAST:::mytheme_graybox()
  pList[[r]] <- p_tmp
  
}
p_all <- drawFigs(pList, layout.dim = c(1,2), common.legend = T, legend.position = 'none')

ggsave(file="iscmeb_pois_reorder_xeHBC2.png", plot = p_all,
       width = 12.1, height =6, units = "in", dpi = 100)
##### 3d ###########
load("timeList_method6_xeHBC2.rds")
names(timeList)
index_method_remove <- c(3)
names_order <- names(col_orderd2[-index_method_remove])
timeList <- timeList[names_order]
timeVec <- Reduce(c, timeList)
names(timeVec) <- names(timeList[names_order])

name_use <- c("ProFAST-P" ,    "ProFAST-G" ,    "SpatialPCA-A",  "PRECAST" )
name_choose <-name_use

level_time <- names(col_orderd2)
df <- data.frame(Method= factor(name_choose, levels=level_time), Time= timeVec[name_choose])

p2 <-ggplot(data=df, aes(x=Method, y=Time,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd2[level_time])+
  theme_classic() + scale_x_discrete(breaks=NULL) + theme(text=element_text(size=20))+ ylab("Time (sec.)") +
  xlab(NULL) + scale_y_log10()+ theme(legend.position = 'none')

load("memListGB_xeHBC3.rds")
memVec <- unlist(memListGB)
new_names <- rename_vec(names(memVec), oldname = "SpatialPCA-A", newname = "SpatialPCA-L")
names(memVec) <- new_names
index_method_remove <- c(3)
name_use <- names(col_orderd3[-index_method_remove])
# name_use <- c("ProFAST-P" ,    "ProFAST-G" ,    "SpatialPCA-A",  "DR-SC" )

df_raw <- data.frame(Method= factor(name_use, levels=name_use), Memory= memVec[name_use]) #

### choosed methods:
name_use <- c("ProFAST-P", "ProFAST-G" ,    "SpatialPCA-L",   "PRECAST" )
df <- data.frame(Method= factor(name_use, levels=name_use), Memory= memVec[name_use]) #+1

p_mem_cho <- ggplot(data=df, aes(x=Method, y=Memory,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd3[name_use])+
  theme_classic() + scale_x_discrete(breaks=NULL) + theme(text=element_text(size=20))+ ylab("Memory (Gigabytes)")+
  xlab(NULL) + theme(legend.position = 'none')
# p_mem_cho <- p_mem_cho+scale_y_log10() # + ylim(c(0, 12)) # 
drawFigs(list(p2, p_mem_cho), layout.dim = c(1,2), common.legend = T, legend.position = 'right', align='hv')


##### 3e ###########
load("hyList_iscmeb_profastP_xeHBC2.rds")
prop_two_slices <- sapply(hyList_iscmeb_pois, function(x) table(x)/length(x))
dat <- as.data.frame(prop_two_slices)
colnames(dat) <- paste0("Slice", 1:2)
model2 <- lm(Slice2~Slice1, data=dat)
s2 <- summary(model2)
R2_2 <- s2$r.squared 
R2_2
p2 <- ggplot(dat, aes(x=Slice1, y=Slice2)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, color="blue") +
  labs(title=NULL, 
       x="Proportion 1", y="Proportion 2") +
  theme(plot.title = element_text(hjust = 0.5))+ theme_classic() + theme(text=element_text(size=18))

write_fig(p2, filename = "prop_twoSlices_linearReg_profastP", width = 5, height = 5)

##### 3f ###########
load('seu_bc2.rds')
col_here <- c("#F2E6AB", "#9C0141")
p3 <- DotPlot(seu_bc2, features=row.names(seu_bc2), cols=col_here, #  idents = ident_here,
              col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic", size=11))+
  ylab("Cluster") + xlab(NULL) #+ theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
p3

##### 3g ###########
load('dat_kegg_all_HBC2.rds')
ggplot(data = dat_kegg_all, mapping = aes_string(x = "term_name", y = "Cluster")) + 
  geom_point(mapping = aes_string(size = "nlog10P"))+
  theme_bw(base_size=15) + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,family='serif'), 
                                 plot.margin=unit(c(1,1,1,4),'lines'),
                                 axis.text.y = element_text(family='serif'))

##### 3h ###########
load("viper_poisson_ARACNe_res.rda")
### Load new order of clusters
load('hyList_iscmeb_profastP_renumberV2_xeHBC2.rds')
for(r in 1:2){
  Idents(viper_res[[r]]) <- factor(hyList_iscmeb_pois_renumberV2[[r]], levels=1:17)
}


pList1 <- list()
pList2 <- list()
seu <- viper_res[[1]]
features = rownames(seu[["dorothea"]])
for(i in 1:length(features)){
  message("i = ", i)
  pList1[[i]] <- VlnPlot(viper_res[[1]], group.by = "cluster", assay = "dorothea", features = features[i],
                         pt.size = 0, cols=color_use, raster =T) + NoLegend() + ylab("Activity score") + xlab("Cluster")+
    theme(axis.text.x = element_text(hjust = NULL,vjust = NULL, angle = 0), text=element_text(size=16))
  
  pList2[[i]] <- VlnPlot(viper_res[[2]], group.by = "cluster", assay = "dorothea", features = features[i],
                         pt.size = 0, cols=color_use, raster =T) + NoLegend() + ylab("Activity score") + xlab("Cluster")+
    theme(axis.text.x = element_text(hjust = NULL,vjust = NULL, angle = 0), text=element_text(size=16))
  
}
pList1[[1]]
library(cowplot)
p11 <- plot_grid(plotlist = pList1, ncol=1)
ggsave(file="activityScore_xenium1.png", plot = p11,
       width = 5, height =20, units = "in", dpi = 200)

p12 <- plot_grid(plotlist = pList2, ncol=1)
ggsave(file="activityScore_xenium2.png", plot = p11,
       width = 5, height =20, units = "in", dpi = 200)

##### 3i ###########

scaleFUN <- function(x) sprintf("%.2f", x)

bar_plot <- function(percentage_long,geom_bar_position,legend_position,color_pal,
                     base_text_size=20){
  ggplot(percentage_long, aes(y = value, x = factor(Cluster), fill = variable)) +        ## global aes
    scale_fill_manual(values= color_pal,name = 'Cell Type')+
    scale_y_continuous(labels=scaleFUN) +
    geom_bar(position=geom_bar_position, stat="identity",width=0.7,color="black") +
    ggtitle(paste("")) +
    theme_bw()+xlab("")+ylab("")+
    theme_classic() +
    theme(plot.title = element_text(size = base_text_size,hjust = 0.5),
          text = element_text(size = base_text_size-2),
          axis.text = element_text(size = base_text_size-4),
          axis.line = element_line(colour = "grey"),
          legend.position = legend_position,
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "gray10"))
}

library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]

pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_cluster <- c(pal1(max_n), pal2(max_n2)[c(1,3,8,9,12)])

load('norm_weightsList_merged_celltypes.rds')
norm_weights <- NULL
for(s in 1:2){
  # s<- 1
  tmp_weights <- norm_weightsList[[s]]
  row.names(tmp_weights) <- paste0(row.names(tmp_weights), "_xeHBC", s)
  norm_weights <- rbind(norm_weights, tmp_weights)
}

library(dplyr)
library(reshape2)
library(ggplot2)
figure_list <- list()
need_scale = T
pos_namesList <- lapply(1:2, function(i) paste0(row.names(posList[[i]]), "_xeHBC", i))   

load("hyList_iscmeb_profastP_renumberV2_xeHBC2.rds")

cluster_id = cbind(unlist(pos_namesList),unlist(hyList_iscmeb_pois_renumberV2))
row.names(cluster_id) <- cluster_id[,1]
cluster_weight = merge(norm_weights,cluster_id[,2],by = 0)
rownames(cluster_weight) = cluster_weight$Row.names
cluster_weight = cluster_weight[,-1]

percentage = as.data.frame(cluster_weight %>% group_by(y) %>% summarise(across(everything(), sum)))
#colnames(percentage) = c( "Cluster",'Immune cell',"CAF","TAM","HPC-like cell","Malignant cell")
colnames(percentage)[1] = c( "Cluster")
percentage_long <- melt(as.data.frame(percentage),id.vars ='Cluster')


scale_flag <- c(TRUE,FALSE)

for(ii in 1:2){
  if (need_scale == scale_flag[ii]){
    percentage_scale = matrix(0,nrow(percentage),ncol(percentage))
    for (tt in 1:nrow(percentage)){
      percentage_scale[tt,] = as.numeric(cbind(percentage$Cluster[tt],percentage[tt,2:ncol(percentage)]/sum(rowSums(percentage[,2:ncol(percentage)]))))
    }
    colnames(percentage_scale) = colnames(percentage)
    #sum(percentage_scale[,2:ncol(percentage)])
    percentage_long <- melt(as.data.frame(percentage_scale),id.vars ='Cluster')
  }
  
  if (need_scale == scale_flag[ii]){
    geom_bar_position = 'stack'
  }else{
    geom_bar_position = 'fill'
  }
  
  level_names <- sort(levels(percentage_long$variable))
  percentage_long$variable <- factor(percentage_long$variable, levels=level_names)
  
  percentage_long$Cluster <- factor(percentage_long$Cluster, levels = 1: length(unique(percentage_long$Cluster)))
  
  if (need_scale != scale_flag[ii]){
    figure_list[[ii]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, 
                                  legend_position='right',
                                  color_pal =cols_cluster, base_text_size = 20) + theme(legend.position = 'none')
  }else{
    figure_list[[ii]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, 
                                  legend_position='none',
                                  color_pal = cols_cluster,
                                  base_text_size = 20)
  }
  
  
}
figure_list[[1]] + theme(legend.position = 'right') #bottom
figure_list[[2]] + theme(legend.position = 'right') #bottom
library(ggpubr)
ggsave(file='./result/combined_percentage_scaled_pois_renumber.png', plot =figure_list[[1]], 
       width = 7.5, height = 5.1, units = "in", bg = 'white', dpi = 500,limitsize = F)

ggsave(file='./result/combined_percentage_unscaled_pois_renumber.png', plot =figure_list[[2]], 
       width = 7.5, height = 5.1, units = "in", bg = 'white', dpi = 500,limitsize = F)

##### 3j ###########
load("posList_xeHBC2.rds")
load("norm_weightsList_merged_celltypes.rds")
load("hyList_iscmeb_profastP_renumberV2_xeHBC2.rds")
norm_weights <- norm_weightsList[[1]]
# Cancer regions
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
cell_type <- "Cancer"

message("cell_type=", cell_type)
  figure_list_malig = list()
  cols_use <- c("#3AB370","white", "#D62728")#c("#D62728","white", "#3AB370")
  for (slice_set in 1:2){
    # slice_set <- 1
    st_coord = as.data.frame(posList[[slice_set]])
    colnames(st_coord) <- c("imagerow", "imagecol")
    
    norm_weights = norm_weightsList[[slice_set]]
    
    cluster_id = cbind(as.matrix(st_coord),hyList_iscmeb_pois_renumberV2[[slice_set]])
    cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
    rownames(cluster_weight) = cluster_weight$Row.names
    cluster_weight = cluster_weight[,-1]
    
    # st_coord$imagerow = max(st_coord$imagerow)-st_coord$imagerow
    
    
    plot_val = as.data.frame(cluster_weight[,cell_type])
    rownames(plot_val) = rownames(cluster_weight)
    colnames(plot_val) = cell_type
    barcodes = rownames(plot_val)
    my_table = st_coord[barcodes, ]
    my_table$Proportion = plot_val[barcodes,]
    my_table$Proportion <- range01(my_table$Proportion)
    med <- quantile(my_table$Proportion, 0.94)
    my_table$Proportion[ my_table$Proportion > med] <- med
    
    
    my_table$cluster = cluster_weight[barcodes,'y']
    ylimit = c(0, 1)
    
    my_table$border_color = '#cccccc'
    my_table$stroke = 0
    my_table <- my_table[order(my_table$Proportion, decreasing = F),]
    
    table(my_table$cluster)
    head(my_table)
    plot <- ggplot(my_table, aes(x = imagerow, y = imagecol)) +  # stroke = my_table$stroke,
      geom_point(aes(fill = Proportion),size = 1.3, alpha = 3,pch=21) +
      #scale_fill_gradientn(colors = my_pal,limits = ylimit) +  # ,colour = my_table$border_color
      scale_fill_gradientn(colors = cols_use) +
      scale_shape_identity() + scale_size_identity() + theme_classic() +
      mytheme_graybox(base_size = 28,  border_color = 'black')
    
    # xlim <- c(min(st_coord$imagecol) - 1, max(st_coord$imagecol) + 1)
    # ylim <- c(min(st_coord$imagerow) - 1, max(st_coord$imagerow) + 1)
    figure_list_malig[[slice_set]] <- plot 
    if (slice_set == 1){
      figure_list_malig[[slice_set]] <- figure_list_malig[[slice_set]] + ggtitle(cell_type)
    }
  }
  
  figure_celltype_weight <- ggpubr::ggarrange(plotlist = figure_list_malig,ncol=2,nrow=1,
                                              align='hv',common.legend = T,
                                              legend = 'right')
write_fig(figure_celltype_weight, filename = paste0(cell_type,"_select_type_gauss.png"),
            width = 10, height = 6, dpi = 200)


##### 3k ###########
load('moransI_xeHBC2.rds')
index_method_remove <- c(3)
p1 <- PRECAST::boxPlot(mr1[,names(col_orderd3[-index_method_remove])], ylabel="Moran's I", 
                       cols=col_orderd3[-index_method_remove]) +
  theme_classic(base_size = 20)
p2 <- PRECAST::boxPlot(mr2[,names(col_orderd3[-index_method_remove])], ylabel="Moran's I", 
                       cols=col_orderd3[-index_method_remove]) +
  theme_classic(base_size = 20)
drawFigs(list(p1, p2), c(1, 2), common.legend = T)

