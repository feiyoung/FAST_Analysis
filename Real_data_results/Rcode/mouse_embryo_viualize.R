rm(list=ls())
### dir_current <- "E:/Research paper/IntTemporalSpatial/AnalysisCode/ProFAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/util_funcs.R"))

setwd(paste0(dir_current, "Real_data_results/dataFiles/ME26/") )

##### 5a ###########
M <- 26
load("MacR2_Mat_method5_E12E16.rds")
#boxplot(MacR2_Mat[-(M+1), ])
name_id <- which(names(col_orderd2) %in% colnames(MacR2_Mat))
names_order <- names(col_orderd2)[name_id][-1] # remove poisson
library(ggplot2)
cols_use <- col_orderd2
p1 <- PRECAST::volinPlot(MacR2_Mat[-(M+1),names_order], ylabel = "MacFadden's R2") +
  scale_fill_manual(values = cols_use[names_order]) + theme(legend.position = 'right') + ylab(NULL)
p1
write_fig(p1, filename = 'MacR2_comp_method4_E12E16', height = 5, width = 4, dpi=500)

##### 5b ###########

M <- 26
load("clusterPF_List_E12E16.rds")
name_id <- which(names(col_orderd2) %in% names(clusPF_iscmeb_List))
names_order <- names(col_orderd2)[name_id][-1] # remove ProFAST-P
ariMat <- sapply(clusPF_iscmeb_List, function(x) x[,1])
nmiMat <- sapply(clusPF_iscmeb_List, function(x) x[,2])
p_ari <- PRECAST::volinPlot(ariMat[-(M+1),names_order], ylabel = "ARI") +
  scale_fill_manual(values = col_orderd2[names_order]) + theme(legend.position = 'none') + ylab(NULL)
p_nmi <- PRECAST::volinPlot(nmiMat[-(M+1),names_order], ylabel = "NMI") +
  scale_fill_manual(values = col_orderd2[names_order]) + theme(legend.position = 'none') + ylab(NULL)
drawFigs(list(p_ari, p_nmi), layout.dim = c(1,2), common.legend = T, align='hv', legend.position = "right")


##### 5c ###########
load('posList1_E12E16.rds')
load('UMAP3_correct_Listall_E12E16.rds')

E_list <- list(5:10, 11:14, 15:21, 22: 26, 27:30)
E_list <- lapply(E_list, function(x) x- 4)
slice_ids <- sapply(E_list, function(x) x[1])
indexList <- get_indexList(posList1)
Methods_use <- c(1,2,4)
names(UMAP3_correct_Listall)
pList_select <- list()
kk <- 1
for(r in slice_ids){ ## each sample
  pList_eachmethod <- list()
  for(j in 1:3){ 
    # j <- 1
    j_method <- Methods_use[j]
    message("r = ", r)
    umap3_tmp <- UMAP3_correct_Listall[[j_method]]
    ptmp <- plot_RGB(posList1[[r]], umap3_tmp[indexList[[r]],], pointsize = 0.8) + scale_y_reverse()+ mytheme_graybox(border_color = "white") # 
    pList_select[[kk]]<- ptmp
    kk <- kk + 1
  }
}


library(cowplot)
p12 <- plot_grid(plotlist =pList_select, nrow=3, ncol=5, byrow = F) # , rel_widths = seq(0.6,1, by=0.1), rel_heights = seq(0.6,1, by=0.1)
ggsave(file=paste0("UMAP_RGB_method3_SpaHeatmap.png"), plot = p12,
       width = 10, height =7, units = "in", dpi = 50)


### clustering assignment
load("yList_E12E16.rds")
yList_str <- lapply(yList, factor, levels=c(unique(unlist(yList))))
tissue_names <- as.character(unique(unlist(yList_str)))
names(tissue_names) <- 1:38
yList <- lapply(yList_str, as.numeric)
load("clusterList_profastG_E12_E16.rds")
load("iscmeb_clusterList_allMethods_E12E16.rds")
names(iscmeb_cluster_allList_methodAll)
cluster_List <- c(list("Annotation"=yList, "ProFAST-G"=clusterList_profastG), 
                  iscmeb_cluster_allList_methodAll[c("PCA",  "multiBatchPCA", "NMF")])

clusterMat <- matrix(NA, nrow=length(unlist(cluster_List[[1]])), ncol=length(cluster_List))
for(i in seq_along(cluster_List)) clusterMat[,i] <- unlist(cluster_List[[i]])
colnames(clusterMat) <- names(cluster_List)

clusterMat_sub <- matrix(as.numeric(clusterMat), ncol= ncol(clusterMat))
colnames(clusterMat_sub) <- colnames(clusterMat)

load("colorList_method5_E12E16.rds")
M <- length(yList)
nMethod <-  ncol(clusterMat_sub)

E_list <- list(5:10, 11:14, 15:21, 22: 26, 27:30)
E_list <- lapply(E_list, function(x) x- 4)
slice_ids <- sapply(E_list, function(x) x[1])
colnames(clusterMat_sub)
Methods_use <- c(1:3, 5)
pList_select <- list()
kk <- 1
rflag <- 0
for(r in slice_ids){ ## each sample
  # r <- 3
  rflag <- rflag + 1
  pList_eachmethod <- list()
  for(j in 1: 4){ ## each method
    # j <- 1
    j_method <- Methods_use[j]
    message("r = ", r)
    cluster_tmp <- factor(clusterMat_sub[indexList[[r]],j_method])
    #xrange <- range(posList1[[r]][,1]) + slice_ids[5-(rflag)+1]*c(-10, 10)
    #yrange <- range(posList1[[r]][,2]) + slice_ids[5-(rflag)+1]*c(-80, 80)
    p_tmp <- plot_scatter(posList1[[r]], meta_data = data.frame(cluster=cluster_tmp),
                          label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                          point_size = 0.5, cols = colorList[[j_method]][as.numeric(levels(cluster_tmp))], point_alpha = 0.7) + 
      PRECAST:::mytheme_graybox(border_color = "white") +  scale_y_reverse() 
    # xlim(xrange) + ylim(yrange) +
    
    pList_select[[kk]]<- p_tmp
    kk <- kk + 1
  }
}

library(cowplot)
p12 <- plot_grid(plotlist =pList_select, nrow=4, ncol=5, byrow = F) # , rel_widths = seq(0.6,1, by=0.1), rel_heights = seq(0.6,1, by=0.1)
ggsave(file=paste0("Allslices1_iSCMEB_cluster_SpaHeatmap.png"), plot = p12,
       width = 10, height =10, units = "in", dpi = 50)



##### 5d ###########
## The data file is too large to be illustrated; so the code can be found in Real_data_analysis folder

##### 5e ###########
library(Seurat)
load('seu_brain_subfeature_E12E16.rds')
pList <- list()
for(j in 1: 10){
  #j <- 1
  message("j = ", j)
  p1 <- VlnPlot(seu_brain_subfeature, features=row.names(seu_brain_subfeature)[j], raster=F, pt.size=0.5)+xlab("Embryo day") + 
    theme(text=element_text(size=16), legend.position = 'none') 
  pList[[j]] <- p1
  
}
idx_main <- c(1, 4, 5, 6)
plist <- pList[idx_main]
library(cowplot)
p12 <- plot_grid(plotlist = plist, nrow=2, ncol=2, align=c("hv"))

ggsave(filename = paste0("brain_marker_main4.png"), plot=p12,
       width = 9.5, height = 8.4, dpi=50)
