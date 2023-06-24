### dir_current <- "E:/Research paper/IntTemporalSpatial/AnalysisCode/ProFAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/util_funcs.R"))
setwd(paste0(dir_current, "Real_data_results/dataFiles/DLPFC12/") )

##### 2a #####
load("MacR2_sortList_DLPFC12.rds")
names(MacR2_sortList)
MacR2List <- MacR2_sortList
names(MacR2List)[4] <-  "SpatialPCA-L"
names(MacR2List)[2] <- "LIGER"
names_order <- setdiff(Methods_ordered3, c("SpatialPCA-F"))
MacR2List <- MacR2List[names_order]                 
MacMat <- Reduce(cbind, MacR2List)
colnames(MacMat) <- names(MacR2List)
new_names <- rename_vec(colnames(MacMat), oldname = "SpatialPCA-A", newname = "SpatialPCA-L")
colnames(MacMat) <- new_names
library(ggplot2)

cols_use <- col_orderd3[names_order]
p1 <- PRECAST::volinPlot(MacMat[,names_order], ylabel = "MacFadden's R2") +
  scale_fill_manual(values = cols_use[names_order]) + theme_classic(base_size = 20) +
  scale_x_discrete(breaks=NULL)
p1

load("time_sortList_DLPFC12.rds")
names(time_sortList)[4] <-  "SpatialPCA-L"
names(time_sortList)[2] <- "LIGER"
names_order <- setdiff(Methods_ordered3, c("SpatialPCA-F") )
df <- data.frame(Method= factor(names_order, levels=names_order), Time= unlist(time_sortList[names_order]))
p2<-ggplot(data=df, aes(x=Method, y=Time,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd3[names_order])+ scale_y_log10()+
  theme_classic(base_size = 20) + scale_x_discrete(breaks=NULL) + ylab("Time (sec.)")+
  xlab(NULL) + theme(legend.position = 'none')

load("memListGB_DLPFC12.rds")
memVec <- unlist(memListGB)
names(memVec) <- rename_vec(names(memVec), oldname = "SpatialPCA-A", newname = "SpatialPCA-L")
df <- data.frame(Method= factor(names_order, levels=names_order), Memory=memVec[names_order]+1)

p3<-ggplot(data=df, aes(x=Method, y=Memory,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd3[names_order])+ 
  theme_classic(base_size = 20)+scale_x_discrete(breaks=NULL) +
  scale_y_discrete(breaks=c(1,3,10,100)) +  ylab("Memory (Gigabytes)")+
  xlab(NULL) + theme(legend.position = 'none') + scale_y_log10()

drawFigs(list(p1, p2, p3), layout.dim = c(1,3), common.legend = T, legend.position = 'right', align='hv')



##### 2b #####
load('iscmeb_cluster_List_DLPFC12.rds')
load("posList1_DLPFC12.rds")
load("yList_DLPFC12.rds")

yList_layer <-  yList
names(cluster_List)

yList <- lapply(yList, factor, levels=c(unique(unlist(yList))))
layer_unorder <- c(unique(unlist(yList)))
yList <- lapply(yList, as.numeric)
cluster_List <- c(Annotation=list(yList), cluster_List)
clusterMat <- matrix(NA, nrow=length(unlist(cluster_List[[1]])), ncol=length(cluster_List))
for(i in seq_along(cluster_List)) clusterMat[,i] <- unlist(cluster_List[[i]])
colnames(clusterMat) <- names(cluster_List)

clusterMat_sub <- matrix(as.numeric(clusterMat), ncol= ncol(clusterMat))
colnames(clusterMat_sub) <- colnames(clusterMat)
apply(clusterMat_sub, 2, max)

indexList <- get_indexList(posList1)
order_names <- colnames(clusterMat_sub)
clusterMat_sub <- clusterMat_sub[, order_names]
## allign colors
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]

pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_cluster <- c(tableau_color_pal()(10), pal1(10)[c(2:4, 6,8:9)], pal2(max_n2)[-3]) # 
names(layer_unorder) <- cols_cluster[1:8]
base_cols_cluster <- 1:8
base_cluster <- clusterMat_sub[,1]
ClusterMat <- clusterMat_sub[,2:ncol(clusterMat_sub)]
align_colors(ClusterMat[,1:2], base_cluster, cols_cluster)
colorList <- c(list(cols_cluster[1:8]), align_colors(ClusterMat, base_cluster, cols_cluster))
colorList <- lapply(colorList, function(x){
  names(x) <- 1:length(x)
  return(x)
})
names(colorList) <- c("Annotation",colnames(ClusterMat))
M <- length(yList)
nMethod <-  ncol(clusterMat_sub)
## Plot sample 10
pList10 <- list()
for(r in 10){ ## each sample
  # r <- 1
  pList_eachmethod <- list()
  for(j in 1: nMethod){ ## each method
    # j <- 1
    
    message("r = ", r)
    cluster_tmp <- factor(clusterMat_sub[indexList[[r]],j])
    p_tmp <- plot_scatter(posList1[[r]], meta_data = data.frame(cluster=cluster_tmp),
                          label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                          point_size = 1.2, cols = colorList[[j]][as.numeric(levels(cluster_tmp))], point_alpha = 0.7) + 
      PRECAST:::mytheme_graybox()
    
    
    pList10[[j]]<- rotate_90_clockwise(p_tmp)
  }
}

names(pList10) <- colnames(clusterMat_sub)
method_chose <- c("Annotation", "ProFAST-P", "SpatialPCA-A", "DR-SC", "scVI","multiBatchPCA")
p12 <- cowplot::plot_grid(plotlist =pList10[method_chose], nrow=1, ncol=6)
ggsave(file="indi10_iSCME_cluster_SpaHeatmap.png", plot = p12,
       width = 15, height =3, units = "in", dpi = 500)

### RGB plot
load('UMAP3List_all_correct_DLPFC12.rds')
load("posList1_DLPFC12.rds")
names(UMAP3List_all_correct)
names(UMAP3List_all_correct)[c(7:8)] <- c("ProFAST-G", "ProFAST-P")
## take out sample 10
UMAP3List_all_correct <- UMAP3List_all_correct[names(col_orderd2)[-3]]
umap10List <- lapply(UMAP3List_all_correct, function(x) x[idx_List[[10]],])
### Plot 151674, r=10
pList_umap_RGB_r10 <- list()
for(i in 1:10){ ## each method
  
  message("i = ", i)
  r <- 10
  ptmp <- plot_RGB(posList1[[r]], umap10List[[i]], pointsize = 1) + mytheme_graybox()
  ptmp <- rotate_90_clockwise(ptmp) 
  pList_umap_RGB_r10[[i]] <- ptmp
}
names(pList_umap_RGB_r10) <- names(umap10List)
subnames <- c("ProFAST-P",  "SpatialPCA-A",  "DR-SC",   "multiBatchPCA", "NMF")
library(cowplot)
p12 <- plot_grid(plotlist =pList_umap_RGB_r10[subnames], nrow=1, ncol=5)
ggsave(file="umap_correct_RGB_r10_DLPFC12.png", plot = p12,
       width = 15, height =3, units = "in", dpi = 200)
##### 2c #####
load("clusterPF_List_iSCMEB_all_DR9_DLPFC12.rds")
clusterPF_List <- iscmeb_clusterPF_allList
names(clusterPF_List)[c(4,7:8)] <- c("SpatialPCA-L", "ProFAST-G", "ProFAST-P")
ariList <- lapply(clusterPF_List, function(x) x[,1])
nmiList <- lapply(clusterPF_List, function(x) x[,2])
ariMat <- Reduce(cbind, ariList); colnames(ariMat) <- names(ariList)
nmiMat <- Reduce(cbind, nmiList); colnames(nmiMat) <- names(ariList)
apply(ariMat, 2, median)
index_method_remove <- c(3,11)
name_order2 <- names(col_orderd3[-index_method_remove])
setdiff(colnames(ariMat), name_order2)
p1 <- PRECAST::volinPlot(ariMat[-13, name_order2], ylabel = "ARI", cols=col_orderd3[name_order2])
p2 <- PRECAST::volinPlot(nmiMat[-13, name_order2], ylabel = "NMI", cols=col_orderd3[name_order2])
drawFigs(list(p1,p2), layout.dim = c(1,2), common.legend = T, legend.position = 'bottom')


##### 2d #####
load('seu_dot_DLPFC12.rds')
col_here <- c("#F2E6AB", "#9C0141") 
library(ggplot2)
p1 <- DotPlot(seu_dot, features=row.names(seu_dot), cols=col_here, #  idents = ident_here,
              col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic"))+
  ylab("Domain") + xlab(NULL) + theme(axis.text.x = element_text(size=12, angle = 25, hjust = 1, family='serif'),
                                      axis.text.y = element_text(size=12, face= "italic", family='serif'))

write_fig(p1, filename = "DotPlot_profastP_combine_DLPFC12",
          width = 12, height = 8.1, dpi=500)

##### 2e #####
## see data analysis folder

##### 2f #####
## Fig. 2f is plotted in Python; see sep_paga_dlpfc12.ipynb script

