### dir_current <- "E:/Research paper/IntTemporalSpatial/AnalysisCode/ProFAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/util_funcs.R"))

setwd(paste0(dir_current, "Real_data_results/dataFiles/HCC4/") )

##### 4a ###########

load("HCC4_posList1.rds")
posList <- posList1
load("clusterList_iscmeb_pois_reordered_K9_HCC4.rds")

pt_size_sample <- 0.2; pt_alpha_sample <- 0.5
pt_size_cluster <- 0.2; pt_alpha_cluster <- 0.5
base_axis_size <- 20

red_colors8 <- c("#8A0B0B",'#D30F0D', "#D82929", "#F43835",  "#E35353","#F98787", "#F9AFB0", "#FDD6D7")
blue_colors3 <- c("#282DEE" ,  "#0259B6", "#545BC3")
cols_clusters_pois <- c(red_colors8[-c(2,5)], blue_colors3)
#cols_clusters_pois <- cols_clusters
names(cols_clusters_pois) <- 1:9

pList2 <- list()
M <- 4
for(r in 1:M){
  
  i_cluster <- sort(as.numeric(unique(cluster_pois_renumber_K9[[r]])))
  meta_data <- data.frame(cluster = factor(cluster_pois_renumber_K9[[r]], levels=i_cluster))
  pList2[[r]] <- plot_scatter(posList[[r]], meta_data, label_name = 'cluster', 
                              base_size=base_axis_size,cols=cols_clusters_pois[i_cluster], point_size = 1.5, 
                              point_alpha = pt_alpha_sample,no_guides=F) + labs(x=NULL, y=NULL)+
    PRECAST:::mytheme_graybox(border_color = 'white') + xlim(c(0,77)) + ylim(c(1,120))
  
}
drawFigs(pList2, layout.dim = c(1,4), common.legend = F, legend.position = 'none')

##### 4b ###########
load("ccorMat_hV_precast_profast_hcc4.rds")
colnames(ccorMat) <- paste0("HCC", 1:4)
row.names(ccorMat) <- paste0(1:15)
median_ccor <- apply(ccorMat, 2, median) # max
mean(median_ccor)
library(ggplot2)
df_ccor <- data.frame(Sample = names(median_ccor), CCor= median_ccor)

p_ccor <- ggplot(data=df_ccor, aes(x=Sample, y=CCor,  fill=Sample)) +
  geom_bar(stat="identity", color='black')+ #scale_fill_manual(values=col_orderd2)+
  theme_classic() + scale_x_discrete(breaks=NULL) + theme(text=element_text(size=20))+ ylab("Median CCor")+
  xlab(NULL) +  theme(legend.position = 'bottom') +guides(fill= guide_legend(nrow = 2))
p_ccor

##### 4c ###########
load("tSNE2_iscmeb_pois_K9_HCC4.rds")
load("sampleID_HCC4.rds")
load("clusterList_iscmeb_pois_reordered_K9_HCC4.rds")

cols_sample <- chooseColors(palettes_name = "Hue n", n_colors = 4, plot_colors = T)
names(cols_sample) <- 1:4
cols_clusters <-  cols_clusters_pois
base_axis_size <- 12; pt_size_sample <- 1
# 
meta_data <- data.frame(sample=factor(sampleID), Domain = factor(unlist(cluster_pois_renumber_K9), levels=1:9))
p1 <- plot_scatter(tSNE2_iscmeb_pois, meta_data, label_name = 'sample', 
                   base_size=base_axis_size,cols=cols_sample, point_size = pt_size_sample, 
                   point_alpha = 1,no_guides=F)# + labs(x=NULL, y=NULL)
p2 <- plot_scatter(tSNE2_iscmeb_pois, meta_data, label_name = 'Domain', 
                   base_size=base_axis_size, cols=cols_clusters, point_size = pt_size_sample, 
                   point_alpha = 1,no_guides=F)
p12 <- drawFigs(list(p1, p2), layout.dim = c(2,1), common.legend = F)
ggsave(file="HCC4_iscmeb_pois_tSNE2.png", plot = p12,
       width = 6.5, height =10, units = "in", dpi = 200)

##### 4d ###########

load('timeVec_HCC4_new.rds')
names(timeVec)[names(timeVec)=='SpatialPCA-A'] <- 'SpatialPCA-L'
name_choose <-names(timeVec) #c('ProFAST_gauss', 'ProFAST_pois1','DR-SC',  'SpatialPCA_fast',   'SpatialPCA_raw' )

name_use <- c("ProFAST-P" ,    "ProFAST-G" ,    "SpatialPCA-L", "SpatialPCA-F",  "PRECAST" )
name_choose <-name_use

level_time <- names(col_orderd3)
df <- data.frame(Method= factor(name_choose, levels=level_time), Time= timeVec[name_choose])
p2 <-ggplot(data=df, aes(x=Method, y=Time,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd3[level_time])+
  theme_classic() + scale_x_discrete(breaks=NULL) + theme(text=element_text(size=20))+ ylab("Time (sec.)") +
  xlab(NULL) + scale_y_log10()+ theme(legend.position = 'none')
p2

load("memListGB_HCC4.rds")
memVec <- unlist(memListGB)
new_names <- rename_vec(names(memVec), oldname = "SpatialPCA-A", newname = "SpatialPCA-L")
names(memVec) <-new_names
name_use <- names(col_orderd3)
### choosed methods:
name_use <- c("ProFAST-P" ,    "ProFAST-G" ,  "SpatialPCA-F",  "SpatialPCA-L",   "PRECAST" )
df <- data.frame(Method= factor(name_use, levels=name_use), Memory= memVec[name_use]) # +1

p_mem_cho <- ggplot(data=df, aes(x=Method, y=Memory,  fill=Method)) +
  geom_bar(stat="identity", color='black')+ scale_fill_manual(values=col_orderd3[name_use])+
  theme_classic() + scale_x_discrete(breaks=NULL) + theme(text=element_text(size=20))+ ylab("Memory (Gigabytes)")+
  xlab(NULL) + ylab(NULL) + theme(legend.position = 'none')
drawFigs(pList=list(p2, p_mem_cho), layout.dim = c(1,2), common.legend = T, legend.position = 'bottom')

##### 4e ###########
load('seu_combine_hcc4.rds')
col_here <- c("#F2E6AB", "#9C0141") 
library(ggplot2)
p1 <- DotPlot(seu_combine_hcc4, features=row.names(seu_combine_hcc4), cols=col_here, #  idents = ident_here,
              col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic"))+
  ylab("Domain") + xlab(NULL) +theme(axis.text.y = element_text(face = "italic"),axis.text.x=element_text(size=16))

write_fig(p1, filename = "DotPlot_profastP_combine",
          width = 13, height = 8.1)

##### 4f ###########
load('Prop_alt_count_pois_version_newK9.rds')
library(ggplot2)
# create vectors of different length
vec1 <- Prop_alt_count[1:3]
vec2 <- Prop_alt_count[4:6]
vec3 <- Prop_alt_count[7:9]
cols_use_two <- c("#9C0141", "#F47646", "#757CCE" )
simutool::colorbar_adj_transparent(cols_use_two, alpha=1)
names(cols_use_two) <- c('TNE1', "TNE2", "Stroma")
# combine vectors into a data frame
data <- data.frame(value = c(vec1, vec2, vec3),
                   group = factor(rep(names(cols_use_two), each=3), levels = names(cols_use_two)) ) 

# plot boxplot
ggplot(data, aes(x = group, y = value, fill=group)) + geom_boxplot(width = 0.4) +
  labs(x = "Vector", y = "Value") + scale_fill_manual(values=cols_use_two)+
  labs(x = "", y = "Mutation count per spot") +
  theme_classic() + theme(text = element_text(size=20))

##### 4g ###########
load('seu_mutate_list_gene5.rds')
gene_mutate_list <- c("CERS2", "ETS2", "OCEL1", "RIF1", "IGSF23")
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))*2-1}
featurePlot <- function(seu, feature, cols, quant=0.2, pt_size=1,title_size=12, to_upper=TRUE){
  require(ggplot2)
  dat <- data.frame(Spatial_1=seu$row, Spatial_2=seu$col)
  dat$Expression <- seu[['RNA']]@scale.data[feature,]
  med <- quantile(seu[['RNA']]@scale.data[feature,], probs =quant)
  if(to_upper){
    fe <- toupper(feature)
  }else{
    fe <- feature
  }
  ggplot(data=dat, aes(x=Spatial_1, y=Spatial_2, color=Expression)) + geom_point(size=pt_size) +
    scale_colour_gradient2(
      low = cols[1],
      mid = "white",
      high = cols[2], midpoint = med) + mytheme_graybox() + 
    ggtitle(fe) + theme(title =element_text(size=title_size, color=1)) + xlim(c(0,77)) + ylim(c(1,120))
  
}

i <- 3
seu <- seu_mutate_list[[i]]
library(Seurat)
seu <- NormalizeData(seu)
scale.data <- pbapply::pbapply(seu[["RNA"]]@data, 1, range01)
seu[["RNA"]]@scale.data <- t(scale.data)
seu_mutate_list[[i]] <- seu

cols <-c("#361A95", "#D62728")
p31 <- featurePlot(seu_mutate_list[[i]], feature='CERS2', cols=cols, quant= 0.3,
                            pt_size = 1, title_size=12)

p32 <- featurePlot(seu_mutate_list[[i]], feature='ETS2', cols=cols, quant= 0.3,
                  pt_size = 1, title_size=12)

drawFigs(list(p31,p32), layout.dim = c(1,2))

##### 4h ###########

load('mIarray_HCC4.rds')
load("moranI_array_HCC4.rds")
pList <- list()
for(r in 1:4){
  # r <- 1
  mMat1 <- mIarray[r,,]
  colnames(mMat1) <- c(attr(moranI_array,"MethodNames"), "scVI", "PRECAST")
  method_use <- setdiff(names(col_orderd2), c("SpatialPCA-F",  "SpatialPCA-A"))
  pList[[r]] <- PRECAST::boxPlot(mMat1[,method_use], cols=col_orderd2[method_use], ylabel = "Moran's I")+
    theme_classic(base_size = 20)
  
}

drawFigs(pList[1:2], layout.dim = c(1,2), common.legend = T, legend.position = 'bottom')
