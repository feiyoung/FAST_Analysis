# Load data and preprocessing ---------------------------------------------

source("./util_funcs.R")
library(Seurat)
library(data.table)

seuList <- list()
for(isample in 1:2){
  # isample <- 2
  message("isample = ", isample)
  data = Read10X(paste0("./hBreast_ffpe_rep", isample, "/cell_feature_matrix_mtx"))
  rna = data[[1]]
  meta = fread(paste0("./hBreast_ffpe_rep",isample, "/cell_info.csv.gz"))
  
  all(colnames(rna) == meta$cell_id)
  rownames(meta) = meta$cell_id
  meta$cell_id <- NULL
  seu <- CreateSeuratObject(counts = rna, meta.data=meta)
  
  seuList[[isample]] <- seu
}
geneNames <- intersect(row.names(seuList[[1]]), row.names(seuList[[2]]))
### weak quality control because the number of genes are only 313.
seulist <- pbapply::pblapply(seuList, filter_spot, min_feature=15) 
num_gene_spot_filter <- sapply(seulist, dim)

seulist <- lapply(seulist, NormalizeData)

posList <- list()
posList1 <- list()
AdjList <- list()
for(m in 1: length(seulist)){
  # m <- 1
  message("m = ", m)
  pos <- cbind(seulist[[m]]$x_centroid, seulist[[m]]$y_centroid)
  posList1[[m]] <- pos + (m-1)*1e6 # to differentiate different sections
  posList[[m]] <- pos 
  Adj <- DR.SC::getAdj_auto(pos, radius.upper = 50)
  AdjList[[m]] <- Adj
}



# Spatial dimension reduction using ProFAST -------------------------------
## The previous name of FAST is called ProFAST
library(ProFAST)

XList_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@data)) 
XList <- lapply(XList_sp, as.matrix)
n_vec <- sapply(XList, nrow)
library("profmem")
options(profmem.threshold = 2000)
hq <- 15
tic <- proc.time()
pm_profastG <- profmem({
  reslist <- ProFAST_run(XList_sp, AdjList, q=hq, maxIter=25)
})
pm_profastG
total_mem_profastG <- sum(pm_profastG$bytes, na.rm = T)
save(pm_profastG, total_mem_profastG, file='total_mem_profastG_xeHBC3.rds')
toc <- proc.time()
time_ProFAST <- toc[3] - tic[3] 

save(reslist, time_ProFAST, file='reslist_ProFASTG_xeHBC2.rds')


XList_count_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@counts))
tic <- proc.time()
pm_profastP <- profmem({
  reslist_hvg_count <- ProFAST_run(XList_count_sp, AdjList,  q=hq, fit.model = 'poisson') 
})
total_mem_profastP <- sum(pm_profastP$bytes, na.rm = T)
save(pm_profastP, total_mem_profastP, file='total_mem_profastP_xeHBC3.rds')
toc <- proc.time()
time_ProFAST_hvg_count <- toc[3] - tic[3]
save(reslist_hvg_count, time_ProFAST_hvg_count, file ='reslist_ProFASTP_xeHBC2.rds' )


# Downstream analysis -----------------------------------------------------
# Embedding alignment and spatial clustering ------------------------------
## output the function in ProFAST because fit.iscmeb() is now embeded in the R package but not exported.
fit.iscmeb <- ProFAST:::fit.iscmeb

sampleID <- get_sampleID(reslist$hV)
library(harmony)
### Choose number of clusters
tic <- proc.time()
hZ_harmony_profastP <- HarmonyMatrix(matlist2mat(reslist_hvg_count$hV), meta_data = data.frame(sample = sampleID),
                                     vars_use = "sample", do_pca = F)
toc <- proc.time()
res_louvain_harmony_profastP <- drLouvain(hZ_harmony_profastP, resolution = 0.5)

hK_pois <- length(unique(res_louvain_harmony_profastP))
tic <- proc.time()
reslist_iscmeb_pois <- fit.iscmeb(
  reslist_hvg_count$hV,
  AdjList,
  K=hK_pois,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 30,
  init.start=1,
  epsLogLik = 1e-05,
  coreNum = 1,
  verbose = TRUE)
toc <- proc.time()
time_iscmeb_pois <- toc[3] - tic[3]
save(reslist_iscmeb_pois, time_iscmeb_pois, file='reslist_iscmeb_profastP_xeHBC2.rds')
hyList_iscmeb_pois <- reslist_iscmeb_pois@idents
save(hyList_iscmeb_pois, file='hyList_iscmeb_profastP_xeHBC2.rds')

# Re-order the cluster number

### Poisson versions:
load('hyList_iscmeb_profastP_xeHBC2.rds')
rawID <- c(2,4,9,7,10, 5,6,8,13,   1, 15, 12, 11, 3,14, 16, 17)
names(rawID) <- 1:17
hyList_iscmeb_pois_renumber <- lapply(hyList_iscmeb_pois, replace_ID, rawID=rawID)
lapply(hyList_iscmeb_pois_renumber, table)
save(hyList_iscmeb_pois_renumber, file='hyList_iscmeb_profastP_renumber_xeHBC2.rds')


# Removing unwanted variation in multi-section SRT data -------------------
library(PRECAST)
data("Human_HK_genes")
gene_intersect <- geneNames
idx_hk <- which(gene_intersect %in% as.character(Human_HK_genes$Gene))
housekeep_genes <- gene_intersect[idx_hk]
save(housekeep_genes, file='housekeep_genes_xeHBC2.rds')

XList_hk <- pbapply::pblapply(seulist, function(x) Matrix::t(x[['RNA']]@data[housekeep_genes,]))
nvec <- sapply(XList, nrow)
HList <- lapply(XList_hk, as.matrix) ## use these six housekeeping genes as negative control.
HList <- lapply(HList, scale, scale=F)
for(r in 1:2) colnames(HList[[r]]) <- paste0(1:ncol(HList[[1]]))

M <- length(XList)
p <- ncol(XList[[1]])
Tm <- matrix(0, M, 1)
RList_pois <- reslist_iscmeb_pois@fitList[[1]]$Rf

correct_genesR <- ProFAST:::correct_genesR ## output the function in ProFAST because correct_genesR() is now embeded in the R package.
tic <- proc.time()
correct_List_pois <- correct_genesR(XList, RList_pois, HList, Tm, AdjList)
toc <- proc.time()
time_all_pois <- toc[3] - tic[3] #

XList_correct_all_pois <- correct_List_pois$XList_correct
## Add row names and colnames (symbol) for each slice
for(m in 1:M){
  row.names( XList_correct_all_pois[[m]]) <- paste0("slice", m, "_", row.names( XList[[m]]) )
}
hX_pois <- matlist2mat(XList_correct_all_pois)
library(Seurat)
library(Matrix)
count <- sparseMatrix(i=1,j=1, x=1, dims=dim(t(hX_pois)))
row.names(count) <- colnames(hX_pois)
colnames(count) <- row.names(hX_pois)
seuAll_count <- CreateSeuratObject(counts = count)
length(unique(colnames(seuAll_count)))
row.names(hX_pois) <- colnames(seuAll_count)
seuAll_count[["RNA"]]@data <- t(hX_pois)
load("hyList_iscmeb_profastP_renumber_xeHBC2.rds")
Idents(seuAll_count) <-  factor(unlist(hyList_iscmeb_pois_renumber), levels=1:17)


# Combined differential expression analysis -------------------------------
dat_degs_pois <- FindAllMarkers(seuAll_count)
save(dat_degs_pois, file='housekeep_dat_degs_choseK17_iscmeb_profastP_xeHBC2.rds')
cutoff <- 0.001#0.05
dim(subset(dat_degs_pois,  p_val_adj<cutoff & avg_log2FC>0.25))
filename <- 'xeHBC2_DEGsList_iSCMEB_poisK17.xlsx'
K_use <- length(unique(dat_degs_pois$cluster))
for(k in 1:K_use){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_pois,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  ## sort by log fold change
  dat_degs_sub3 <- dat_degs_sub3[order(dat_degs_sub3$avg_log2FC,decreasing = T),]
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Cluster', k),
                   append = T,row.names = F)
}


# KEGG pathway analysis ---------------------------------------------------

library(gprofiler2)
table(dat_degs_pois$cluster)
termList_cluster <- list()
for(k in 1: 17){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_pois,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  que1 <- toupper(dat_degs_sub3$gene)
  gostres <- gost(query = que1,
                  organism = "hsapiens", correction_method='fdr', user_threshold =1)
  termList_cluster[[k]] <- gostres
}
save(termList_cluster, file='profiler_termListPoisClusterDEGs_xeHBC2.rds')

get_top_pathway <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id", "term_size","p_value")])
  }
  
  return(df_sub)
}

### Focus on KEGG pathways
dat_kegg <- NULL
source_set <- 'KEGG'
for(ii in 1:length(termList_cluster)){
  ## ii <- 1
  message("ii=", ii)
  gostres2 <- termList_cluster[[ii]]
  max(gostres2$result$term_size)
  dat_tmp <- subset(gostres2$result, term_size<500 & p_value<0.05)
  table(dat_tmp$source)
  dat1 <- get_top_pathway(dat_tmp, ntop=5, source_set = source_set)
  dat1$nlog10P <- -log10(dat1$p_value)
  dat1_sub <- subset(dat1[order(dat1$nlog10P),],  nchar(term_name)< 60)
  if(nrow(dat1_sub)>0){
    dat1_sub$Cluster <- paste0("Cluster", ii)
    dat_kegg <- rbind(dat_kegg, dat1_sub)
  }
  
}
KEGG_id <- dat_kegg$term_id

## Take out all these KEGG pathways for each cluster
dat_kegg_all <- NULL
for(ii in 1:length(termList_cluster)){
  ## ii <- 5
  message("ii=", ii)
  gostres2 <- termList_cluster[[ii]]
  max(gostres2$result$term_size)
  dat_tmp <- subset(gostres2$result, term_size<500& p_value<0.05)
  table(dat_tmp$source)
  dat1 <- subset(dat_tmp, term_id %in% KEGG_id)
  # dat1 <- rbind(dat1, subset(dat_tmp[,c("source", "term_name", "term_id","p_value")], term_id %in% extra_term_idList[[ii]]))
  dat1$nlog10P <- -log10(dat1$p_value)
  dat1_sub <- subset(dat1[order(dat1$nlog10P),],  nchar(term_name)< 60)
  if(nrow(dat1_sub)>0){
    dat1_sub$Cluster <- ii
    dat_kegg_all <- rbind(dat_kegg_all, dat1_sub)
  }
}




#  Cell type deconvolution analysis ----------------------------------------------------------------------

library(Seurat)
expdat <- ReadMtx('count_matrix_sparse.mtx', cells='count_matrix_barcodes.tsv',
                  features = 'count_matrix_genes.tsv', feature.column = 1)
library(data.table)
metadat <- fread("metadata.csv")
table(metadat$celltype_subset) # 49 cluster
table(metadat$celltype_major)  # 9 clusters
metadat <- as.data.frame(metadat)
row.names(metadat) <- metadat$V1
sc_data = CreateSeuratObject(counts = expdat,
                             meta.data = metadat)
head(sc_data)
######### run RCTD
output_dir = 'RCTD_deconvolution_result'

cell_type_columns = 'celltype_subset'

for (i in 1: length(seulist)){
  # i <- 1
  message("i = ", i)
  st_s = seulist[[i]]
  output_name = 'GSE176078_st'
  st_count = st_s@assays$RNA@counts
  st_coord = st_s@meta.data[,c('x_centroid','y_centroid')]
  output_name = paste(output_name,i,sep="")
  RCTD_run(sc_data,cell_type_columns,st_count,st_coord,output_dir,output_name)   
}



# Transcription factor analysis -------------------------------------------
run_viper_Seurat <- function(input, regulons, options = list(), tidy = FALSE, assay_key = "RNA") {
  if (tidy) {
    tidy <- FALSE
    warning("The argument 'tidy' cannot be TRUE for 'Seurat' objects. ",
            "'tidy' is set to FALSE")
  }
  mat <- as.matrix(Seurat::GetAssayData(input, assay = assay_key, slot = "data"))
  
  tf_activities <- do.call(viper::viper,
                           c(list(eset = mat, regulon = regulons), options))
  
  # include TF activities in Seurat object
  dorothea_assay <- Seurat::CreateAssayObject(data = tf_activities)
  Seurat::Key(dorothea_assay) <- "dorothea_"
  input[["dorothea"]] <- dorothea_assay
  
  return(input)
}

my_run_viper <- function(seulist, network, path, plot_filename_prefix, filename_viper_res, filename_viper_res_tf) {
  # run viper
  nT = length(seulist)
  viper_res <- lapply(1:nT, function(i) {
    run_viper_Seurat(seulist[[i]], network, assay_key = "RNA")
  })
  
  # violin plot
  setwd(path)
  lapply(1:nT, function(i) {
    seu <- viper_res[[i]]
    features = rownames(seu[["dorothea"]])
    plot_ncol = ceiling(sqrt(length(features)))
    p0 = VlnPlot(seu, group.by = "cluster", assay = "dorothea", features = features, pt.size = 0.1, ncol = plot_ncol) + 
      NoLegend()
    ggsave(file=paste0(plot_filename_prefix, i, ".png"), plot = p0, width = 12, height = 12)
  })
  
  save(viper_res, file = paste0(filename_viper_res, ".rda"))
  
  # active tfs
  viper_res_tf = lapply(1:2, function(i) {
    features = rownames(viper_res[[i]][["dorothea"]])
    network[features]
  })
  
  save(viper_res_tf, file = paste0(filename_viper_res_tf, ".rda"))
}


library(SeuratDisk)
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(minet)
library(viper)

geneNames <- intersect(row.names(seuList[[1]]), row.names(seuList[[2]]))
seulist <- pbapply::pblapply(seuList, filter_spot, min_feature=15) 
seulist <- lapply(1:2, function(i) {
  seu <- seulist[[i]] %>% 
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData() %>%
    AddMetaData(metadata = as.integer(hyList_iscmeb_pois_renumber[[i]]), col.name = "cluster")
  seu
})
library(aracne.networks)
data(regulonbrca)
require('org.Hs.eg.db')
# 'org.Hs.eg.db' for human
# 'org.Mm.eg.db' for mouse
regulonbrca = lapply(regulonbrca, function(obj) {
  tmp <- names(obj$tfmode)
  names(obj$tfmode) <- GeneAnswers::getSymbols(tmp, 'org.Hs.eg.db')
  obj
})
tmp <- names(regulonbrca)
names(regulonbrca) <- GeneAnswers::getSymbols(tmp, 'org.Hs.eg.db')
# run viper with regulonbrca

path = "./"
plot_filename_prefix = "viper_poisson_ARACNe"
filename_viper_res = "viper_poisson_ARACNe_res"
filename_viper_res_tf = "viper_poisson_ARACNe_res_tf"
my_run_viper(seulist, regulonbrca, path, plot_filename_prefix, filename_viper_res, filename_viper_res_tf)

