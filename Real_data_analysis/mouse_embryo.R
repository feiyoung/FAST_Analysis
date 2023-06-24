
# Load data and preprocessing ---------------------------------------------

library(PRECAST)
library(Seurat)
source("./util_funcs.R")

###### Start analyzing
file_vec <- c('E12.5_E1S1.MOSTA', 'E12.5_E1S2.MOSTA', 'E12.5_E1S3.MOSTA' ,'E12.5_E1S4.MOSTA',  'E12.5_E1S5.MOSTA', 'E12.5_E2S1.MOSTA',
              'E13.5_E1S1.MOSTA', 'E13.5_E1S2.MOSTA', 'E13.5_E1S3.MOSTA', 'E13.5_E1S4.MOSTA',
              'E14.5_E1S1.MOSTA', 'E14.5_E1S2.MOSTA','E14.5_E1S3.MOSTA', 'E14.5_E1S4.MOSTA', 'E14.5_E1S5.MOSTA', 'E14.5_E2S1.MOSTA',
              'E14.5_E2S2.MOSTA', 'E15.5_E1S2.MOSTA','E15.5_E1S3.MOSTA', 'E15.5_E1S4.MOSTA', 'E15.5_E2S1.MOSTA', 'E15.5_E1S1.MOSTA',
              "E16.5_E1S1.MOSTA", "E16.5_E1S2.MOSTA", "E16.5_E1S3.MOSTA", "E16.5_E1S4.MOSTA")

sample_names <- sub(".MOSTA", "", file_vec)
seuList_filter <- list()
raw_numberList <- list()
for(m in 1: length(file_vec)){
  # m <- 1
  message("m = ", m)
  dir_name <- file_vec[m]
  seu_tmp <- read_embryo_data(dir_name)
  sp_count <- seu_tmp[["RNA"]]@counts 
  raw_numberList[[m]] <- dim(sp_count)
  ## Filter counts:
  min_spots <- 20
  gene_flag <- rowSums(sp_count>0)>min_spots
  sp_count_filter <-  sp_count[gene_flag, ]
  min_genes <- 20
  spot_flag <- colSums(sp_count>0) > min_genes
  sp_count_filter <- sp_count_filter[, spot_flag]
  seuList_filter[[m]] <- seu_tmp[row.names(sp_count_filter), colnames(sp_count_filter)]
}

save(raw_numberList, file='raw_numberList_E12E16.rds')
raw_numMat <- Reduce(cbind, raw_numberList)
rowMeans(raw_numMat)
rowSums(raw_numMat)

M <- length(seuList_filter)
HVGList <- list(); timeVec_HVGs <- rep(NA, M)
for(m in 1:M){
  # m <- 1
  tic <- proc.time()
  seu <- seuList_filter[[m]]
  seu <-  NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  toc <- proc.time()
  var.features <- seu@assays$RNA@var.features
  HVGList[[m]] <- var.features
  timeVec_HVGs[m] <- toc[3] - tic[3]
}


save(timeVec_HVGs, HVGList, file = paste0("E12E16_HVGList.rds"))

gene_hvg_2000 <- selectIntFeatures(seuList_filter,HVGList) # 
save(gene_hvg_2000, file='gene_hvg_2000_E12E16.rds')

## Get the barcodes for each slides
barcodelist <- lapply(seulist, colnames)
save(barcodelist, file='barcodelist_E12E16.rds')
load('barcodelist_E12E16.rds')
seuList_filter <-  pbapply::pblapply(1:length(barcodelist), function(j) seuList_filter[[j]][, barcodelist[[j]]] )
gene_intersect <- Reduce(intersect, lapply(seuList_filter, row.names))
library(PRECAST)
data("Mouse_HK_genes")
idx_hk <- which(gene_intersect %in% as.character(Mouse_HK_genes$Gene))
housekeep_genes <- gene_intersect[idx_hk]

seulist_HK <- pbapply::pblapply(seuList_filter, function(x) x[(housekeep_genes),])
saveRDS(seulist_HK, file='seulist_HK.RDS')

## Use SPARK-X to select the gene without time-varying values for all slices
seulist_HK <- pbapply::pblapply(seulist_HK, DR.SC::FindSVGs, nfeatures = nrow(seulist_HK[[1]]))
geneList_noSpa <- lapply(seulist_HK, function(x) {
  
  order_idx <- order(x[['RNA']]@meta.features$adjusted.pval.SVGs, decreasing = T) # rank the gene with largest p-value the first
  genes <- row.names(x[['RNA']]@meta.features)[order_idx]
  return(genes[1:200])
}) ## return a list with ordered gene vector
num_genes_eachbatch <- sapply(geneList_noSpa, length)

HKfeatures <- selectHKFeatures(seulist_HK, geneList_noSpa, HKFeatures = 200)
HKsymbols <- housekeep_genes[(housekeep_genes) %in% HKfeatures]
save(HKsymbols, file='housekeep_noSpa_genes_E12E16.rds')



yList <- lapply(seulist, function(x) x$annotation)
posList <- list()
for(m in 1: length(seulist)){
  message("m = ", m)
  pos <- cbind(seulist[[m]]$row, seulist[[m]]$col)
  posList[[m]] <- pos + (m-1)*1e6
}
save(yList, file='yList_E12E16.rds')
save(posList, file='posList_E12E16.rds')

AdjList <- list()
for(m in 1: length(posList)){
  message("m = ", m)
  pos <-  posList[[m]] 
  AdjList[[m]] <- DR.SC::getAdj_auto(pos, radius.upper = 4)
}


# Spatial dimension reduction using ProFAST -------------------------------
library(ProFAST)
XList_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@data)) # save XList, XList_count, posList, yList, basic information.
XList <- lapply(XList_sp, as.matrix)
n_vec <- sapply(XList, nrow)
hq <- 15
tic <- proc.time()
reslist <- ProFAST_run(XList, AdjList, q=hq, maxIter=25)
toc <- proc.time()
time_ProFAST <- toc[3] - tic[3] 
save(reslist, time_ProFAST, file='reslist_ProFAST_gauss_E12E16.rds')
sum(sapply(reslist$hV, nrow))
(R2_pro_hvg_log <- evaluate_DR_PF(reslist$hV, yList))
save(R2_pro_hvg_log, file='MacR2_spaFactor_gauss_E12E16.rds')


# Downstream analysis -----------------------------------------------------
# Embedding alignment and spatial clustering ------------------------------
## output the function in ProFAST because fit.iscmeb() is now embeded in the R package but not exported.
fit.iscmeb <- ProFAST:::fit.iscmeb


sampleID <- get_sampleID(reslist$hV)
library(harmony)
### Choose number of clusters
tic <- proc.time()
hZ_harmony_profastG <- HarmonyMatrix(matlist2mat(reslist$hV), meta_data = data.frame(sample = sampleID),
                                     vars_use = "sample", do_pca = F)
toc <- proc.time()
res_louvain_harmony_profastG <- drLouvain(hZ_harmony_profastG, resolution = 0.2)

hK<- length(unique(res_louvain_harmony_profastG))
tic <- proc.time()
reslist_iscmeb_profastG <- fit.iscmeb(
  reslist$hV,
  AdjList,
  K=hK,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 25,
  init.start=1,
  seed = 1,
  epsLogLik = 1e-05,
  coreNum = 1,
  verbose = TRUE)
toc <- proc.time()
save(reslist_iscmeb_profastG, file='reslist_iscmeb_profastG_E12E16.rds')

# Removing unwanted variation in multi-section SRT data -------------------
library(PRECAST)
library(Seurat)
load('housekeep_noSpa_genes_E12E16.rds')
XList_hk <- pbapply::pblapply(seulist_HK, function(x) x[['RNA']]@data[HKsymbols,])
nvec <- sapply(XList_hk, ncol)
all(nvec == sapply(posList, nrow))
XList_hk <- pbapply::pblapply(XList_hk, function(x) base::t(as.matrix(x)))
Xmat_hk <- matlist2mat(XList_hk)
princ <- approxPCA(Xmat_hk, q=10)
HList <- mat2list(princ$PCs, nvec)
rm(Xmat_hk)
RList <-  reslist_iscmeb_profastG@fitList[[1]]$Rf
all(sapply(RList, nrow) == sapply(HList, nrow))
save(RList, file='RList_for_batch_correct_genes_E12E16.rds')

load('RList_for_batch_correct_genes_E12E16.rds')
load('HList_for_batch_correct_genes_E12E16.rds')
XList <- lapply(XList_sp, as.matrix)
M <- length(XList)
p <- ncol(XList[[1]])
E_list <- list(5:10, 11:14, 15:21, 22: 26, 27:30)
Tm <- NULL
for(r in 1:length(E_list)){
  Tm <- c(Tm, rep(12.5+(r-1), length(E_list[[r]])))
}
Tm <- matrix(Tm, length(Tm), 1)

correct_genes_subsampleR <- ProFAST:::correct_genes_subsampleR ## output the function in ProFAST because correct_genes_subsampleR() is now embeded in the R package.
tic <- proc.time()
correct_List_pois <- correct_genes_subsampleR(XList, RList, HList, Tm, AdjList, subsample_rate = 0.1)
toc <- proc.time()
time_all_pois <- toc[3] - tic[3] #

XList_correct_all_pois <- correct_List_pois$XList_correct

hX <- matlist2mat(XList_correct_all)
library(Seurat)
library(Matrix)
count <- sparseMatrix(i=1,j=1, x=1, dims=dim(t(hX)))
row.names(count) <- colnames(hX)
colnames(count) <- row.names(hX)
seuAll <- CreateSeuratObject(counts = count)
length(unique(colnames(seuAll)))
row.names(hX) <- colnames(seuAll)
seuAll[["RNA"]]@data <- t(hX)
Idents(seuAll) <- factor(unlist(clusterList_profastG), levels=1:20)
saveRDS(seuAll, file='seuAll_gauss_E12E16.RDS')

# Combined differential expression analysis -------------------------------
dat_degs <- FindAllMarkers(seuAll, only.pos = T)

save(dat_degs, file='dat_degs_choseK20_iscmeb_gauss_E12E16.rds')
cutoff <- 0.001#0.05
dim(subset(dat_degs,  p_val_adj<cutoff & avg_log2FC>0.25))
filename <- 'E12E16_DEGsList_iSCMEB_gaussK20_reorder.xlsx'
K_use <- length(unique(dat_degs$cluster))
for(k in 1:K_use){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  ## sort by log fold change
  dat_degs_sub3 <- dat_degs_sub3[order(dat_degs_sub3$avg_log2FC,decreasing = T),]
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Cluster', k),
                   append = T,row.names = F)
}

# Detect the temporal genes using linear regression -----------------------
### For brain
idx_brain <- which(Idents(seuAll) %in% c(1:6))
seu_brain <- seuAll[,idx_brain]
## linear regression
regMat <- matrix(NA, nrow=nrow(seu_brain), ncol=3)
row.names(regMat) <- row.names(seu_brain)
colnames(regMat) <- c("R2", "adj_R2", "pval")
for(j_gene in 1: nrow(seu_brain)){
  #j_gene <- 1
  message("j_gene = ", j_gene)
  lm1 <- lm(seu_brain[["RNA"]]@data[j_gene,] ~ seu_brain@meta.data$Tm_num)
  dat <- summary(lm1)
  dat$coefficients
  regMat[j_gene, ] <-  c(R2=dat$r.squared, adj_R2= dat$adj.r.squared, pval= dat$coefficients[2,4])
  
}

regMat <- regMat[order(regMat[,1], decreasing = T), ]
save(regMat, file='regMat_Brain_E12E16.rds')
