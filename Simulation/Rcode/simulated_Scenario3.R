######### We simulated three sections  based on a single section of the human DLPFC Visium dataset (sample ID: 151672) with eight spatial domains
### Scenario 3 (batch effect=high, biological effect= high)
# setwd("E:/Research paper/IntTemporalSpatial/AnalysisCode/ProFAST_Analysis/")

rm(list=ls())
source("./Real_data_results/Rcode/util_funcs.R")
load("./Simulation/simulated_datasets/batch_facScale0.6_de_prop0.3_seulist.rds")

library(ProFAST)
library(Seurat)

yList <- lapply(seulist, function(x) as.character(x$Group))
posList1 <- list()
AdjList <- list()
## To differetiate different sections by adding a large value for the coordinates
## since SpatialPCA and DR-SC, using the spatial coordinates, are designed to handle one section.
for(m in 1: length(seulist)){
  message("m = ", m)
  pos <- cbind(seulist[[m]]$row, seulist[[m]]$col)
  AdjList[[m]] <- PRECAST::getAdj_reg(pos, platform = "Visium")
  posList1[[m]] <- pos + (m-1)*1e6
}


seulist <-  pbapply::pblapply(seulist, NormalizeData)
XList_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@data))
nvec <- sapply(seulist, ncol)


library(ProFAST)
hq <- 15

## SpaFactor-G
tic <- proc.time()
reslist_fac_all <- ProFAST_run(XList_sp, AdjList, q=hq)
toc <- proc.time()
time_SpaMFactor <- toc[3] - tic[3]
(R2_fac_all <- evaluate_DR_PF2(reslist_fac_all$hV, yList))
# [1] 0.9319603 0.9321006 0.9503799

XList_count_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@counts))
## SpaFactor-P
tic <- proc.time()
reslist_fac_pois <- ProFAST_run(XList_count_sp, AdjList, q=hq, fit.model = 'poisson')
toc <- proc.time()
time_SpaMFactor_pois <- toc[3] - tic[3]
(R2_fac_pois <- evaluate_DR_PF2(reslist_fac_pois$hV, yList))
# [1] 0.9696653 0.9699204 0.9828730

### PCA
XList_all <- lapply(XList_sp, as.matrix)
Xmat <- matlist2mat(XList_all)
tic <- proc.time()
princ <- DR.SC:::wpca(Xmat, q= hq, weighted = F)
hV_pcaList <- mat2list(princ$PCs, nvec=nvec)
toc <- proc.time()
time_pca <- toc[3] - tic[3]
(R2_pca <- evaluate_DR_PF2(hV_pcaList, yList))
## [1] 0.2693746 0.2648712 0.2820139
### multibatchPCA
require(batchelor)
tic_bPCA <- proc.time()
res_bPCA <- multiBatchPCA(lapply(XList_all, t), d= hq)
toc_bPCA <- proc.time()
time_bPCA <- toc_bPCA[3] - tic_bPCA[3]
(R2_bPCA <- evaluate_DR_PF2(res_bPCA@listData, yList))
## [1] 0.2671312 0.2642469 0.2829737

