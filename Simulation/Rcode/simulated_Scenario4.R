######### We simulated three sections  based on a single section of the human DLPFC Visium dataset (sample ID: 151672) with eight spatial domains
### Scenario 4 (batch effect=high, biological effect= low): the most difficult scenario
rm(list=ls())

load("./Simulation/simulated_datasets/batch_facScale0.6_de_prop0.15_seulist.rds")
source("./Real_data_results/Rcode/util_funcs.R")

library(FAST)
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


library(FAST)
hq <- 15

## SpaFactor-G
tic <- proc.time()
reslist_fac_all <- FAST_run(XList_sp, AdjList, q=hq)
toc <- proc.time()
time_SpaMFactor <- toc[3] - tic[3]
(R2_fac_all <- evaluate_DR_PF2(reslist_fac_all$hV, yList))
# [1] 0.4952630 0.5624238 0.5270935

XList_count_sp <- lapply(seulist, function(x) Matrix::t(x[["RNA"]]@counts))
## SpaFactor-P
tic <- proc.time()
reslist_fac_pois <- FAST_run(XList_count_sp, AdjList, q=hq, fit.model = 'poisson')
toc <- proc.time()
time_SpaMFactor_pois <- toc[3] - tic[3]
(R2_fac_pois <- evaluate_DR_PF2(reslist_fac_pois$hV, yList))
# [1] 0.8768253 0.9061791 0.9321143

### PCA
XList_all <- lapply(XList_sp, as.matrix)
Xmat <- matlist2mat(XList_all)
tic <- proc.time()
princ <- DR.SC:::wpca(Xmat, q= hq, weighted = F)
hV_pcaList <- mat2list(princ$PCs, nvec=nvec)
toc <- proc.time()
time_pca <- toc[3] - tic[3]
(R2_pca <- evaluate_DR_PF2(hV_pcaList, yList))
## [1] 0.03922998 0.04282111 0.03659799
### multibatchPCA
require(batchelor)
tic_bPCA <- proc.time()
res_bPCA <- multiBatchPCA(lapply(XList_all, t), d= hq)
toc_bPCA <- proc.time()
time_bPCA <- toc_bPCA[3] - tic_bPCA[3]
(R2_bPCA <- evaluate_DR_PF2(res_bPCA@listData, yList))
## [1] 0.04013854 0.04414979 0.03809931



