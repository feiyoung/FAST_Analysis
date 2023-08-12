# Load data and preprocessing ---------------------------------------------

source("./util_funcs.R")
name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                            151671, 151672, 151673, 151674, 151675, 151676))
n_sample <- length(name_ID12)
# dir_file <- "./braindata/"
url_brainA <- "https://github.com/feiyoung/DR-SC.Analysis/raw/main/data/DLPFC_data/"; url_brainB <- ".rds"

library(DR.SC)
library(Seurat)
library(SingleCellExperiment)
HVGList <- list()
timeVec <- rep(NA, n_sample)
seuList_filter <- list()
for(r in 1: n_sample){
  # r <- 1
  message('r = ', r)
  id <- name_ID12[r]
  dlpfc <- readRDS(url(paste0(url_brainA, name_ID12[r],url_brainB) ))
  sp_count <- dlpfc@assays@data$counts
  meta.data <- as.data.frame(colData(dlpfc))
  ## Filter counts:
  min_spots <- 20
  gene_flag <- rowSums(sp_count>0)>min_spots
  sp_count_filter <-  sp_count[gene_flag, ]
  min_genes <- 20
  spot_flag <- colSums(sp_count>0) > min_genes
  sp_count_filter <- sp_count_filter[, spot_flag]
  
  tic <- proc.time()
  seu <- CreateSeuratObject(counts=sp_count_filter, meta.data = meta.data[spot_flag,])
  seu <-  NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  toc <- proc.time()
  var.features <- seu@assays$RNA@var.features
  seuList_filter[[r]] <- seu
  HVGList[[r]] <- var.features
  timeVec[r] <- toc[3] - tic[3]
}
save(timeVec, HVGList, file = paste0("DLPFC12_HVGList.rds"))

gene_intersect <- Reduce(intersect, lapply(seuList_filter, row.names))
## transfer to symbols
gene_symbols <- transferGeneNames(gene_intersect, species="Human")
gene_symbols <- toupper(gene_symbols)
names(gene_symbols) <- gene_intersect
save(gene_symbols, file='gene_symbols_DLPFC12.rds')

## select the housekeeping genes
library(PRECAST)
data("Human_HK_genes")
idx_hk <- which(gene_symbols %in% as.character(Human_HK_genes$Gene))
housekeep_genes <- gene_symbols[idx_hk]
seulist_HK <- pbapply::pblapply(seuList_filter, function(x) x[names(housekeep_genes),])
## Select the housekeeping genes without lowest spatial variation.
seulist_HK <- pbapply::pblapply(seulist_HK, DR.SC::FindSVGs, nfeatures = nrow(seulist_HK[[1]]))
geneList_noSpa <- lapply(seulist_HK, function(x) {
  
  order_idx <- order(x[['RNA']]@meta.features$adjusted.pval.SVGs, decreasing = T) # rank the gene with largest p-value the first
  genes <- row.names(x[['RNA']]@meta.features)[order_idx]
  # idx <- which(x[['RNA']]@meta.features$adjusted.pval.SVGs[order_idx]>0.05)
  return(genes[1:200])
}) ## return a list with ordered gene vector
num_genes_eachbatch <- sapply(geneList_noSpa, length)

HKfeatures <- selectHKFeatures(seulist_HK, geneList_noSpa, HKFeatures = 200) ## select top 200
HKsymbols <- housekeep_genes[names(housekeep_genes) %in% HKfeatures]
save(HKsymbols, file='housekeep_noSpa_top200genes_DLPFC12.rds')


## select the integrated HVGs
gene_hvg_2000 <- selectIntFeatures(seuList_filter,HVGList) # 
save(gene_hvg_2000, file='gene_hvg_2000_DLPFC12.rds')

seulist_hvg <- lapply(seuList_filter, function(x) x[gene_hvg_2000, ])
yList <- lapply(seulist_hvg, function(x) x$layer_guess_reordered)
for(r in 1:length(yList)){
  y_true <- as.character(yList[[r]])
  y_true[is.na(y_true)] <- "NA"
  yList[[r]] <- y_true 
}
posList <- list()
AdjList <- list()
## To differetiate different sections by adding a large value for the coordinates
## since SpatialPCA and DR-SC, used the spatial coordinates, are designed to handle one section.
for(m in 1: length(seulist_hvg)){
  message("m = ", m)
  pos <- cbind(seulist_hvg[[m]]$row, seulist_hvg[[m]]$col)
  AdjList[[m]] <- PRECAST::getAdj_reg(pos, platform = "Visium")
  posList[[m]] <- pos + (m-1)* 1e6 
  
}
save(yList, file='yList_DLPFC12.rds')
save(posList, file='posList_DLPFC12.rds')




# Spatial dimension reduction using ProFAST -------------------------------
## The previous name of FAST is called ProFAST

library(Seurat)
seulist_hvg <- lapply(seulist_hvg, NormalizeData)
XList1 <- lapply(seulist_hvg, function(x) Matrix::t(x[["RNA"]]@data))
sum(sapply(XList1, nrow))
XList <- lapply(XList1, as.matrix)

library(ProFAST)
hq <- 15

### Gaussian version 
tic <- proc.time()
pm_profastG <- profmem({
  reslist_profastG <- ProFAST_run(XList1, AdjList = AdjList, q=hq)
})
total_mem_profastG <- sum(pm_profastG$bytes, na.rm = T)
save(pm_profastG, total_mem_profastG, file='total_mem_profastG_DLPFC12.rds')

toc <- proc.time()
time_profastG <- toc[3] - tic[3]
save(reslist_profastG,time_profastG, file="reslist_hvg_profastG_DLPFC12.rds")
(R2_profastG <- evaluate_DR_PF2(reslist_profastG$hV, yList))
save(R2_profastG, file='MacR2_hvg_profastG_DLPFC12.rds')


### Poisson version
XList_count <- lapply(seulist_hvg, function(x) Matrix::t(x[["RNA"]]@counts))
tic <- proc.time()
pm_profastP <- profmem({
  reslist_profastP <- ProFAST_run(XList_count, AdjList = AdjList, fit.model = "poisson", 
                                  q=hq)
})
total_mem_profastP <- sum(pm_profastP$bytes, na.rm = T)
save(pm_profastP, total_mem_profastP, file='total_mem_profastP_DLPFC12.rds')
toc <- proc.time()
time_profastP <- toc[3] - tic[3]
save(reslist_profastP,time_profastP, file="reslist_hvg_profastP_DLPFC12.rds")
(R2_profastP <- evaluate_DR_PF2(reslist_profastP$hV, yList))




# Downstream analysis -----------------------------------------------------


# Embedding alignment and spatial clustering ------------------------------
## output the function in ProFAST because fit.iscmeb() is now embeded in the R package but not exported.
fit.iscmeb <- ProFAST:::fit.iscmeb

embed_uncorrectList <- list()
embed_uncorrectList[["ProFAST-G"]] <- reslist_profastG$hV
embed_uncorrectList[["ProFAST-P"]] <- reslist_profastP$hV
sampleID <- get_sampleID(posList)
hamonyallList <- list()
library(harmony)
methods2 <- c("ProFAST-G","ProFAST-P")
for(im in methods2){
  message("im = ", im)
  tic <- proc.time()
  PCs_all <- matlist2mat(embed_uncorrectList[[im]])
  hZ_harmony <- HarmonyMatrix(PCs_all, meta_data = data.frame(sample = sampleID),
                              vars_use = "sample", do_pca = F)
  toc <- proc.time()
  hamonyallList[[im]] <- hZ_harmony
}

### Determine the number of clusters for ProFAST-P
res_louvain_profastP <- drLouvain(hamonyallList[["ProFAST-P"]], resolution = 0.2)
res_louvain_profastG <- drLouvain(hamonyallList[["ProFAST-G"]], resolution = 0.2)

hK_profastP  <- length(unique(res_louvain_profastP))
tic <- proc.time()
reslist_iscmeb_profastP <- fit.iscmeb(
  reslist_profastP$hV,
  AdjList,
  K=hK_profastP,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 25, # 25
  epsLogLik = 1e-05,
  verbose = TRUE)
toc <- proc.time()
time_iscmeb_profastP <- toc[3] - tic[3]
save(reslist_iscmeb_profastP,time_iscmeb_profastP, file='reslist_iscmeb_profastP_DLPFC12.rds')

hK_profastG <- length(unique(res_louvain_profastG)) 
tic <- proc.time()
reslist_iscmeb_profastG <- fit.iscmeb(
  reslist_profastG$hV,
  AdjList,
  K=hK_profastG,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 25,
  init.start=5,
  epsLogLik = 1e-05,
  coreNum = 1,
  verbose = TRUE)
toc <- proc.time()
time_iscmeb_profastG <- toc[3] - tic[3]
save(reslist_iscmeb_profastG, time_iscmeb_profastG, file='reslist_iscmeb_profastG_DLPFC12.rds')


# Removing unwanted variation in multi-section SRT data -------------------
load('gene_symbols_DLPFC12.rds')
str(gene_symbols)
HKsymbols
XList_hk <- pbapply::pblapply(seulist_HK, function(x) x[['RNA']]@data[names(HKsymbols),])
nvec <- sapply(XList_hk, nrow)
XList_hk <- pbapply::pblapply(XList_hk, function(x) base::t(as.matrix(x)))
Xmat_hk <- matlist2mat(XList_hk)
princ <- approxPCA(Xmat_hk, q=10)
HList <- mat2list(princ$PCs, nvec)
rm(Xmat_hk)

M <- length(XList)
p <- ncol(XList[[1]])
Tm <- matrix(0, M, 1)
correct_genesR <- ProFAST:::correct_genesR ## output the function in ProFAST because correct_genesR() is now embeded in the R package.

load('reslist_iscmeb_profastP_DLPFC12.rds')
### Poisson version:
RList_pois <- reslist_iscmeb_profastP@fitList[[1]]$Rf
tic <- proc.time()
correct_List_pois <- correct_genesR(XList, RList_pois, HList, Tm, AdjList)
toc <- proc.time()
time_all_pois <- toc[3] - tic[3] 

XList_correct_all_pois <- correct_List_pois$XList_correct
## Add row names and colnames (symbol) for each slice
for(m in 1:M){
  
  row.names( XList_correct_all_pois[[m]]) <- paste0("slice", m, "_", row.names( XList[[m]]) )
  colnames(XList_correct_all_pois[[m]]) <- unname(gene_symbols[colnames(XList[[m]])])
}
hX_pois <- matlist2mat(XList_correct_all_pois)
library(Seurat)
library(Matrix)
count <- sparseMatrix(i=1,j=1, x=1, dims=dim(t(hX_pois)))
row.names(count) <- colnames(hX_pois)
colnames(count) <- row.names(hX_pois)
seuAll <- CreateSeuratObject(counts = count)
length(unique(colnames(seuAll)))
row.names(hX_pois) <- colnames(seuAll)
seuAll[["RNA"]]@data <- t(hX_pois)
Idents(seuAll) <- factor(unlist(reslist_iscmeb_profastP@idents))
save(seuAll, file='seuAll_profastP.rds')



# Combined differential expression analysis -------------------------------
dat_degs_pois <- FindAllMarkers(seuAll)

save(dat_degs_pois, file='housekeep_noSpa_dat_degs_allgenescorrect_iscmeb_profastP.rds')
cutoff <- 0.001#0.05
dim(subset(dat_degs_pois,   p_val_adj<cutoff & avg_log2FC>0.25))

filename <- 'DLPFC12_DEGsList_iSCMEB_pois_orderedV2.xlsx'

K_use <- length(unique(dat_degs_pois$cluster))
for(k in 1:K_use){
  # k <- 8
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_pois,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  ## sort by log fold change
  dat_degs_sub3 <- dat_degs_sub3[order(dat_degs_sub3$avg_log2FC,decreasing = T),]
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),
                   append = T,row.names = F)
}




# Cell-cell interaction analysis ------------------------------------------

## use the renamed domains
str(cluster_iscmeb_pois_rename)
cluster_iscmeb_pois_rename2 <- lapply(cluster_iscmeb_pois_rename, function(x){
  x[x %in% c("Layer2/3_1", 'Layer2/3_2')] <- "Layer2/3"
  x[x=="Layer5"] <- "WM"
  return(x)
})

color_use_here <- c('#4E79A7', '#EDC948', '#F28E2B', '#B07AA1', '#76B7B2', '#E15759')
names(color_use_here) <- sort(unique(unlist(cluster_iscmeb_pois_rename2)))
# color_use_here <- color_use[sort(names(rawID_pois))]
library(CellChat)

data.input =  t(hX_pois)# normalized data matrix
meta = data.frame(labels=factor(unlist(cluster_iscmeb_pois_rename2))) # a dataframe with rownames containing cell mata data
row.names(meta) <- row.names(hX_pois)
head(meta)
# Prepare input data for CelChat analysis
dim(data.input)
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) #
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
## Write the pathway information into csv to look up
write.csv(CellChatDB$interaction, file='cellchat_pathway_info.csv')
# set the used database in the object
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

count_backup <- cellchat@net$count
count_visua <- count_backup
count_visua[count_visua<quantile(count_visua, 0.05)] <- quantile(count_visua, 0.05)
cellchat@net$count <- count_visua
het1 <- netVisual_heatmap(cellchat, signaling = NULL,measure = "count",
                          font.size = 12, font.size.title = 14, 
                          color.use = color_use_here) # color.heatmap = "Reds"
p1 = grid.grabExpr(draw(het1))
p <- plot_grid(p1)
# Save the plot as a high-resolution figure
ggsave("./cell_cell_interaction_figs/interaction_quantity_heatmap.png", p, dpi = 300, 
       width = 7, height = 8, units = "in")
weight_backup <- cellchat@net$weight
weight_visua <- weight_backup#
weight_visua[weight_visua>quantile(weight_visua, 0.97)] <- quantile(weight_visua, 0.97)
weight_visua <- weight_visua- min(weight_visua)/ (max(weight_visua)-min(weight_visua))
cellchat@net$weight <- weight_visua
het2 <- netVisual_heatmap(cellchat, signaling = NULL,measure = "weight", font.size = 12, 
                          font.size.title = 14, color.use = color_use_here)
p2 = grid.grabExpr(draw(het2))
p2 <- plot_grid(p2)
# Save the plot as a high-resolution figure
ggsave("./cell_cell_interaction_figs/interaction_strength_heatmap.png", p2, dpi = 300, 
       width = 7, height = 7, units = "in")


# PAGA trajectory analysis ------------------------------------------------

## see the sep_paga_dlpfc12.ipynb script.



