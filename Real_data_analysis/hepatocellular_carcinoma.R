# Load data and preprocessing ---------------------------------------------

source("./util_funcs.R")

library(Seurat)
library(Matrix)
library(dplyr)

setwd("./HCC4Data/")
seuList <- list()
for(iter in 1: 4){
  # iter <- 1
  cat('input HCC data', iter, '\n')
  # load and read data
  hcc <- read10XVisium(paste0("./HCC", iter))
  seuList[[iter]] <- hcc
  rm(hcc)
}

seuList <- pbapply::pblapply(seuList, filter_spot, min_feature=20) 
seuList <- pbapply::pblapply(seuList,filter_gene,min_spots=20)

HVGList <- pbapply::pblapply(seuList, function(x) {
  x <- FindVariableFeatures(x)
  x@assays$RNA@var.features}
)
save(HVGList, file='HVGList2000_HCC4.rds')

gene_intersect <- Reduce(intersect, lapply(seuList, row.names))
library(PRECAST)
data("Human_HK_genes")
idx_hk <- which(gene_intersect %in% as.character(Human_HK_genes$Gene))
housekeep_genes <- gene_intersect[idx_hk]
save(housekeep_genes, file='housekeep_genes_allList_HCC4.rds')
seulist_HK <- pbapply::pblapply(seuList, function(x) x[(housekeep_genes),])
save(seulist_HK, file='seulist_HK_HCC4.rds')




seulist_HK <- pbapply::pblapply(seulist_HK, DR.SC::FindSVGs, nfeatures = 2160)

geneList_noSpa <- lapply(seulist_HK, function(x) {
  
  order_idx <- order(x[['RNA']]@meta.features$adjusted.pval.SVGs, decreasing = T) # rank the gene with largest p-value the first
  genes <- row.names(x[['RNA']]@meta.features)[order_idx]
  # idx <- which(x[['RNA']]@meta.features$adjusted.pval.SVGs[order_idx]>0.05)
  return(genes[1:200])
}) ## return a list with ordered gene vector

HKfeatures <- selectHKFeatures(seulist_HK, geneList_noSpa, HKFeatures = 200)
HKsymbols <- housekeep_genes[housekeep_genes %in% HKfeatures]
save(HKsymbols, file='housekeep_noSpa_genes_HCC4.rds')

gene_hvg_2000 <- selectIntFeatures(seuList,HVGList) # 
save(gene_hvg_2000, file='gene_hvg_2000_HCC4.rds')

seulist_hvg <- lapply(seuList, function(x) x[gene_hvg_2000, ])
seulist_hvg <- pbapply::pblapply(seulist_hvg, filter_spot, min_feature=15) 
sapply(seulist_hvg, ncol)


posList <- list()
posList1 <- list()
AdjList <- list()
for(m in 1: length(seulist_hvg)){
  message("m = ", m)
  pos <- cbind(seulist_hvg[[m]]$row, seulist_hvg[[m]]$col)
  AdjList[[m]] <- PRECAST::getAdj_reg(pos, platform = "Visium")
  posList[[m]] <- pos + (m-1)*1e6 # to differentiate different sections
  posList1[[m]] <- pos
  
}
save(posList1, file = 'HCC4_posList1.rds')


# Spatial dimension reduction using ProFAST -------------------------------
#### The previous name of FAST is called ProFAST
library(ProFAST)
library(Seurat)
seulist_hvg <- lapply(seulist_hvg, NormalizeData)
XList_sp <- lapply(seulist_hvg, function(x) Matrix::t(x[["RNA"]]@data))
XList <- lapply(XList_sp, as.matrix)

hq <- 15
library("profmem")
options(profmem.threshold = 2000)
library(ProFAST)
pm_profastG <- profmem({
  tic <- proc.time()
  reslist_hvg_log <- ProFAST_run(XList_sp, AdjList, q=hq)
  toc <- proc.time()
})
total_mem_profastG <- sum(pm_profastG$bytes, na.rm = T)
save(pm_profastG, total_mem_profastG, file='total_mem_profastG_HCC4.rds')

time_ProFAST_hvg_log <- toc[3] - tic[3]
save(reslist_hvg_log,time_ProFAST_hvg_log, file="ProFAST_reslist_hvg_log_HCC4.rds")

## Poisson model
XList_count <- lapply(seulist_hvg, function(x) Matrix::t(x[["RNA"]]@counts))
tic <- proc.time()
pm_profastP <- profmem({
  reslist_hvg_count <- ProFAST_run(XList_count, AdjList,  q=hq, fit.model = 'poisson') 
})
total_mem_profastP <- sum(pm_profastP$bytes, na.rm = T)
save(pm_profastP, total_mem_profastP, file='total_mem_profastP_HCC4.rds')
toc <- proc.time()
time_ProFAST_hvg_count <- toc[3] - tic[3]
save(reslist_hvg_count ,time_ProFAST_hvg_count, file="ProFAST_reslist_hvg_count1_HCC4.rds")


# Downstream analysis -----------------------------------------------------
# Embedding alignment and spatial clustering ------------------------------
## output the function in ProFAST because fit.iscmeb() is now embeded in the R package but not exported.
fit.iscmeb <- ProFAST:::fit.iscmeb
### Use BIC in iSC-MEB to select the cluster number. nine clusters are identified
tic <- proc.time()
hK_count <- 4:12
tic <- proc.time()
reslist_iscmeb_pois_cho <- fit.iscmeb(
  reslist_hvg_count$hV,
  AdjList,
  K=hK_count,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 30,
  init.start=1,
  seed = 1,
  epsLogLik = 1e-05,
  coreNum = 4,c_penalty=0.5, 
  verbose = TRUE)
toc <- proc.time()
time_iscmeb_pois_cho <- toc[3] - tic[3]
## The optimal number of clusters is 9.
save(reslist_iscmeb_pois_cho,time_iscmeb_pois_cho, file='reslist_iscmeb_profastP_HCC4_2025.rds')
lapply(reslist_iscmeb_pois_cho@idents, table) 

### The results of reslist_iscmeb_pois_cho@idents is the same as the original reslist_iscmeb_pois_K9.
reslist_iscmeb_pois_K9 <- reslist_iscmeb_pois_cho 
lapply(reslist_iscmeb_pois_K9@idents, table)
clusterList_iscmeb_pois_K9 <- reslist_iscmeb_pois_K9@idents
save(reslist_iscmeb_pois_K9, file='reslist_iscmeb_pois_K9_HCC4.rds')
save(clusterList_iscmeb_pois_K9, file='clusterList_iscmeb_pois_K9_HCC4.rds')

rawID <- c(6,4,8,5,3,2,1, 9, 7)
names(rawID) <- 1:9
cluster_pois_renumber_K9 <- lapply(clusterList_iscmeb_pois_K9, replace_ID, rawID=rawID)
lapply(cluster_pois_renumber_K9, table)
save(cluster_pois_renumber_K9, file='clusterList_iscmeb_pois_reordered_K9_HCC4.rds')

# Removing unwanted variation in multi-section SRT data -------------------
XList_hk <- pbapply::pblapply(seulist_HK, function(x) x[['RNA']]@data[HKsymbols,])
nvec <- sapply(XList, nrow)
XList_hk <- pbapply::pblapply(XList_hk, function(x) base::t(as.matrix(x)))
Xmat_hk <- matlist2mat(XList_hk)
princ <- approxPCA(Xmat_hk, q=10)
HList <- mat2list(princ$PCs, nvec)
rm(Xmat_hk)
M <- length(XList)
Tm <- matrix(0, nrow=M, ncol=1)
RList_pois <- reslist_iscmeb_pois_K9@fitList[[1]]$Rf

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
Idents(seuAll_count) <- factor(unlist(cluster_pois_renumber_K9), levels=1:9)
save(seuAll_count, file='seuAll_pois_K9_HCC4.rds')



# Combined differential expression analysis -------------------------------
dat_degs_pois <- FindAllMarkers(seuAll_count)
save(dat_degs_pois, file='housekeep_noSpa_dat_degs_choseK9_iscmeb_pois1_HCC4.rds')

filename <- 'HCC4_DEGsList_iSCMEB_pois9.xlsx'
cutoff <- 0.001#
K_use <- length(unique(dat_degs_pois$cluster))
for(k in 1:K_use){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_pois,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  ## sort by log fold change
  dat_degs_sub3 <- dat_degs_sub3[order(dat_degs_sub3$avg_log2FC,decreasing = T),]
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),
                   append = T,row.names = F)
}


# Somatic mutation analysis -----------------------------------------------
# Find the HCC related SNPs
library(data.table)
dat <- fread("data_mutations.txt") ##  1-based coordinate system used by this data, thus there is no need to add 1 when we create the bed file.
dim(dat)
dat[1:5,1:10]
## Convert the data frame to a GRanges object using the GRanges() function from the GenomicRanges package. This will enable you to work with genomic intervals in R. For example, if your text file has a header row with column names "chr", "start", and "end", you could use the following command:
library(GenomicRanges)
mygranges <- GRanges(seqnames=dat$Chromosome, ranges=IRanges(start=dat$ Start_Position, end=dat$End_Position))
mybed <- data.frame(seqnames=paste0("chr",mygranges@seqnames), start=mygranges@ranges@start, end=mygranges@ranges@start, name=".", score=0, strand="+")

write.table(mybed, file="data_mutations.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


hg38_bed <- fread(file='hglft_genome_40cd7_188a60.bed')

dat_fail_conversed <- read.table("failed_SNPs.txt")

idx <- which(dat$Start_Position%in% dat_fail_conversed$V2)
dat_filter <- dat[-idx, ]
### combine the hg37 and hg38
dat_combine <- cbind(dat_filter[,1:16], hg38_bed)
save(dat_filter, file='HCC_mutation_related_variates.rds')
save(dat_combine, file='HCC_convert_basic_data.rds')


### Invetigate HCC1 to HCC4
library(data.table)
mutate_allList <- list()
for(r in 1:4){
  message("our r = ", r)
  setwd(paste0("/data/lishang/home/chaixr/project/ForLiuWei/gatk_snp_calling_bulk_HCC1-4/pileup_summary/HCC",r))
  
  files_dir <- list.files()
  idx <- grep("\\d$", files_dir, fixed=F)
  files_dir <- files_dir[idx]
  
  ## read barcode and biallelic SNPs for each spot
  mutateList <- list()
  counter <- 0
  for(file1_dir in (files_dir)){ # The unique name is chr+position
    # file1_dir <- files_dir[1]
    message('counter=',counter,", file1_dir = ", file1_dir)
    file_new <- gsub("_",".",file1_dir, fixed=T)
    mutate_dat1 <- fread(file=paste0(file1_dir, '/',"pileups.table.", file_new))
    barcode1 <- read.table(file=paste0(file1_dir, '/',"barcode.txt"), header = F)
    # read all sites then subset the variates in this spot.
    allsites_dat <- fread(file=paste0(file1_dir, '/',"pileups.table.allsites.sorted.", file_new))
    idx_mutate <- which((paste0(allsites_dat$V1,"_",allsites_dat$V2) %in% paste0(mutate_dat1$contig,"_", mutate_dat1$position)))
    flag_equal <- all(allsites_dat[idx_mutate,]$V2 == mutate_dat1$position)
    
    
    if(nrow(barcode1)==1 && ncol(barcode1)==1 && flag_equal){
      mutate_dat1$ref_allele <- allsites_dat[idx_mutate,]$V3
      mutate_dat1$alt_allele <- allsites_dat[idx_mutate,]$V4
      message("Number of good spots is ", counter)
      mutateList[[barcode1[1,1]]] <- as.data.frame(mutate_dat1)
    }else{
      stop("The barcode has problems in ", file_new)
    }
    counter <- counter + 1
  }
  
  mutate_allList[[r]] <- mutateList
}

save(mutate_allList, file='mutate_allList.rds')

library(LDlinkR)

load("HCC_convert_basic_data.rds")
idx <- which(dat_combine$dbSNP_RS!= 'novel')
dat_sub <- dat_combine[idx, ]

datLD_List <- list()
for(i in seq_along(dat_sub$dbSNP_RS)){ # 
  message("i = i", i)
  out <- try({dat <- LDproxy(dat_sub$dbSNP_RS[i], pop='ALL', r2d='r2', token='30a7ab5ff7be', genome_build='grch38')}, silent = T)
  if(class(out)=='try-error') dat <- NULL
  datLD_List[[i]] <- dat
}
names(datLD_List) <- dat_sub$dbSNP_RS
save(datLD_List, file='datLD_RSList_database_human38.rds')
name_snp <- NULL
for(i in seq_along(datLD_List)){
  
  message("i = ", i )
  if(length(nrow(datLD_List[[i]])>1)){
    if( nrow(datLD_List[[i]])>1){
      name_snp <- c(name_snp, names(datLD_List)[i])
    }
  }
  
}

datLD_snp <- datLD_List[name_snp]
save(datLD_snp, file='datLD_snp_database_human38.rds')

LDsnp <- Reduce(rbind, datLD_snp)
LDsnp_sub <- subset(LDsnp, R2>0.2)
inter_coords <- intersect(paste0(mutate_dat_allsites$V1, ":",mutate_dat_allsites$V2), LDsnp_sub$Coord)
idx <- match(inter_coords, LDsnp$Coord)

LD_snp_use <- LDsnp[idx, ]
row.names(LD_snp_use) <- inter_coords

## evaluate the freq, proportion in each cluster.
metric_names <- c("whether_variate", "alt_count", 'ref_count')
variates_hcc_arrayList <- list()
mergedatList <- list(); ss <- 1
for(r in 1:4){
  # r <- 1
  mutateList <-  mutate_allList[[r]]
  num_detected_spots <- length(mutateList)
  variates_hcc1_array <- array(0, dim=c(length(inter_coords),  num_detected_spots,  length(metric_names)))
  num_nonzero_hcc1_variation <- 0
  for(i in 1: num_detected_spots){ # 
    # i <- 1
    message("i = ", i)
    tmpdat <- mutateList[[i]]
    row.names(tmpdat) <- paste0(mutateList[[i]]$contig,":" ,mutateList[[i]]$position)
    indx1 <- which(inter_coords %in% paste0(mutateList[[i]]$contig,":" ,mutateList[[i]]$position))
    indx2 <-  which(paste0(mutateList[[i]]$contig,":" ,mutateList[[i]]$position) %in% inter_coords)
    
    
    if(length(indx1)>0){
      dat_tmp2 <- mutateList[[i]][indx2,]
      row.names(dat_tmp2) <- paste0(dat_tmp2$contig,":" , dat_tmp2$position)
      dat_tmp3 <- LD_snp_use[inter_coords[indx1],]
      dat_merge <- merge(dat_tmp3, dat_tmp2, by=0)
      for(ii in seq_along(indx1)){
        # ii <- 1
        if(dat_merge$Alleles[ii]== paste0("(",dat_merge$ref_allele[ii],"/",dat_merge$alt_allele[ii],")")){
          num_nonzero_hcc1_variation <- num_nonzero_hcc1_variation + 1
          variates_hcc1_array[indx1[ii],i, ] <- as.vector(c(1, unlist(as.vector( tmpdat[inter_coords[indx1[ii]], c("ref_count", "alt_count")]))))
        }
      }
      
      mergedatList[[ss]] <- dat_merge; ss <- ss + 1
    }
  }
  colnames(variates_hcc1_array) <- names(mutateList)
  row.names(variates_hcc1_array) <- inter_coords
  apply(variates_hcc1_array, c(1,3), sum)
  
  variates_hcc_arrayList[[r]] <- variates_hcc1_array
  
}

variates_hcc_barcode_arrayList <- list()
for(r in 1:4){
  message("r = ", r)
  variates_hcc_barcode_arrayList[[r]] <-  variates_hcc_arrayList[[r]][,barcode_hccList[[r]],]
}
save(variates_hcc_barcode_arrayList, file='variates_hcc_barcode_arrayList_LDsnps.rds')


load(file='cluster_and_barcode_hccList.rds')
load("variates_hcc_barcode_arrayList_LDsnps.rds")

cluster_level_hccAllList <- list()

for(r in 1:4){
  # r <- 1
  cluster_hcc4 <- cluster_idrsc_renumber[[r]]
  clusters_unique <- sort(unique(cluster_hcc4))
  variates_hcc4_barcode_array <- variates_hcc_barcode_arrayList[[r]]
  cluster_level_hcc4List <-list()
  for(ic in seq_along(clusters_unique)){
    # ic <- 1
    idx <- which(cluster_hcc4==clusters_unique[ic])
    if(length(idx)>1){
      cluster_level_hcc4List[[ic]] <- rbind(slice_level_altcount_hcc4 = rowSums(variates_hcc4_barcode_array[,idx,3]),
                                            slice_level_refcount_hcc4 = rowSums(variates_hcc4_barcode_array[,idx,2]),
                                            slice_level_bothcount_hcc4 = rowSums(variates_hcc4_barcode_array[,idx,3]+variates_hcc4_barcode_array[,idx,2]))
      
    }else{
      cluster_level_hcc4List[[ic]] <- c(slice_level_altcount_hcc4 = sum(variates_hcc4_barcode_array[,idx,3]),
                                        slice_level_refcount_hcc4 = sum(variates_hcc4_barcode_array[,idx,2]),
                                        slice_level_bothcount_hcc4 = sum(variates_hcc4_barcode_array[,idx,3]+variates_hcc4_barcode_array[,idx,2]))
    }
    
  }
  
  names(cluster_level_hcc4List) <- clusters_unique
  
  
  cluster_level_hccAllList[[r]] <- cluster_level_hcc4List
}

cluster_level_sumList <- list()
s <- 0
for(k in 1:9){
  s <- 0
  for(r in 1:4){
    if(!is.null(cluster_level_hccAllList[[r]][[as.character(k)]])){
      if(s == 0){
        cluster_level_sumList[[k]] <- cluster_level_hccAllList[[r]][[as.character(k)]]
      }else{
        cluster_level_sumList[[k]] <- cluster_level_sumList[[k]] + cluster_level_hccAllList[[r]][[as.character(k)]]
      }
      
      s <- s+1
    }
    
  }
}
names(cluster_level_sumList) <- 1:9
colnames(cluster_level_sumList[[ic]])
num_spots <- table(unlist(cluster_idrsc_renumber))
pList <- list()
for(k in 1:9){
  #k <- 
  ic <- 3
  ic_name <- names(cluster_level_sumList)[ic]
  tmp <- cluster_level_sumList[[ic]][1, ]/ num_spots[ic]
  tmp[order(tmp, decreasing = T)]
  pList[[k]] <- barPlot(cluster_level_sumList[[ic]][1, ]/ num_spots[ic], ylabel = 'alt_count/Num_spots') + # , cols=cols
    ggtitle(paste0("cluster ", ic)) + theme(legend.position = 'none')
  
}

