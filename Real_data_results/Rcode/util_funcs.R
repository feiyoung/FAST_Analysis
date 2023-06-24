## Created Date: 2022-09-02
## By Wei Liu
#------------ Methods' name and order
Methods_ordered <- c("SpaFactor-P","SpaFactor-G","SpatialPCA-F", "SpatialPCA-A", "DR-SC", "PCA", "multiBatchPCA",  "LIGER")
col_orderd <- c("#E04D50FF", "#F08A21FF","deepskyblue2","#4374A5FF", "#2AB673FF", "#DFE0EEFF", "#70B5B0FF", "#FCDDDEFF")
names(col_orderd) <- Methods_ordered

Methods_ordered2 <- c("ProFAST-P","ProFAST-G","SpatialPCA-F", "SpatialPCA-A", "PRECAST", "DR-SC", "scVI", "PCA", "multiBatchPCA", "NMF", "LIGER")
col_orderd2 <- c("#E04D50FF", "#F08A21FF","deepskyblue2","#4374A5FF", "#FF9896", "#2AB673FF", "#DFCDE4","#DFE0EEFF", "#70B5B0FF", "#FACB12", "#FCDDDEFF")
names(col_orderd2) <- Methods_ordered2

SpatialPCA_A_newname <- "SpatialPCA-L"

Methods_ordered3 <- c("ProFAST-P","ProFAST-G","SpatialPCA-F", "SpatialPCA-L", "PRECAST", "DR-SC", "scVI", "PCA", "multiBatchPCA", "NMF", "LIGER")
col_orderd3 <- c("#E04D50FF", "#F08A21FF","deepskyblue2","#4374A5FF", "#FF9896", "#2AB673FF", "#DFCDE4","#DFE0EEFF", "#70B5B0FF", "#FACB12", "#FCDDDEFF")
names(col_orderd3) <- Methods_ordered3

.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

unit_MB <- 8*1024^2

rename_vec <- function(namevec, oldname, newname){
  idx <- which(namevec==oldname)
  namevec[idx] <- newname
  return(namevec)
}
rename_vec(Methods_ordered2, oldname = "SpatialPCA-A", newname = "SpatialPCA-L")

#--------------Our methods:

replace_ID <- function(x, rawID){
  # x <- cluster_idrsc[[1]]
  uni_x <- unique(x)
  if(!all(uni_x %in% rawID) ) stop("Some element of x is not in rawID!")
  
  idx <- which(rawID %in% uni_x)
  sub_rawID <- rawID[idx]
  y <- rep(NA, length(x))
  
  for(k in 1:length(sub_rawID)){
    y[x==sub_rawID[k]] <- names(sub_rawID)[k]
  }
  return(y)
}


## Downstream function
get_top_pathway1 <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }
  
  return(df_sub)
}

barPlot_enrich <- function(top_dat, source='Ont', term_name="Term", nlog10P='nlog10P',
                           bar_width=0.8, base_size=20, font_family='serif', cols= ggthemes::canva_pal()(4)){
  # source='Ont'; term_name="Term"; nlog10P='nlog10P'
  require(ggplot2) # y=term_name,
  order_idx <- order(top_dat[,nlog10P])
  top_dat <- top_dat[order_idx,]
  top_dat[, term_name] <- factor(top_dat[, term_name], levels=top_dat[order_idx,term_name])
  p1 <- ggplot(data=top_dat, aes_string(x=term_name,y=nlog10P, fill=source)) +
    scale_fill_manual(values=cols)+
    geom_bar(position = "dodge", stat="identity",width =bar_width)+ coord_flip() +
    theme_classic() + theme(text=element_text(size=base_size, family=font_family)) 
  return(p1)
}

correct_genesR <- function(XList, RList, HList, Tm, AdjList, maxIter=30, epsELBO=1e-4, verbose=TRUE){
  
  library(SpaMFactor)
  M <- length(XList);  p <- ncol(XList[[1]])
  K <- ncol(RList[[1]]); q <- ncol(HList[[1]]); d <- ncol(Tm)
  Xmat <- matlist2mat(XList)
  R <- matlist2mat(RList)
  colnames(R) <- NULL
  H <- matlist2mat(HList)
  colnames(H) <- NULL
  nvec <- sapply(RList, nrow)
  TmList <- list()
  for(m in 1:M){
    TmList[[m]] <- matrix(Tm[m,], nrow=nvec[m], ncol=ncol(Tm), byrow=T)
  }
  TM <- matlist2mat(TmList)
  
  lm1 <- lm(Xmat~ R+H+TM+0)
  coef_all <- coef(lm1)
  rm(R, H, Xmat, lm1)
  alphaj_int <- coef_all[paste0("R", 1:K),]
  gammaj_int <- coef_all[paste0("H", 1:q),]
  if(d == 1){
    zetaj_int <- matrix(coef_all[paste0("TM"), ], nrow=1)
  }else{
    zetaj_int <- coef_all[paste0("TM", 1:d)]
  }
  if(sum(Tm)<1e-20){
    if(ncol(Tm)>1){
      stop("Tm must be a one-column matrix when it is full-zero!")
    }else{
      zetaj_int <- matrix(1, 1, p)
    }
  }
  
  
  sigmaj_int <- matrix(1, M, p)
  psij_int <- matrix(1,M, p)
  # maxIter <- 30; epsELBO <- 1e-4; verbose<- TRUE
  reslist <- correct_genes(XList, RList, HList, Tm, Adjlist=AdjList, sigmaj_int, psij_int, 
                              alphaj_int, gammaj_int, zetaj_int, maxIter, epsELBO, verbose) 
  # str(reslist)
  correct_list <- list(XList_correct=lapply(1:M, function(j) XList[[j]] - HList[[j]] %*% reslist$gamma),
                       gammaj=reslist$gamma, alphaj=reslist$alpha, zetaj=reslist$zeta)
  
  
  # for(m in 1:M){
  #   if(!any(sapply(XList, function(x) is.null(row.names(x))))){
  #     
  #     row.names( correct_list$XList_correct[[m]]) <- paste0("slice", m, "_", row.names(XList[[m]]) )
  #   }
  #   if(!any(sapply(XList, function(x) is.null(colnames(x))))){
  #        colnames(correct_list$XList_correct[[m]]) <- colnames(XList[[m]])
  #   }
  # }
  
  return(correct_list)
  
}

get_top_pathway <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id", "term_size","p_value")])
  }
  
  return(df_sub)
}

selectHKFeatures <- function(seulist, HKFeatureList, HKFeatures=200){
  ## This function is used for selecting common informative features
  if(length(seulist) != length(HKFeatureList)) stop("The length of suelist and HKFeatureList must be equal!")
  if(length(seulist) ==1){
    if(length(HKFeatureList[[1]]) >= HKFeatures){
      genelist <- HKFeatureList[[1]][1:HKFeatures]
    }else{
      genelist <- HKFeatureList[[1]]
      warning("The IntFeature is larger than the  number of elements in FeatureList!")
    }
    return(genelist)
  } 
  geneUnion <- unique(unlist(HKFeatureList))
  ## ensure each seuobject has the genes in geneUnion
  
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(lapply(seulist, function(x){
    assay <- DefaultAssay(x)
    geneUnion[Matrix::rowSums(x[[assay]]@counts[geneUnion,])==0]
  })))
  
  
  #geneUnion[pbapply::pbapply(x@assays$RNA@counts[geneUnion,],1, sd)==0])))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  
  # sort by number of datasets that identified this gene as the gene without spatial variation
  nsample <- length(seulist)
  numVec <- rep(0, length(gene_Var))
  rankMat <-matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], HKFeatureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(HKFeatureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
    
  }
  
  cutNum <- sort(numVec, decreasing = T)[min(HKFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(HKFeatures, length(numVec)) - length(genelist1)
  
  gene2 <- gene_Var[numVec==cutNum]
  
  
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)
  
  return(genelist)
}

#----------------------------------------functions for memory usage recording

myprofmem <- function (expr, envir = parent.frame(), substitute = TRUE, threshold = getOption("profmem.threshold", 
                                                                                 0L))
{
  require(profmem)
  profmem_stack = getFromNamespace("profmem_stack", "profmem")
  if (substitute) 
    expr <- substitute(expr)
  profmem_begin(threshold = threshold)
  error <- NULL
  value <- tryCatch({
    eval(expr, envir = envir, enclos = baseenv())
  }, error = function(ex) {
    error <<- ex
    NULL
  })
  
  pm <- myprofmem_end()
  attr(pm, "expression") <- expr
  attr(pm, "value") <- value
  attr(pm, "error") <- error
  pm
}


myprofmem_suspend <- function () 
{
  Rprofmem("")
  profmem_stack = getFromNamespace("profmem_stack", "profmem")
  profmem:::profmem_stack("suspend")
  if (profmem:::profmem_stack("depth") == 0) 
    return()
  drop <- length(sys.calls()) + 4L
  pathname <- profmem:::profmem_pathname()
  data <- myreadRprofmem(pathname, drop = drop)
  attr(data, "threshold") <- profmem_stack("threshold")
  profmem:::profmem_stack("append", data)
  invisible()
}
myreadRprofmem <- function (pathname, as = c("Rprofmem", "fixed", "raw"), 
          drop = 0L, ...){
  # as = "Rprofmem"
  profmem_stack = getFromNamespace("profmem_stack", "profmem")
  profmem:::stop_if_not(file_test("-f", pathname))
  as <- match.arg(as)
  drop <- as.integer(drop)
  profmem:::stop_if_not(length(drop) == 1, drop >= 0)
  bfr <- readLines(pathname, warn = FALSE)
  if (as == "raw") 
    return(bfr)
  pattern <- "^(new page|[0-9]+)[ ]?:(new page|[0-9]+)[ ]?:"
  while (any(grepl(pattern, bfr))) {
    bfr <- gsub(pattern, "\\1 :\n\\2 :", bfr)
    bfr <- unlist(strsplit(bfr, split = "\n", fixed = TRUE))
  }
  if (as == "fixed") 
    return(bfr)
  if (getOption("profmem.debug", FALSE)) 
    print(bfr)
  pattern <- "^([0-9]+|new page)[ ]?:(.*)"
  
  # ## I write new code:
  # mylist <- list()
  # for(i in 1: length(bfr)){
  #   # i<-1 
  #   message("i = ", i, "/", length(bfr))
  #   x <- bfr[i]
  #   bytes <- gsub(pattern, "\\1", x)
  #   what <- rep("alloc", times = length(x))
  #   idxs <- which(bytes == "new page")
  #   if (length(idxs) > 0) {
  #     what[idxs] <- "new page"
  #     bytes[idxs] <- ""
  #   }
  #   bytes <- as.numeric(bytes)
  #   trace <- gsub(pattern, "\\2", x)
  #   trace <- gsub("\" \"", "\", \"", trace, fixed = TRUE)
  #   trace <- sprintf("c(%s)", trace)
  #   # trace <- eval(parse(text = trace), enclos = baseenv())
  #   # trace <- trace[seq_len(max(0L, length(trace) - drop))]
  #   mylist[[i]] <- list(what = what, bytes = bytes, trace = trace)
  # }
  
  res <- lapply(bfr, FUN = function(x) {
    bytes <- gsub(pattern, "\\1", x)
    what <- rep("alloc", times = length(x))
    idxs <- which(bytes == "new page")
    if (length(idxs) > 0) {
      what[idxs] <- "new page"
      bytes[idxs] <- ""
    }
    bytes <- as.numeric(bytes)
    trace <- gsub(pattern, "\\2", x)
    trace <- gsub("\" \"", "\", \"", trace, fixed = TRUE)
    trace <- sprintf("c(%s)", trace)
    # trace <- eval(parse(text = trace), enclos = baseenv())
    # trace <- trace[seq_len(max(0L, length(trace) - drop))]
    list(what = what, bytes = bytes, trace = trace)
  })
  if (length(res) == 0) {
    what <- character(0L)
    bytes <- numeric(0L)
    traces <- list()
  }
  else {
    what <- unlist(lapply(res, FUN = function(x) x$what), 
                   use.names = FALSE)
    bytes <- unlist(lapply(res, FUN = function(x) x$bytes), 
                    use.names = FALSE)
    traces <- lapply(res, FUN = function(x) x$trace)
  }
  res <- data.frame(what = what, bytes = bytes, stringsAsFactors = FALSE)
  res$trace <- traces
  bfr <- bytes <- traces <- NULL
  class(res) <- c("Rprofmem", class(res))
  profmem:::stop_if_not(all(c("what", "bytes", "trace") %in% 
                    names(res)))
  res
}

myprofmem_end <- function () 
{
profmem_stack = getFromNamespace("profmem_stack", "profmem")
  myprofmem_suspend()
  depth <- profmem:::profmem_stack("depth")
  if (depth == 0) {
    stop("Did you forget to call profmem_begin()?")
  }
  data <-  profmem:::profmem_stack("pop")
  profmem_resume()
  data
}

#---------------------------------------Memory usage recording


RCTD_structure <- function(sc_obj, clust_vr) {
  
  sc_obj[["Name"]] = sc_obj@meta.data[, clust_vr]
  
  # Cell type dictionary between cluster and cell-type
  ct <- unique(sc_obj@meta.data[, clust_vr])
  df_ct <- data.frame("Cluster" = 1:length(ct),
                      "Name" = ct)
  metadata <- sc_obj@meta.data %>%
    # Rownames to columns must be before left join since after it the rownames are erased
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(df_ct, by = c("Name" = "Name")) %>%
    # Change names to 鈥渂arcode�? 鈥渃luster�? 鈥渘UMI�?   
    mutate(
    cluster = Cluster,
  nUMI = nCount_RNA
  ) %>%
  dplyr::select(barcode, cluster, nUMI)

expr_mtrx <- sc_obj@assays$RNA@counts

return(list("meta_data" = metadata,
            "cell_type_dict" = df_ct,
            "dge" = expr_mtrx))
}


RCTD_run <- function(sc_data,cell_type_columns,st_count,st_coord,output_dir,output_name){   
  #sc_data: seurat object; Rownames should be genes and colnames represent cells/barcodes names.
  #st_data: Rownames should be genes and colnames represent barcodes/pixel names.
  library(tidyr)
  library(dplyr)
  library(spacexr)
  library(Seurat)
   library(Matrix)
  
  print('Prepare SC data')
  sc_ls <- RCTD_structure(sc_obj = sc_data,clust_vr = cell_type_columns)   #use which column as cell type
  meta_data <- sc_ls[[1]]
  sc_counts <- sc_ls[[3]]
  cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- as.numeric(meta_data$nUMI); names(nUMI) <- meta_data$barcode # create nUMI named list
 
  reference <- Reference(sc_counts, cell_types, nUMI)
 
  print('Prepare ST data')
  puck <- SpatialRNA(st_coord, st_count)  # rownames are barcodes/pixel names
 
  print('Run RCTD')
  # Creating and running RCTD
  myRCTD <- create.RCTD(spatialRNA = puck,
                        reference = reference,
                        max_cores = 1, UMI_min = 0, UMI_max = 1e20,counts_MIN =0,
                        CELL_MIN = 18)
  myRCTD <- run.RCTD(RCTD = myRCTD,
                     doublet_mode = 'doublet')
  deconv_results <- myRCTD@results
  
  print('Save result')
  dir.create(file.path(output_dir), showWarnings = FALSE)
  saveRDS(object = myRCTD,paste(output_dir,"/deconv_results_",output_name,".RDS",sep=""))
  
  #normalize the cell type proportions to sum to 1.
  rctd_deconv <- sweep(deconv_results$weights, 1, rowSums(deconv_results$weights), '/') 
  rctd_deconv = as.matrix(rctd_deconv)
  rctd_deconv[rctd_deconv<1e-3] = 0
  colnames(rctd_deconv) = sc_ls[[2]][,'Name']
  
  ### Save results
  write.csv(rctd_deconv, file = paste(output_dir,"/output_weights_",output_name,".csv",sep=""))
}



align_colors <- function(ClusterMat,base_cluster, cols_cluster){
  # j <- 1
  if(length(cols_cluster) < max(apply(ClusterMat, 2, max))) stop("Number of cols is not enough!")
  givecolors <- function(mapID, base_cols_cluster, other_cols){
    m <- ncol(mapID)
    color_ID <- unname(base_cols_cluster)
    row_clusterID <- as.numeric(names(mapID))
    row_clusterCols <- rep("", m)
    col_used_count <- rep(0, length(base_cols_cluster))
    tmp_base_cluster <- NULL
    for(i in 1: m){
      # i <- 1
      
      if((mapID[1,i] %in% color_ID) &&  (max(mapID[2, mapID[1,] == mapID[1,i]]) == mapID[2,i])){
        col_chose <- names(base_cols_cluster)[mapID[1, i]]
        row_clusterCols[i] <- col_chose
        tmp_base_cluster <- c(tmp_base_cluster, mapID[1, i])
        color_ID <- setdiff(color_ID, mapID[1, i])
        
      }else{
        row_clusterCols[i] <- other_cols[1]
        other_cols <- other_cols[-1]
      }
      
      col_used_count[mapID[1, i]] <- col_used_count[mapID[1, i]] + 1
      
      
    }
    if(any(is.na(row_clusterCols))) row_clusterCols[is.na(row_clusterCols)] <- names(base_cols_cluster)[color_ID]
    return(row_clusterCols)
  }
  base1 <- sort(unique(base_cluster))
  n_base <- length(base1)
  base_cols_cluster <- base1
  names(base_cols_cluster) <- cols_cluster[1:n_base]
  other_cols <- cols_cluster[-(1:n_base)]
  colorList <- list()
  for(j in 1:ncol(ClusterMat)){
    # j <- 1
    stab <- table(ClusterMat[,j], base_cluster)
    mapID <- rbind(apply(stab, 1, which.max), apply(stab, 1, max))
    colorList[[j]] <- givecolors(mapID, base_cols_cluster, other_cols)
  }
  return(colorList)
}




# Functions for simulation-----------------

# To make CellAssign works well, we redefine a simulation function: genSimDataFromRealData <- function()
genSimDataFromRealData <- function(seu,  sim_seed=1, ngenes = 2000,  batch_facLoc=0.1, batch_facScale = 0.1, 
                                   de_facLoc = 0, de_facScale = 1, de_prop = 0.5, 
                                   findMarkers=FALSE, do.plotBatch=FALSE){
  
  library(SingleCellExperiment)
  require(Seurat)
  library(splatter)
  pos = seu@meta.data[,c("row", "col")]
  y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  
  cnts = as.matrix(seu[["RNA"]]@counts)
  init_params <- splatEstimate(cnts)
  
  num_mk_per_ct = 5
  
  C = 8
  I = NULL
  N = 7000
  L = 3
  
  de_propvec = rep(de_prop,8)
  debug = FALSE
  
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  group_prob <- as.vector(table(y)/length(y))
  
  
  params <- setParams(
    init_params, 
    batchCells = rep(N, L), # 3N here represents a large number such that
    # we have sufficient cells of each type to be 
    # allocated to the spatial transcriptomics data
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale, ## 
    nGenes = ngenes,
    group.prob = group_prob,
    out.prob = 0,
    de.prob = de_propvec,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)
  
  sim_groups <- splatSimulate(
    params = params, 
    method = "groups", 
    verbose = FALSE)
  
  library(Seurat)
  library(scuttle)
  sim_groups <- logNormCounts(sim_groups)
  
  if(do.plotBatch){
    library(scater)
    sim_groups <- runPCA(sim_groups)
    plt <- plotPCA(sim_groups, colour_by = "Batch")
    require(ggplot2)
    if(!dir.exists("Figs")){
      dir.create("Figs")
    }
    filename <- paste0("Batch_",as.character(Sys.time()))
    ggsave(file=paste0("./Figs/",filename,".png"), plot = plt,
           width =7, height =5.5, units = "in", dpi = 200)
  }
  
  if(findMarkers){
    seu = as.Seurat(sim_groups)
    seu = SCTransform(seu, assay = "originalexp")
    Idents(seu) = seu@meta.data$Group
    all.markers = FindAllMarkers(seu, assay = "SCT", logfc.threshold = 0.1)
    
    library(dplyr)
    all.markers %>%
      group_by(cluster) %>%
      top_n(n = num_mk_per_ct, wt = avg_log2FC) -> top5
    
  }else{
    all.markers <- top5 <- NULL
  }
  
  
  meta_data_all <- colData(sim_groups)
  # y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  # pos <- cbind(seu$row, seu$col)
  Groupy = paste0("Group", y)
  K <- length(unique(y))
  nbatch <- 3
  
  idxList_sim <- list()
  idxList_pos <- list(1:length(y), which(pos[,1]< quantile(pos[,1], 0.9)),
                      which(pos[,2]< quantile(pos[,2], 0.9)))
  ### Ensure the locations of each spot
  posList <- lapply(idxList_pos, function(idx) pos[idx, ])
  yList <- lapply(idxList_pos, function(idx) y[idx])
  for(r in 1:nbatch){
    message("r = ", r)
    y_tmp <- yList[[r]]
    num_each_celltype = table(y_tmp)
    K_tmp <- length(unique(y_tmp))
    idx1 = rep(0, length(y_tmp))
    for (i in 1:K_tmp){
      idx1[y_tmp==i] = which(meta_data_all[,3] == paste0("Group",i) & meta_data_all[,2]==paste0("Batch", r))[1:num_each_celltype[i]] 
    }
    idxList_sim[[r]] <- idx1 ## align with posList
    
  }
  
  sceList_sim <- lapply(idxList_sim, function(idx) sim_groups[,idx])
  ## Add spatial coordinates
  sceList_sim <- lapply(1: nbatch, function(r){
    sce <- sceList_sim[[r]]
    colData(sce)$row <- posList[[r]][,1]
    colData(sce)$col <- posList[[r]][,2]
    return(sce)
  })
  
  sce2seurat <- function(sce){
    # sce <- sceList_sim[[1]]
    require(SingleCellExperiment)
    count <- counts(sce)
    meta_data <- as.data.frame(colData(sce))
    require(Seurat)
    seu <- CreateSeuratObject(counts=count, meta.data = meta_data)
    
    return(seu)
  }
  
  seuList_use <- pbapply::pblapply(sceList_sim, sce2seurat)
  
  out_seu = list(top5 = top5, seuList = seuList_use, all.markers = all.markers)
  
  #out = list(top5 = top5, sim_groups = sim_groups, all.markers = all.markers)
  return(out_seu)
  
}



## Generate batch data based a real data
generate_seuList_marker <- function(seu, sim_seed=1, ngenes = 2000, batch_facLoc=0.1, batch_facScale = 0.1,
                             de_facLoc=0.1, de_facScale=0.4, de_prop = 0.5, findMarker=TRUE){
  
  library(SingleCellExperiment)
  library(splatter)
  ## read position and annotation label from real data dlpfc
  pos = seu@meta.data[,c("row", "col")]
  
  y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  cnts = as.matrix(seu[["RNA"]]@counts)
  init_params <- splatEstimate(cnts)
  
  
  #batch_facLoc = 0.1 ## Batch effects in location
  #batch_facScale = 0.1
  C = length(unique(y)) # the number of spatial Domains
  I = NULL
  N = 5000
  L = 3 ## 3 batches
  num_mk_per_ct = 5
  debug = FALSE
  
  # 1.simulate count data
  
  group_prob <- as.vector(table(y)/length(y))
  
  
  
  
  
  
  ## Generate data with batch effects
  params <- setParams(
    init_params, 
    batchCells = rep(N, L), # 3N here represents a large number such that
    # we have sufficient cells of each type to be 
    # allocated to the spatial transcriptomics data
    batch.rmEffect = FALSE,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale,
    nGenes = ngenes,
    group.prob = group_prob,
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)
  
  
  sim_groups <- splatSimulate(
    params = params, 
    method = "groups", 
    verbose = FALSE)
  
  ### Find marker genes
  if(findMarker){
    ## use one data batch to find markers
    library(Seurat)
    library(scuttle)
    sim_group1 <- sim_groups[,colData(sim_groups)$Batch=='Batch1']
    meta_data_all <- colData(sim_group1)
    require(SingleCellExperiment)
    count <- counts(sim_group1)
    seu1 = CreateSeuratObject(counts=as.sparse(count), meta.data=as.data.frame(meta_data_all))
    seu1 = NormalizeData(seu1)
    Idents(seu1) = seu1@meta.data$Group
    all.markers = FindAllMarkers(seu1, assay = "RNA")
    library(dplyr)
    all.markers %>%
      group_by(cluster) %>%
      top_n(n = num_mk_per_ct, wt = avg_log2FC) -> top5
  }else{
    all.markers <- top5 <- NULL
  }
  
  
  
  ## reorder the sample to  be the same as y
  meta_data_all <- colData(sim_groups)
  y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  Groupy = paste0("Group", y)
  K <- length(unique(y))
  nbatch <- 3
  pos <- cbind(seu$row, seu$col)
  idxList_sim <- list()
  idxList_pos <- list(1:length(y), which(pos[,1]< quantile(pos[,1], 0.9)),
                      which(pos[,2]< quantile(pos[,2], 0.9)))
  ### Ensure the locations of each spot
  posList <- lapply(idxList_pos, function(idx) pos[idx, ])
  yList <- lapply(idxList_pos, function(idx) y[idx])
  for(r in 1:nbatch){
    message("r = ", r)
    y_tmp <- yList[[r]]
    num_each_celltype = table(y_tmp)
    K_tmp <- length(unique(y_tmp))
    idx1 = rep(0, length(y_tmp))
    for (i in 1:K_tmp){
      idx1[y_tmp==i] = which(meta_data_all[,3] == paste0("Group",i) & meta_data_all[,2]==paste0("Batch", r))[1:num_each_celltype[i]] 
    }
    idxList_sim[[r]] <- idx1 ## align with posList
    
  }
  
  sceList_sim <- lapply(idxList_sim, function(idx) sim_groups[,idx])
  ## Add spatial coordinates
  sceList_sim <- lapply(1: nbatch, function(r){
    sce <- sceList_sim[[r]]
    colData(sce)$row <- posList[[r]][,1]
    colData(sce)$col <- posList[[r]][,2]
    return(sce)
  })
  
  sce2seurat <- function(sce){
    # sce <- sceList_sim[[1]]
    require(SingleCellExperiment)
    count <- counts(sce)
    meta_data <- as.data.frame(colData(sce))
    require(Seurat)
    seu <- CreateSeuratObject(counts=count, meta.data = meta_data)
    
    return(seu)
  }
  
  seuList_use <- pbapply::pblapply(sceList_sim, sce2seurat)
  
  out_seu = list(top5 = top5, seuList = seuList_use, all.markers = all.markers)
  return(out_seu)
}



merge_genes_seu2 <- function(seu_marker, seu_nonmarker, assay="RNA"){
  require(Seurat)
  count1 <- seu_marker[[assay]]@counts
  row.names(count1) <- paste0(row.names(count1), ".marker")
  
  count2 <- seu_nonmarker[[assay]]@counts
  row.names(count2) <- paste0(row.names(count2), ".nonmarker")
  count_combine <- rbind(count1, count2)
  
  seu_combine <- CreateSeuratObject(counts=count_combine, meta.data = seu_marker@meta.data)
  return(seu_combine)
}
calMarkerIndMatrix <- function(markerVec, add_unknown=TRUE){
  
  # markerVec <- marker_genes
  if(is.null(markerVec)){
    stop("Check markerVec! markerVec must be a named vector!")
  }
  celltypes <- unique(names(markerVec))
  K <- length(celltypes)
  marker_unique <- unname(markerVec[!duplicated(markerVec)])
  d <- length(marker_unique)
  D <- length(markerVec)
  A <- matrix(0, d, K)
  row.names(A) <- marker_unique
  colnames(A) <- celltypes
  for(j in 1:d){
    
    sj <- markerVec[markerVec== marker_unique[j]]
    A[j, names(sj)] <- 1
    
  }
  if(add_unknown){
    A <- cbind(A, "Unknown"=0)
  }
  return(A)
}


generate_count <- function(seu, sim_seed=1,J = 2000, batch_facLoc=0.1, batch_facScale = 0.1){
  
  library(SingleCellExperiment)
  library(splatter)
  ## read position and annotation label from real data dlpfc
  
  pos = seu@meta.data[,c("row", "col")]
  y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  
  cnts = as.matrix(seu[["RNA"]]@counts)
  init_params <- splatEstimate(cnts)
  
  
  #batch_facLoc = 0.1 ## Batch effects in location
  #batch_facScale = 0.1
  C = length(unique(y)) # the number of spatial Domains
  I = NULL
  N = 7000
  L = 3 ## 3 batches
  
  debug = FALSE
  
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  group_prob <- as.vector(table(y)/length(y))
  
  
  params <- setParams(
    init_params, 
    batchCells = rep(N, L), # 3N here represents a large number such that
    # we have sufficient cells of each type to be 
    # allocated to the spatial transcriptomics data
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale,
    nGenes = J,
    group.prob = group_prob,
    seed = sim_seed)
  
  sim_groups <- splatSimulate(
    params = params, 
    method = "groups", 
    verbose = FALSE)
  return(sim_groups)
  
}
generate_seuList1 <- function(seu, sim_seed=1, ngenes = 2000, batch_facLoc=0.1, batch_facScale = 0.1){
  
  sim_groups <- generate_count(seu, sim_seed = sim_seed, J= ngenes, 
                               batch_facLoc= batch_facLoc, batch_facScale = batch_facScale)
  ## reorder the sample to  be the same as y
  meta_data_all <- colData(sim_groups)
  y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  Groupy = paste0("Group", y)
  nbatch <- 3
  pos <- cbind(seu$row, seu$col)
  idxList_sim <- list()
  idxList_pos <- list(1:length(y), which(pos[,1]< quantile(pos[,1], 0.9)), which(pos[,2]< quantile(pos[,2], 0.9)))
  posList <- lapply(idxList_pos, function(idx) pos[idx, ])
  yList <- lapply(idxList_pos, function(idx) y[idx])
  for(r in 1:nbatch){
    message("r = ", r)
    y_tmp <- yList[[r]]
    num_each_celltype = table(y_tmp)
    idx1 = rep(0, length(y_tmp))
    for (i in 1:7){
      idx1[y_tmp==i] = which(meta_data_all[,3] == paste0("Group",i) & meta_data_all[,2]==paste0("Batch", r))[1:num_each_celltype[i]] 
    }
    idxList_sim[[r]] <- idx1 ## align with posList
    
  }
  
  sceList_sim <- lapply(idxList_sim, function(idx) sim_groups[,idx])
  ## Add spatial coordinates
  sceList_sim <- lapply(1: nbatch, function(r){
    sce <- sceList_sim[[r]]
    colData(sce)$row <- posList[[r]][,1]
    colData(sce)$col <- posList[[r]][,2]
    return(sce)
  })
  
  sce2seurat <- function(sce){
    # sce <- sceList_sim[[1]]
    require(SingleCellExperiment)
    count <- counts(sce)
    meta_data <- as.data.frame(colData(sce))
    require(Seurat)
    seu <- CreateSeuratObject(counts=count, meta.data = meta_data)
    
    return(seu)
  }
  
  seuList_use <- pbapply::pblapply(sceList_sim, sce2seurat)
  
  return(seuList_use)
}




### generate latent feature V from CAR model for each sample
cor.mat <- function(p, rho, type='toeplitz'){
  
  mat <- diag(p)
  if(type=='toeplitz'){
    for(i in 2:p){
      for(j in 1:i){
        mat[i,j] <- mat[j,i] <- rho^(abs(i-j))
      }
    }
  }
  if(type=='identity'){
    mat[mat==0] <- rho
  }
  return(mat)
}
cov.mat <- function(sdvec, rho, type='toeplitz'){
  p <- length(sdvec)
  cormat <- cor.mat(p, rho, type)
  covmat <- matrix(0, p, p)
  for(i in 1:p){
    for(j in 1:p){
      covmat[i,j] <- sdvec[i]*sdvec[j]*cormat[i,j]
    }
  }
  return(covmat)
}
normalizeMatrix <- function(adj){
  require(Matrix)
  s <- colSums(adj)
  n <- nrow(adj)
  adj_norm <- adj
  for(i in 1:n){
    adj_norm[i, i:n] <- adj[i, i:n]/ (2*s[i])
    adj_norm[i:n, i] <-  adj_norm[i, i:n]
  }
  return(adj_norm)
}


# used colors  ------------------------------------------------------------

get_col_clusters <- function(n){
  library(ggthemes)
  palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
  names(palettes)
  pal1 <- tableau_color_pal("Classic 20")
  max_n <- attr(pal1, "max_n")
  pal2 <- tableau_color_pal("Classic Blue-Red 12")
  max_n2 <- attr(pal2, "max_n")
  col_cluster <- c(pal1(max_n), pal2(max_n2))
  if(n > length(col_cluster)) stop(paste0("n must be less than or equal to ", length(col_cluster)))
  col_clusters <- col_cluster[1:n]
  names(col_clusters) <- 1:n
  return(col_clusters)
}



get_col_samples <- function(n){
  cols <- PRECAST:::gg_color_hue(n)
  names(cols) <- 1:n
  return(cols);
}
  

# Visualization: Plot functions ----------------------------------------------------------

drawFigs <- function(pList, layout.dim = NULL, common.legend=FALSE,legend.position='right',  ...){
  if(!is.list(pList)) stop('drawFigs: pList must be a list!')

  if(is.null(layout.dim) && length(pList)>1){
    layout.dim <- c(2, round(length(pList)/2) )
  }
  if(is.null(layout.dim) && length(pList) == 1){
    layout.dim <- c(1,1)
  }
  ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                             nrow = layout.dim[1], common.legend = common.legend,
                             legend = legend.position, ...)

}
chooseColors <- function(palettes_name= c("Nature 10", "CNS 10","Light 13", "Classic 20", "Blink 23", "Hue n"), n_colors = 7,
                         alpha=1, plot_colors=FALSE){
  
  palettes_name <- match.arg(palettes_name)
  colors <- if(palettes_name == "Classic 20"){
    if(n_colors>20) stop("n_colors can not exceed 20 for Classic 20 schame!")
    require(ggthemes)
    # palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
    pal1 <- tableau_color_pal(palettes_name)
    pal1(n_colors)
  }else if(palettes_name == "Nature 10"){
    if(n_colors>10) stop("n_colors can not exceed 10 for Nature 10 schame!")
    cols <- c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",
              "#70B5B0", "#DFE0EE" ,"#DFCDE4", "#FACB12", "#f9decf")
    idx <- round(seq(1,10, length=n_colors))
    cols[idx]
  }else if(palettes_name == "CNS 10"){
    if(n_colors>10) stop("n_colors can not exceed 10 for Nature 10 schame!")
    cols <-CNS10 <- c("#9C0141", "#D53F4D", "#F47646", "#FBB969", "#FEE99A", 
                      "#F6FAAD", "#CAE89C", "#87D0A5", "#469DB4", "#4E62AA")
    idx <- round(seq(1,10, length=n_colors))
    cols[idx]
  }else if(palettes_name == "Blink 23"){
    if(n_colors>23) stop("n_colors can not exceed 23 for Blink 23 schame!")
    cols <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#FE0092", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
              "#D35400", "#00eefd", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
              "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#007CC8")
    idx <- round(seq(1,23, length=n_colors))
    cols[idx]
  }else if(palettes_name == "Light 13"){
    if(n_colors>13) stop("n_colors can not exceed 13 for Light 13 schame!")
    cols <-c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
              "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
              "#FF9896","#91D1C2", "#C7E9C0" ,
              "#6B6ECF", "#7B4173" )
    idx <- round(seq(1,13, length=n_colors))
    cols[idx]
  }else if(palettes_name == "Hue n"){
    .gg_color_hue(n_colors)
  }else{
    stop(paste0("chooseColors: check palettes_name! Unsupported palettes_name: ", palettes_name))
  }
  require(colorspace)
  colors_new = adjust_transparency(colors,   alpha = alpha)
  if(plot_colors){
    barplot(rep(1, length(colors_new)), axes = FALSE, space = 0, col = colors_new)
  }
  
  return(colors_new)
}

# chooseColors("CNS 10", 5, plot_colors = T)

plot_scatter <- function (
  embed_use, meta_data, label_name, xy_names=c('tSNE1', 'tSNE2'), 
  title_name = NULL, no_guides = FALSE, 
  cols = NULL, filename=NULL, y_reverse=FALSE, x_reverse=FALSE,
  point_size = 0.5, point_alpha=1, point_shape=15,
  base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
  legend_pos='right', legend_dir='vertical') {
  require(dplyr)
  require(ggthemes)
  require(ggrepel)
  require(data.table)
  require(ggplot2)
  plt_df <- embed_use %>% data.frame() %>% cbind(meta_data) %>% 
    dplyr::sample_frac(1L)
  plt_df$given_name <- plt_df[[label_name]]
  
  if(!is.null(names(cols))){
    
    cluster <- as.vector(plt_df$given_name)
    uni_cluster <- levels(plt_df$given_name)
    
    
    
    cols  <- cols[uni_cluster]
    cols <- unname(cols)
  }else{
    stop("cols must be named vector!")
  }
  
  
  
  plt <- plt_df %>% ggplot(aes_string(colnames(plt_df)[1],colnames(plt_df)[2], col = label_name, 
                                      fill = label_name)) + #  + theme_tufte(base_size = base_size, ticks= show_ticks)
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.direction = legend_dir, legend.position = legend_pos,
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col))+
    guides(color = guide_legend(override.aes = list(stroke = 1, 
                                                    alpha = 1, shape = 16, size = 4)), alpha = "none") + 
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(x = xy_names[1], 
                                                         y = xy_names[2])
  if (do_points) 
    plt <- plt + geom_point( size = point_size, alpha=point_alpha, shape=point_shape)
  if(!is.null(title_name)) plt <- plt +ggplot2::ggtitle(title_name)
  if (do_density) 
    plt <- plt + geom_density_2d()
  if (no_guides) 
    plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
  if(y_reverse)
    plt <- plt + scale_y_reverse()
  if(x_reverse)
    plt <- plt + scale_x_reverse()
  if(!is.null(filename)){
    if(!dir.exists("Figs")){
      dir.create("Figs")
    }
    ggsave(file=paste0("./Figs/",filename,".png"), plot = plt,
           width =7, height =5.5, units = "in", dpi = 200)
  }
  return(plt)
}
# plot_scatter(embed_use, meta_data, label_name='cluster', cols=cols)
data_save <- function(obj, filename='filename', dir_name='data'){
  if(!dir.exists(dir_name))
    dir.create(dir_name)
  saveRDS(obj, file=paste0("./", dir_name, "/", filename))
}





write_fig <- function(plt, filename=NULL, y_reverse=FALSE, x_reverse=FALSE,width =7, height =5.5, dpi = 200){
  
  require(ggplot2)
  if(y_reverse)
    plt <- plt + scale_y_reverse()
  if(x_reverse)
    plt <- plt + scale_x_reverse()
  if(!is.null(filename)){
    if(!dir.exists("Figs")){
      dir.create("Figs")
    }
    ggsave(file=paste0("./Figs/",filename,".png"), plot = plt,
           width =width, height =height, units = "in", dpi = dpi)
  }
}
# dat <- PRECAST:::gendata_seulist()
# posList <- lapply(dat, function(x) cbind(x$row, x$col))
# clusterList <- lapply(dat, function(x) x$true_cluster)
# clusterList[[1]][clusterList[[1]]==1 | clusterList[[1]]==4 | clusterList[[1]]==5] <- 2
# clusterList[[2]][clusterList[[2]]==3] <- 1
# posList[[3]] <- posList[[1]]
# clusterList[[3]] <- dat[[1]]$true_cluster
combinePlot <- function(posList, clusterList, item='Cluster', point_size=2,text_size=12, 
                        cols=NULL,font_family='', border_col="gray10", y_reverse=FALSE, x_reverse=FALSE,
                        fill_col='white', ncol=2, combine = TRUE, title_name="Sample"){
  
  # item='Cluster'; point_size=2;text_size=12;
  # cols=NULL;font_family=''; border_col="gray10";
  # fill_col='white'; ncol=2; combine = TRUE; title_name="Sample"
  ## Check arguments input
  plot_scatter <- function (
    embed_use, meta_data, label_name, xy_names=c('tSNE1', 'tSNE2'), 
    title_name = NULL, no_guides = FALSE, 
    cols = NULL, filename=NULL, y_reverse=FALSE, x_reverse=FALSE,
    point_size = 0.5, point_alpha=1, point_shape=15,
    base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
    legend_pos='right', legend_dir='vertical') {
    require(dplyr)
    require(ggthemes)
    require(ggrepel)
    require(data.table)
    require(ggplot2)
    plt_df <- embed_use %>% data.frame() %>% cbind(meta_data) %>% 
      dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    
    if(!is.null(names(cols))){
      
      cluster <- as.vector(plt_df$given_name)
      uni_cluster <- levels(factor(cluster))
      idx <- which(names(cols) %in% uni_cluster) # sort in the order of colnames
      cols  <- cols[idx]
      
      cols <- unname(cols[uni_cluster])
    }else{
      stop("cols must be named vector!")
    }
    
    
    
    plt <- plt_df %>% ggplot(aes_string(colnames(plt_df)[1],colnames(plt_df)[2], col = label_name, 
                                        fill = label_name)) + #  + theme_tufte(base_size = base_size, ticks= show_ticks)
      theme(axis.text.x=element_text(size=base_size, color=1),
            axis.text.y=element_text(size=base_size, color=1),
            axis.title.x = element_text(size=base_size+2, color='black'),
            axis.title.y = element_text(size=base_size+2, color='black'),
            strip.text =  element_text(size=base_size, color='black'),
            strip.background = element_rect(
              linetype = 'solid', color='gray3'
            ),
            legend.direction = legend_dir, legend.position = legend_pos,
            legend.text=element_text(size=base_size+1),
            legend.title=element_text(size=base_size+2),
            panel.background= element_rect(fill = 'white', color=border_col))+
      guides(color = guide_legend(override.aes = list(stroke = 1, 
                                                      alpha = 1, shape = 16, size = 4)), alpha = "none") + 
      scale_color_manual(values = cols) + scale_fill_manual(values = cols) + 
      theme(plot.title = element_text(hjust = 0.5)) + labs(x = xy_names[1], 
                                                           y = xy_names[2])
    if (do_points) 
      plt <- plt + geom_point( size = point_size, alpha=point_alpha, shape=point_shape)
    if(!is.null(title_name)) plt <- plt +ggplot2::ggtitle(title_name)
    if (do_density) 
      plt <- plt + geom_density_2d()
    if (no_guides) 
      plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
    if(y_reverse)
      plt <- plt + scale_y_reverse()
    if(x_reverse)
      plt <- plt + scale_x_reverse()
    if(!is.null(filename)){
      if(!dir.exists("Figs")){
        dir.create("Figs")
      }
      ggsave(file=paste0("./Figs/",filename,".png"), plot = plt,
             width =7, height =5.5, units = "in", dpi = 200)
    }
    return(plt)
  }
  # plot_scatter(embed_use, meta_data, label_name='cluste
  sampleID <- PRECAST:::get_sampleID(posList)
  posall <- PRECAST:::matlist2mat(posList)
  meta.data <- data.frame(row=posall[,1], col=posall[,2],
                          Sample=factor(sampleID), Cluster=factor(unlist(clusterList)))
  
  row.names(meta.data) <- paste0("spot", 1:nrow(meta.data))
  
  if(item %in% colnames(meta.data)){
    if(!is.factor(meta.data[, item])) meta.data[, item] <- factor(meta.data[, item])
    
    meta.data$tmp_item_id <- as.numeric(meta.data[, item])
  }
  if(is.null(cols)){
    # to determine the number of colors
    
    nn <- length(unique(meta.data[,item]))
    cols <- PRECAST:::gg_color_hue(nn)
    names(cols) <- 1:nn
  }
  
  
  if(!is.vector(cols) && item != "RGB_UMAP" && item!="RBG_tSNE")
    stop("Check argument: cols! it must be a vector object.")
  
  
  ###Finish  Check  of arguments 
  
  batch_vec <- sort(unique(sampleID))
  pList <- list()
  k <- 1
  item_null_flag <- FALSE
  for(batchi in batch_vec){
    # batchi <- (batch_vec[1])
    meta_data_sub <- subset(meta.data, Sample==batchi)
    
    
    embed_use <- meta_data_sub[,c("row", "col")]
    if(item %in% colnames(meta_data_sub)){
      
      sort_id <- sort(unique(meta_data_sub[, 'tmp_item_id']))
      p1 <- plot_scatter(embed_use, meta_data_sub[item], label_name=item, 
                         point_size=point_size, cols =cols) # [sort_id]
    }
    # simutool::colorbar_adj_transparent(cols[sort_id],1)
    # meta_data <- meta_data_sub[item]; label_name=item;point_size=point_size;cols =cols[sort_id]
    # height <- rep(1, length(sort_id))
    # barplot(height, col = cols[sort_id])
    
    p1 <- p1 + PRECAST:::mytheme_graybox(base_size = text_size, base_family = font_family, bg_fill = fill_col,
                                         border_color = border_col) 
    if(!is.null(title_name)){
      p1 <- p1 + ggtitle(label=paste0(title_name, batchi))
    }
    if(y_reverse)
      p1  <- p1  + scale_y_reverse()
    if(x_reverse)
      p1  <- p1  + scale_x_reverse()
    
    pList[[k]] <- p1
    k <- k + 1
    if(item_null_flag){
      item <- NULL
    }
  }
  if(combine){
    
    pList <- patchwork::wrap_plots(pList, ncol=ncol)
  }
  return(pList)
}

rotate_90_clockwise <- function(pl){
  require(ggplot2)
  pl + coord_flip() + scale_x_reverse()
}



# Metric functions---------------------------------------------------------
## ## Batch effect removal
## average ilsi scores (Harmony, 2019, Nature Method)

ilsi_avg_scores_allspots <- function(embeddings,  category){
  require(scPOP)
  metadata <- data.frame(category=category)
  lisi_scores_pro <- lisi(embeddings, meta_data = metadata, 'category')
  lisi_scores_pro$category # return a vector
}






F1_score_silho_allspots <- function(embeddings, celltype, sampleID){
  #require(scPOP)
  metadata <- data.frame(celltype=celltype, sampleID = sampleID)
  dd <- dist(embeddings)
  sh_scores_pro <- sapply(names(metadata), function(x) {
    cluster::silhouette(as.numeric(as.factor(metadata[[x]])), 
                        dd)[,3]
  })
  
  sh_ct <- (1+sh_scores_pro[,1])/2 # larger is better
  sh_si <- (1+sh_scores_pro[,2])/2 # smaller is better
  f1_score <- (2* (1-sh_si)*sh_ct) / ((1-sh_si) + sh_ct)
  return(f1_score)
}


## average silhouette width (ASW, review, 2021, Nature Method)
silhouette_avg_width <- function(embeddings,  category){
  require(scPOP)
  metadata <- data.frame(category=category)
  # tSNE <- scater::calculateTSNE(t(embeddings))
  sh_scores_pro <- silhouette_width(embeddings, meta.data = metadata, 'category')
  sh_scores_pro
}

lisi_avg_scores <- function(embeddings,  category){
  require(scPOP)
  metadata <- data.frame(category=category)
  lisi_scores_pro <- lisi(embeddings, meta_data = metadata, 'category')
  mean(lisi_scores_pro$category)
}

evaluate_DR_cLISI <- function(hZList, yList){
  cLISI_pro_hV <- rep(NA, length(yList))
  for(r in 1:length(yList)){
    #r <- 1
    cLISI_pro_hV[r] <- lisi_avg_scores(hZList[[r]], yList[[r]])
  }
  return(cLISI_pro_hV)
  
}


F1_score_silho <- function(embeddings, celltype, sampleID){
  require(scPOP)
  metadata <- data.frame(celltype=celltype, sampleID = sampleID)
  sh_scores_pro <- silhouette_width(embeddings, meta.data = metadata, c('celltype', "sampleID") )
  sh_ct <- (1+sh_scores_pro[1])/2 # larger is better
  sh_si <- (1+sh_scores_pro[2])/2 # smaller is better
  f1_score <- (2* (1-sh_si)*sh_ct) / ((1-sh_si) + sh_ct)
  return(f1_score)
}

evaluate_IntegrationPF <- function(embeds, clusterVec, sampleIDVec){
  
  cLISI <- lisi_avg_scores(embeds, clusterVec)
  iLISI <- lisi_avg_scores(embeds, sampleIDVec)
  F1score <- F1_score_silho(embeds, clusterVec, sampleIDVec)
  
  return(c(cLISI=cLISI, iLISI=iLISI, F1score=F1score))
}




cluster_metric <- function(hy, y, type='ARI'){
  
  require(mclust)
  require(aricode)
  switch(type, 
         ARI= adjustedRandIndex(hy, y),
         NMI = NMI(as.vector(hy), y))
}
## compute multiple methods' clustering performance using ARI and NMI
evaluate_clusterPF_mat <- function(anno_useList, cluster_compMat){
  
  
  
  indexList <- get_indexList(anno_useList)
  
  n_compMethods <- ncol(cluster_compMat)
  
  nsample <- length(anno_useList)
  ariMat <- matrix(NA, nrow=nsample+1, ncol=n_compMethods)
  if(!is.null(colnames(cluster_compMat))){
    colnames(ariMat) <- colnames(cluster_compMat)
  }else{
    colnames(ariMat) <- paste0("Method", 1:n_compMethods)
  }
  nmiMat <- ariMat
  for(i in 1:n_compMethods){
    clusterList <- vec2list(cluster_compMat[,i], nvec = sapply(anno_useList, length))
    ariMat[, i] <- c(sapply(1:nsample, function(r) cluster_metric(clusterList[[r]], anno_useList[[r]])),
                     cluster_metric(unlist(clusterList), unlist(anno_useList)))
    nmiMat[, i] <- c(sapply(1:nsample, function(r) cluster_metric(clusterList[[r]], anno_useList[[r]], type='NMI')),
                     cluster_metric(unlist(clusterList), unlist(anno_useList), type='NMI'))
  }
  
  return(list(ariMat, nmiMat))
}

evaluate_clusterPF <- function(hyList, yList){
  
  cbind(ARI=c(sapply(1:length(yList), function(r) cluster_metric(hyList[[r]], yList[[r]])),
        cluster_metric(unlist(hyList), unlist(yList))),
        NMI=c(sapply(1:length(yList), function(r) cluster_metric(hyList[[r]], yList[[r]], type="NMI")),
        cluster_metric(unlist(hyList), unlist(yList), type="NMI")) )
  
}

evaluate_DR_PF2 <- function(hZList, yList){
  MacVec_pro_hV <- rep(NA, length(yList))
  for(r in 1:length(yList)){
    #r <- 1
    MacVec_pro_hV[r] <- get_r2_mcfadden(hZList[[r]], yList[[r]])
  }
  
  MacVec_pro_hV
}


evaluate_DR_PF <- function(hZList, yList){
  MacVec_pro_hV <- rep(NA, length(yList))
  for(r in 1:length(yList)){
    #r <- 1
    MacVec_pro_hV[r] <- get_r2_mcfadden(hZList[[r]], yList[[r]])
  }
  MacR2_all <- get_r2_mcfadden(matlist2mat(hZList), unlist(yList))
  
  c(MacVec_pro_hV, MacR2_all)
}

acc_fun <- function(y1, y2){
    n1 <- length(unique(y1))
    n2 <- length(unique(y2))
    if(n1<n2){ ## ensure n1> n2
      a <- y1
      y1 <- y2
      y2 <- a
      n1 <- length(unique(y1))
      n2 <- length(unique(y2))
    }
    cm <- as.matrix(table(Actual = y1, Predicted = y2))
    rnames <-row.names(cm)
    cnames <- colnames(cm)
    union_names <- union(rnames, cnames)
    n <- length(union_names)
    cm_new <- matrix(0, n, n)
    row.names(cm_new) <- colnames(cm_new) <- union_names
    for(r in 1:n2){
      cm_new[rnames,cnames[r]] <- cm[rnames,cnames[r]]
    }
    
    sum(diag(cm_new)) / length(y1)
  }
 
kappa_fun <- function(y1, y2){
    require(irr)
    dat <- data.frame(y1, y2)
    k_res <- kappa2(dat)
    k_res$value
  }
meanF1_fun <- function(y1, y2){
    ## y1 <- hannoList[[r]]; y2 <- annoList[[r]]
    n1 <- length(unique(y1))
    n2 <- length(unique(y2))
    if(n1<n2){ ## ensure n1> n2
      a <- y1
      y1 <- y2
      y2 <- a
      n1 <- length(unique(y1))
      n2 <- length(unique(y2))
    }
    cm <- as.matrix(table(Actual = y1, Predicted = y2))
    rnames <-row.names(cm)
    cnames <- colnames(cm)
    union_names <- union(rnames, cnames)
    n <- length(union_names)
    cm_new <- matrix(0, n, n)
    row.names(cm_new) <- colnames(cm_new) <- union_names
    for(r in 1:n2){
      cm_new[rnames,cnames[r]] <- cm[rnames,cnames[r]]
    }
    
    n = sum(cm_new) # number of instances
    nc = nrow(cm_new) # number of classes
    diag_cm = diag(cm_new) # number of correctly classified instances per class 
    rowsums = rowSums(cm_new) # number of instances per class
    colsums = colSums(cm_new) # number of predictions per class
    rowsums[rowsums==0] <- 1e-7
    colsums[colsums==0] <- 1e-7
    precision = diag_cm / colsums 
    recall = diag_cm / rowsums 
    f1 = 2 * precision * recall / (precision + recall) 
    mean(f1, na.rm=T)
  }
  
evaluate_annotationPF <- function(hannoList, annoList, 
                                  unknown_label="Unknown",
                                  remove.unknown=FALSE){
  M <- length(annoList)
  if(remove.unknown){
    for(i in 1:M){
      idx <- (hannoList[[i]] != unknown_label) & (annoList[[i]] != unknown_label)
      hannoList[[i]] <- hannoList[[i]][idx]
      annoList[[i]] <- annoList[[i]][idx]
    }
    
  }
  
  
  acc_vec <- c(sapply(1:M, function(r) acc_fun(hannoList[[r]], annoList[[r]])), 
               acc_fun(unlist(hannoList), unlist(annoList)))
  
  kappa_vec <- c(sapply(1:M, function(r) kappa_fun(hannoList[[r]], annoList[[r]])), 
                 kappa_fun(unlist(hannoList), unlist(annoList)))
  
  mF1_vec <- c(sapply(1:M, function(r) meanF1_fun(hannoList[[r]], annoList[[r]])), 
               meanF1_fun(unlist(hannoList), unlist(annoList)))
  
  list(accuracy=acc_vec, kappa=kappa_vec, meanF1=mF1_vec)
  
}



get_metrics_Evals_inSimu <- function(xList, only.combine=TRUE){
    names_metric <- names(xList)
    p <- length(xList[[1]])
    metricMat <- Reduce(cbind,xList)
    colnames(metricMat) <- names_metric
    if(only.combine){
      return(metricMat[p,])
    }else{
      return(metricMat)
    }
 }
  
# Read data ---------------------------------------------------------------


read_embryo_data <- function(dir_name, sample_name=NULL){
  count <- Matrix::readMM(paste0("./",dir_name, "/sparse_matrix.mtx"))
  obs <- data.table::fread(file=paste0("./",dir_name, "/obs.csv"))
  vars <-  data.table::fread(paste0("./",dir_name, "/var_column1.csv"))
  vars <- as.data.frame(vars)
  colnames(count) <- vars[,1]
  obs <- as.data.frame(obs)
  if(!is.null(sample_name)){
    obs[,1] <- paste0(sample_name, "_", obs[,1])
  }
  row.names(count) <- obs[,1]
  row.names(obs) <- obs[,1]
  pos <- data.table::fread(paste0("./",dir_name, "/position.csv"))
  pos <- as.data.frame(pos)
  obs$row <- pos[,1]; obs$col <- pos[,2]
  row.names(pos) <- row.names(count)
  seu_tmp2 <- Seurat::CreateSeuratObject(counts=Matrix::t(count), meta.data = obs)
  return(seu_tmp2)
}



read_embryo_metadata <- function(dir_name){
  
  pos <- data.table::fread(paste0("./",dir_name, "/position.csv"))
  pos <- as.data.frame(pos)
  
  return(pos)
}


# matrix transfer ---------------------------------------------------------
filter_gene <- function(seu, min_spots=20, assay= NULL){
  
  if(is.null(assay)) assay <- DefaultAssay(seu)
  if(sum(dim(seu[[assay]]@counts))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@counts>0)>min_spots
    return(seu[names(gene_flag[unname(gene_flag)]), ])
  }else if(sum(dim(seu[[assay]]@data))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@data>0)>min_spots
    return(seu[names(gene_flag[unname(gene_flag)]), ])
  }else{
    stop("filter_gene: Seuat object must provide slots count or data in assay!")
  }
  
  
}

filter_spot <- function(seu, min_feature=0, assay=NULL){ # each spots at least include 1 non-zero features
  
  if(is.null(assay)) assay <- DefaultAssay(seu)
  col_name <- paste0("nFeature_",assay)
  idx <- seu@meta.data[,col_name] > min_feature
  seu[, idx]
  # subset(seu, subset = nFeature_RNA > min_feature)
}
firstup <- function(x) {
    ## First letter use upper capital
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}
replace_geneName <- function(seu, new_geneNames, assay=NULL){
  
  if(is.null(assay)) assay <- DefaultAssay(seu)
 
  if(length(new_geneNames) != nrow(seu))
    stop("The length of `new_geneNames` must be equal to nrow of `seu`!")
  counts <- seu[[assay]]@counts
  row.names(counts) <- new_geneNames
  meta.data <- seu@meta.data
  seu <- CreateSeuratObject(counts=counts, meta.data = meta.data)
  return(seu)
}
  
transferGeneNames <- function(genelist, now_name = "ensembl", to_name="symbol",
         species="Human", Method='biomaRt'){
  
  firstup <- function(x) {
    ## First letter use upper capital
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  if(! toupper(species) %in% c("HUMAN", "MOUSE")) stop("Check species: the current version only support Human and Mouse!")
  transferredNames <- switch (toupper(species),
                              HUMAN = {
                                if(tolower(Method)=='eg.db'){
                                  require(org.Hs.eg.db)
                                  mapIds(org.Hs.eg.db, keys = genelist,
                                         keytype = toupper(now_name), column=toupper(to_name))
                                }else if(tolower(Method)=='biomart'){
                                  require(biomaRt)
                                  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                                  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                  values=genelist,mart= mart)
                                  
                                  idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
                                  G_list <- G_list[idx_in_genelst,]
                                  idx_dup <- which(duplicated(G_list$ensembl_gene_id))
                                  G_list <- G_list[-idx_dup,]
                                  row.names(G_list) <- G_list$ensembl_gene_id
                                  symbol_list <- G_list[genelist,]$hgnc_symbol
                                  
                                  symbol_list[symbol_list==''] <- NA
                                  symbol_list
                                  
                                }else{
                                  stop("Check Method: the current version only support biomaRt and eg.db!")
                                }
                                
                                
                              },
                              MOUSE= {
                                if(tolower(Method)=='eg.db'){
                                  require(org.Mm.eg.db)
                                  mapIds(org.Mm.eg.db, keys = genelist,
                                         keytype = toupper(now_name), column=toupper(to_name))
                                }else if(tolower(Method)=='biomart'){
                                  require(biomaRt)
                                  mart <- useDataset(" mmusculus_gene_ensembl", useMart("ensembl"))
                                  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                  values=genelist,mart= mart)
                                  
                                  idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
                                  G_list <- G_list[idx_in_genelst,]
                                  idx_dup <- which(duplicated(G_list$ensembl_gene_id))
                                  G_list <- G_list[-idx_dup,]
                                  row.names(G_list) <- G_list$ensembl_gene_id
                                  symbol_list <- G_list[genelist,]$hgnc_symbol
                                  
                                  symbol_list[symbol_list==''] <- NA
                                  symbol_list
                                }else{
                                  stop("Check Method: the current version only support biomaRt and eg.db!")
                                }
                                
                              }
  )
  
  if(toupper(to_name) == 'SYMBOL')
    transferredNames <- firstup(transferredNames)
  
  flag_na <- is.na(transferredNames)
  if(any(flag_na))
    transferredNames[flag_na] <- genelist[flag_na]
  
  return(transferredNames)
}


vec2list <- function(y_int, nvec){
  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}
selectIntFeatures <- function(seulist, spaFeatureList, IntFeatures=2000){
  ## This function is used for selecting common informative features
  if(length(seulist) != length(spaFeatureList)) stop("The length of suelist and spaFeatureList must be equal!")
  if(length(seulist) ==1){
    if(length(spaFeatureList[[1]]) >= IntFeatures){
      genelist <- spaFeatureList[[1]][1:IntFeatures]
    }else{
      genelist <- spaFeatureList[[1]]
      warning("The IntFeature is larger than the  number of elements in FeatureList!")
    }
    return(genelist)
  } 
  if(any(sapply(spaFeatureList, length)< IntFeatures))
    stop("Feature list exists number of features less than IntFeatures!")
  geneUnion <- unique(unlist(lapply(spaFeatureList, function(x) x[1:IntFeatures])))
  ## ensure each seuobject has the genes in geneUnion
  gene_delete <- unique(unlist(lapply(seulist, function(x) setdiff(geneUnion, row.names(x)))))
  geneUnion <- setdiff(geneUnion, gene_delete)
  
  
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(pbapply::pblapply(seulist, function(x) 
    geneUnion[Matrix::rowSums(x@assays$RNA@counts[geneUnion,])==0])))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  
  # sort by number of datasets that identified this gene as SVG.
  nsample <- length(seulist)
  numVec <- rep(0, length(gene_Var))
  rankMat <-matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], spaFeatureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(spaFeatureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
    
  }
  
  cutNum <- sort(numVec, decreasing = T)[min(IntFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(IntFeatures, length(numVec)) - length(genelist1)
  
  gene2 <- gene_Var[numVec==cutNum]
  ### select top 2000 genes that rank 
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)
  
  return(genelist)
}

approxPCA <- function(X, q){ ## speed the computation for initial values.
  require(irlba) 
  n <- nrow(X)
  svdX  <- irlba(A =X, nv = q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}
mat2list <- function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

matlist2mat <- function (XList) 
{
  r_max <- length(XList)
  X0 <- XList[[1]]
  if (r_max > 1) {
    for (r in 2:r_max) {
      X0 <- rbind(X0, XList[[r]])
    }
  }
  return(X0)
}


get_indexList <- function(alist){
  nsample <- length(alist)
  nr <- 0
  indexList <- list()
  for(i in 1:nsample){
    if(is.matrix(alist[[i]])){
      indexList[[i]] <- (nr+1):(nrow(alist[[i]] )+nr)
      nr <- nr + nrow(alist[[i]] )
    }else{
      indexList[[i]] <- (nr+1):(length(alist[[i]] )+nr)
      nr <- nr + length(alist[[i]] )
    }
    
    
  }
  return(indexList)
}

get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}
drLouvain <- function(hZ, resolution=0.8){
  ### Louvain cluster based on estimated integrative low-dimensional embeddings. 
  
  require(Seurat)
  n <- nrow(hZ); q <- ncol(hZ)
  row.names(hZ) <- paste0("spot", 1:n)
  colnames(hZ) <- paste0("gene", 1:q)
  seu <- CreateSeuratObject(counts= t(hZ), assay='RNA')
  DefaultAssay(seu) <- "RNA"
  pca1 <- CreateDimReducObject(embeddings = hZ, key = "PC")
  seu@reductions$"pca" <- pca1
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:q)
  seu <- FindClusters(seu, resolution = resolution)
  return(seu$seurat_clusters)
}




addnames <- function(XList){
  
  for(j in 1:length(XList)){
    row.names(XList[[j]]) <- paste0("cell",j,"_", 1:nrow(XList[[j]]))
    colnames(XList[[j]]) <- paste0("gene", 1:ncol(XList[[j]]))
    
  }
  return(XList)
}

get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}


#----------------------compared methods:

cca_mnn_seurat <- function(XList, q=15, flag_LogCount=FALSE){
  ### Compare with CCA + Louvain in Seurat
  
  XList1 <- addnames(XList)
  require(Seurat)
  if(flag_LogCount){
    assayList <- lapply(XList1, function(x) CreateAssayObject(data=t(x)))
    seuList <- list()
    for(r in 1:length(XList)){
      tmpMat <- t(XList1[[r]])
      seu_tmp <- CreateSeuratObject(counts= tmpMat, assay='tmp')
      seu_tmp[['RNA']] <- assayList[[r]]
      DefaultAssay(seu_tmp) <- "RNA"
      seuList[[r]] <- seu_tmp
      
    }
  }else{
    seuList <- list()
    for(r in 1:length(XList)){
      tmpMat <- t(XList1[[r]]) / 100
      seu_tmp <- CreateSeuratObject(counts= tmpMat, assay='RNA')
      seuList[[r]] <- seu_tmp
      
    }
  }
  
  # FindVariableFeatures(seuList[[1]], nfeatures = 5000,selection.method='dispersion')
  
  seuList13 <- lapply(X = seuList, FUN = function(x) {
    x <- FindVariableFeatures(x, nfeatures = 5000,selection.method='dispersion' )
  })
  rm(seuList)
  
  
  
  features <- SelectIntegrationFeatures(object.list = seuList13)
  immune.anchors <- FindIntegrationAnchors(object.list = seuList13, anchor.features = features)
  # this command creates an 'integrated' data assay
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  rm(seuList13)
  sampleID <- get_sampleID(XList)
  immune.combined$sample <- sampleID
  DefaultAssay(immune.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = q, verbose = FALSE)
  hZ_suerat <- immune.combined[['pca']]@cell.embeddings
  
  return(list(hZ=hZ_suerat,  exprs=immune.combined[['integrated']]@data ))
}




SeuratV3Cluster <- function(hZ){
  
  n <- nrow(hZ); q <- ncol(hZ)
  row.names(hZ) <- paste0("spot", 1:n)
  colnames(hZ) <- paste0("gene", 1:q)
  
  seu <- CreateSeuratObject(counts= t(hZ), assay='RNA')
  DefaultAssay(seu) <- "RNA"
  pca1 <- CreateDimReducObject(embeddings = hZ, key = "PC")
  seu@reductions$"pca" <- pca1
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:q)
  seu <- FindClusters(seu)
  return(seu$seurat_clusters)
}

PCA_stat <- function(X, q){
  require(stats)
  n <- nrow(X)
  princ <- princomp(X)
  PCs <- princ$scores[,1:q]
  loadings <- princ$loadings[,1:q]
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}


SpatialPCA_run <- function(XList, posList,  SpatialPCnum=15, fast=TRUE, bandwidthtype = c("Silverman","SJ"), flag_LogCount=FALSE){
    addnames <- function(XList){
      
      for(j in 1:length(XList)){
        row.names(XList[[j]]) <- paste0("cell",j,"_", 1:nrow(XList[[j]]))
        colnames(XList[[j]]) <- paste0("gene", 1:ncol(XList[[j]]))
        
      }
      return(XList)
    }
	bandwidthtype <- match.arg(bandwidthtype)
    XList <- addnames(XList)
    Xmat <- t(matlist2mat(XList))
    
    require(SpatialPCA)
    location <- matlist2mat(posList) ## location should plus a big number for each sample.
    row.names(location) <- colnames(Xmat)
    # here the column names of countmat and rownames of location should be matched
    if(flag_LogCount){
      stereo_seq = CreateSpatialPCAObject(counts=abs(Xmat), 
                                          location=location, project = "SpatialPCA",
                                          customGenelist=row.names(Xmat),
                                          min.loctions = 0, min.features=0)
      stereo_seq@normalized_expr <- Xmat
    }else{
      stereo_seq = CreateSpatialPCAObject(counts=Xmat, 
                                          location=location, project = "SpatialPCA",
                                          customGenelist=row.names(Xmat),
                                          min.loctions = 0, min.features=0)
    }
    
    
    stereo_seq = SpatialPCA_buildKernel(stereo_seq, kerneltype="gaussian", bandwidthtype=bandwidthtype,
                                        bandwidth.set.by.user=NULL,sparseKernel=TRUE,
                                        sparseKernel_tol=1e-20,sparseKernel_ncore=10)
    stereo_seq = SpatialPCA_EstimateLoading(stereo_seq,fast=fast,SpatialPCnum=SpatialPCnum)
    stereo_seq = SpatialPCA_SpatialPCs(stereo_seq, fast=fast)
    
    VList_SpaPCA <- mat2list(t(as.matrix(stereo_seq@SpatialPCs)), nvec = sapply(XList, nrow))
    
    return(VList_SpaPCA)
  }
 
 
 
Run_BASS <- function(cntList, posList, R, nPC=15,  C=20, doLogNormalize=TRUE, scaleFeature=FALSE){
  
  
  require(BASS)
  tic_bass_raw <- proc.time()
  set.seed(0)
  # Set up BASS object
  BASS_raw <- createBASSObject(X=cntList, xy=posList, C = C, R = R,
                               beta_method = "SW", init_method = "mclust", 
                               nsample = 10000)
  BASS_raw <- BASS.preprocess(BASS_raw, doLogNormalize = doLogNormalize,
                              geneSelect = "hvgs", doPCA = TRUE, 
                              scaleFeature = scaleFeature, nPC = nPC)
  # Run BASS algorithm
  BASS_raw <- BASS.run(BASS_raw)
  BASS_raw <- BASS.postprocess(BASS_raw)
  zlabels_raw <- BASS_raw@results$z # spatial domain labels
  
  toc_bass_raw <- proc.time()
  time_bass_raw <- toc_bass_raw[3] - tic_bass_raw[3]
  
  BASS_raw@results$elapsed_time <- time_bass_raw
  return(BASS_raw)
}

get_fastMNN <- function(XList, q= 15){
  ### Compare with fastMNN
  require(batchelor)
  XList1 <- lapply(XList, t)
  out <- fastMNN(XList1, d= q)
  # Corrected values for use in clustering, etc.
  #str(reducedDim(out)) 
  # Extracting corrected expression values for gene 10.
  # summary(assay(out)[10,])
  list(hZ = reducedDim(out), exprs = assay(out))
}

liger_run <- function(countList, q=15, name="sample", thresh = 1e-06, max.iters = 30){
  require(rliger)
  for(r in 1:length(countList))
    colnames(countList[[r]]) <- paste0(name, r, colnames(countList[[r]]))
  names(countList) <- paste0(name, 1:length(countList))
  
  tic <- proc.time()
  osm.liger <- createLiger(countList, remove.missing = F)
  osm.liger <- rliger:::normalize(osm.liger)
  ## We have used top 2000 HVGs, thus there is no need to select genes 
  ## osm.liger <- selectGenes(osm.liger, unshared = TRUE, unshared.datasets = list(2), unshared.thresh= 0.4)
  osm.liger@var.genes <- row.names(countList[[1]])
  osm.liger <- scaleNotCenter(osm.liger)
  osm.liger <- optimizeALS(osm.liger, k = q, thresh=thresh, max.iters=max.iters)
  osm.liger <- quantile_norm(osm.liger)
  toc <- proc.time()
  time_liger <- toc[3] - tic[3]
  
  # osm.liger <- louvainCluster(osm.liger)
  # cluster_liger <- osm.liger@clusters
  hZ_liger <- osm.liger@H.norm
  return(list(hZ=mat2list(hZ_liger, nvec = sapply(countList, ncol)), time_used = time_liger))
}



dat <- NULL
get_r2_mcfadden <- function(embeds, y){
  library(nnet)
  library(performance)
  y <- as.numeric(as.factor(y))
  hq <- ncol(embeds)
  dat <- as.data.frame(cbind(y=y, x=embeds))
  dat$y <- factor(dat$y)
  name <-  c('y', paste0('V', 1:hq))
  names(dat) <-name
  formu <- paste0("y~")
  for(i in 1:hq){
    if(i < hq){
      formu <- paste(formu, name[i+1], seq='+')
    }else{
      formu <- paste(formu, name[i+1], seq='')
    }
    
  }
  model1 <- nnet::multinom(as.formula(formu), data = dat)
  R2 <- r2_mcfadden(model1)
  return(R2$R2_adjusted)
}



### Functions for plots

barPlot_real <- function(vec, ylabel='ARI', cols=NULL,...){
  require(ggplot2)
  
  
  ## filter vec
  N <- length(vec)
  vec_use <-vec[!is.na(vec)]
  
  
  df_use <- data.frame(value=vec_use, 
                       Method=names(vec_use))
  df_use$Method <- factor(df_use$Method, levels=names(vec_use))
  
  
  
  
  ## CCor
  p1 <- ggplot(df_use, aes(x=Method, y=value, fill=Method)) + 
    geom_bar(position = "dodge", stat="identity",width = 1, ...) + # , ...
    #geom_errorbar( aes(ymin=value-sd, ymax=value+sd), width=0.4, colour="orange",  size=1.3, position=position_dodge(.9)) + 
    #facet_grid(beta~Error , scales="fixed",labeller = label_bquote(beta == .(beta))) 
    labs(y=ylabel, x=NULL)+ 
    scale_x_discrete(breaks = NULL) 
  
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
}

library(ggplot2)
theme_classic_me <- theme_classic() + theme(text=element_text(size=20))


mytheme_graybox <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                             base_rect_size = base_size/22, border_color = 'gray10', bg_fill='white')
{
  half_line <- base_size/2
  t <- theme(panel.background = element_rect(fill = bg_fill,
                                             colour = NA), panel.border = element_rect(fill = NA,
                                                                                       colour = border_color),
             #line = element_blank(), #rect = element_blank(),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE,
             panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t
}


plot_RGB <- function(position, embed_3d, pointsize=2,textsize=15){
  
  # suppressMessages(require(ggplot2))
  
  info = as.data.frame(position)
  colnames(info) = c("sdimx","sdimy")
  
  
  r = (embed_3d[,1]-min(embed_3d[,1]))/(max(embed_3d[,1])-min(embed_3d[,1]))
  g = (embed_3d[,2]-min(embed_3d[,2]))/(max(embed_3d[,2])-min(embed_3d[,2]))
  b = (embed_3d[,3]-min(embed_3d[,3]))/(max(embed_3d[,3])-min(embed_3d[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
    geom_point(size=pointsize) +
    scale_color_identity()+
    theme_void()+
    theme(plot.title = element_text(size = textsize),
          text = element_text(size = textsize),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 22) ,
          legend.position = "bottom")
  
  p1
}

get_trans_colors <- function(color, num=2){
  require(colorspace)
  alphavec <- seq(0.2, 1, length=num)
  color_vec <- rep(NA, num)
  for(i in 1:num)
    color_vec[i] = adjust_transparency(color,   alpha = alphavec[i])
  return(color_vec)
}

colorbar_adj_transparent <- function(colors, alpha=0.6, plot=T){
  require(colorspace)
  ramp.list = adjust_transparency(colors,   alpha = alpha)
  print(ramp.list)
  if(plot==T){
    barplot(rep(1, length(ramp.list)), axes = FALSE, space = 0, col = ramp.list)
  }
  return(ramp.list)
}