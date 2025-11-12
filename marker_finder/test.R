
setwd("~/Desktop/devil_reviews/marker_finder")
rm(list = ls())
library(Seurat)
library(tidyverse)

d1 = ReadMtx(mtx = "data/BE1run12/HTB178/matrix.mtx.gz", 
            features = "data/BE1run12/HTB178/features.tsv.gz", 
            cells = "data/BE1run12/HTB178/barcodes.tsv.gz")

d2 = ReadMtx(mtx = "data/BE1run12/DV90/matrix.mtx.gz", 
             features = "data/BE1run12/DV90/features.tsv.gz", 
             cells = "data/BE1run12/DV90/barcodes.tsv.gz")

d = cbind(d1, d2)

seurat_object <- CreateSeuratObject(counts = d)
NPC = 20
cluster_res=0.1

seurat_object <- NormalizeData(seurat_object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = NPC) %>%
  FindNeighbors(dims = 1:NPC) %>%
  FindClusters(resolution = cluster_res) %>%
  RunUMAP(dims=1:NPC)

new_clusters <- as.numeric(seurat_object$seurat_clusters)
new_clusters <- factor(new_clusters, levels = sort(unique(new_clusters)))
seurat_object$seurat_clusters <- new_clusters

DimPlot(seurat_object)
counts <- as.matrix(seurat_object@assays$RNA$counts)

whole_res <- dplyr::tibble()

c = 2
SEED = 123
for (c in unique(seurat_object$seurat_clusters)) {
  
  print(c)
  
  idx_cluster <- which(seurat_object$seurat_clusters == c)
  idx_others <- which(!(seurat_object$seurat_clusters == c))
  
  if (length(idx_others) > length(idx_cluster)) {
    set.seed(SEED)
    #idx_others <- sample(idx_others, as.integer(length(idx_others) * .05), replace = F)
    idx_others <- sample(idx_others, length(idx_cluster), replace = F)
  }
  
  cell_idx <- c(idx_cluster, idx_others)
  # clusters <- as.numeric(as.factor(seurat_object$donor))
  design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_object$seurat_clusters == c))
  
  # First filter
  dm <- design_matrix[cell_idx,]
  cc <- counts[,cell_idx]
  # clusters <- clusters[cell_idx]
  
  # Second filter
  cell_idx <- which((colSums(cc) > 20) == TRUE)
  dm <- dm[cell_idx,]
  cc <- cc[,cell_idx]
  # clusters <- clusters[cell_idx]
  
  # Third filter
  gene_idx = which((rowSums(cc) > 20) == TRUE)
  cc <- cc[gene_idx,]
  
  rownames(dm) <- colnames(cc)
  
  if (method == 'devil') {
    s <- Sys.time()
    fit <- devil::fit_devil(cc, dm, size_factors = T, overdispersion = T, init_overdispersion = 100, offset = 1e-6, verbose = TRUE, tolerance = 1e-3, max_iter = 100, parallel.cores = 1)
    e <- Sys.time()
    
    fit$beta %>% view()
    
    res <- devil::test_de(fit, contrast = c(0,1), max_lfc = 50) %>% dplyr::mutate(cluster = c)
    
    #res <- devil::test_de(fit, contrast = c(0,1), clusters = 1:length(idxs), max_lfc = Inf) %>% dplyr::mutate(cluster = c)
  } else if (method == "glmGamPoi") {
    s <- Sys.time()
    fit <- glmGamPoi::glm_gp(cc, dm, size_factors = "normed_sum", verbose = T)
    e <- Sys.time()
    #fit <- glmGamPoi::glm_gp(cc, dm, size_factors = FALSE, verbose = T)
    res <- glmGamPoi::test_de(fit, contrast = c(0,1))
    res <- res %>% dplyr::as_tibble() %>% dplyr::select(name, pval, adj_pval, lfc) %>% dplyr::mutate(cluster = c)
  } else if (method == "nebula") {
    s <- Sys.time()
    sf <- devil:::calculate_sf(cc)
    data_g = nebula::group_cell(count=cc,id=clusters,pred=dm)
    fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, ncore = 1, mincp = 0, cpc = 0, offset = sf)
    e <- Sys.time()
    #fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, ncore = 1, mincp = 0, cpc = 0)
    res <- dplyr::tibble(
      name = fit$summary$gene,
      pval = fit$summary$p_groupTRUE,
      adj_pval = p.adjust(fit$summary$p_groupTRUE, "BH"),
      lfc=fit$summary$logFC_groupTRUE
    ) %>% dplyr::mutate(cluster = c)
  } else {
    stop("method not recognized")
  }
  
  res <- res %>% dplyr::mutate(delta_time = e - s)
  whole_res <- dplyr::bind_rows(whole_res, res)
}
whole_res
}