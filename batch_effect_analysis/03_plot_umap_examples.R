
rm(list = ls())
# -------------------- setup --------------------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
  library(Matrix)
  library(tidyverse)
  library(devil)
  library(edgeR)
  # library(glmGamPoi)
  # library(limma)
  library(nebula)
  library(Seurat)
})

set.seed(1)
source("utils_methods.R")
source("utils_sim.R")

param_grid = readRDS("data/param_grid.rds")
idxs = param_grid %>% dplyr::mutate(idx = row_number()) %>% 
  dplyr::select(batch_regime, effect_regime, scale_regime, idx)  %>% 
  dplyr::group_by(batch_regime, effect_regime, scale_regime) %>% 
  dplyr::slice_head(n = 1) %>% 
  dplyr::pull(idx)

idx = 141

for (idx in idxs) {
  print(idx)
  
  row = param_grid[idx,]
  
  sim = sim_treatment_clusters(
    ngenes                   = row$ngenes,
    cells_per_patient        = row$cells_per_patient,
    n_batches                = row$n_batches,
    patients_per_batch       = row$patients_per_batch,
    batch_lfc_sd             = row$batch_lfc_sd,
    libsize_cv               = row$libsize_cv,
    disp_mult_sd             = row$disp_mult_sd,
    assignment               = row$assignment,
    link_cluster_to_treatment= row$link_cluster_to_treatment,
    de_prop_treat            = row$de_prop_treat,
    lfc_treat_loc            = row$lfc_treat_loc,
    lfc_treat_sd             = row$lfc_treat_sd,
    de_prop_cluster          = row$de_prop_cluster,
    seed                     = row$seed
  )
  
  Y    <- sim$counts
  meta <- sim$meta
  Y = as.matrix(Y)
  
  # -------------------- Seurat UMAP --------------------
  # make sure we have cell IDs to match Y
  # assume meta has a column called "cell" that matches colnames(Y)
  if ("cell" %in% colnames(meta)) {
    # reorder meta to match the columns of Y
    meta <- meta[match(colnames(Y), meta$cell), ]
    stopifnot(all(meta$cell == colnames(Y)))
    rownames(meta) <- meta$cell
  } else {
    # fallback: enforce a shared ID if none provided
    colnames(Y) <- paste0("cell_", seq_len(ncol(Y)))
    rownames(meta) <- colnames(Y)
  }
  
  # Create Seurat object with aligned metadata
  seu <- CreateSeuratObject(
    counts   = Y,
    meta.data = meta
  )
  
  seu$batch = meta$batch
  seu$patient = meta$patient
  seu$condition = meta$condition
  
  # Standard Seurat pipeline
  # seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  # seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  # seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
  # seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  #seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
  seu = RunTSNE(seu)
  
  # Extract UMAP embeddings and combine with metadata
  umap_df <- Embeddings(seu, "tsne") %>%
    as.data.frame()
  colnames(umap_df) <- c("TSNE-1", "TSNE-2")
  
  plot_df <- seu@meta.data %>%
    dplyr::bind_cols(umap_df)
  
  if (row$scale_regime == "small") plot_df = plot_df %>% dplyr::slice_sample(prop = .5)
  if (row$scale_regime == "large") plot_df = plot_df %>% dplyr::slice_sample(prop = .1)
  
  p_umap <- ggplot(plot_df, aes(x = `TSNE-1`, y = `TSNE-2`, 
                                colour = condition, shape  = batch), 
                   environment = emptyenv()) +
    geom_point(size = 2) +
    theme_bw()
  
  dir.create("results/umaps", recursive = T)
  saveRDS(p_umap, 
          file = paste0("results/umaps/Batch_", row$batch_regime, "_effect_", row$effect_regime, "_scale_", row$scale_regime, ".RDS"))
  
}

umaps = lapply(idxs, function(idx) {
  row = param_grid[idx,]
  
  sim = sim_treatment_clusters(
    ngenes                   = row$ngenes,
    cells_per_patient        = row$cells_per_patient,
    n_batches                = row$n_batches,
    patients_per_batch       = row$patients_per_batch,
    batch_lfc_sd             = row$batch_lfc_sd,
    libsize_cv               = row$libsize_cv,
    disp_mult_sd             = row$disp_mult_sd,
    assignment               = row$assignment,
    link_cluster_to_treatment= row$link_cluster_to_treatment,
    de_prop_treat            = row$de_prop_treat,
    lfc_treat_loc            = row$lfc_treat_loc,
    lfc_treat_sd             = row$lfc_treat_sd,
    de_prop_cluster          = row$de_prop_cluster,
    seed                     = row$seed
  )
  
  Y    <- sim$counts
  meta <- sim$meta
  Y = as.matrix(Y)
  
  # -------------------- Seurat UMAP --------------------
  # make sure we have cell IDs to match Y
  # assume meta has a column called "cell" that matches colnames(Y)
  if ("cell" %in% colnames(meta)) {
    # reorder meta to match the columns of Y
    meta <- meta[match(colnames(Y), meta$cell), ]
    stopifnot(all(meta$cell == colnames(Y)))
    rownames(meta) <- meta$cell
  } else {
    # fallback: enforce a shared ID if none provided
    colnames(Y) <- paste0("cell_", seq_len(ncol(Y)))
    rownames(meta) <- colnames(Y)
  }
  
  # Create Seurat object with aligned metadata
  seu <- CreateSeuratObject(
    counts   = Y,
    meta.data = meta
  )
  
  seu$batch = meta$batch
  seu$patient = meta$patient
  seu$condition = meta$condition
  
  # Standard Seurat pipeline
  # seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  # seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  # seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
  # seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  #seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
  seu = RunTSNE(seu)
  
  # Extract UMAP embeddings and combine with metadata
  umap_df <- Embeddings(seu, "tsne") %>%
    as.data.frame()
  colnames(umap_df) <- c("TSNE-1", "TSNE-2")
  
  plot_df <- seu@meta.data %>%
    dplyr::bind_cols(umap_df)
  
  if (row$scale_regime == "small") plot_df = plot_df %>% dplyr::slice_sample(prop = .5)
  if (row$scale_regime == "large") plot_df = plot_df %>% dplyr::slice_sample(prop = .1)
  
  p_umap <- ggplot(plot_df, aes(x = `TSNE-1`, y = `TSNE-2`, 
                                colour = condition, shape  = batch), 
                   environment = emptyenv()) +
    geom_point(size = 2) +
    theme_bw()
  
  dir.create("results/umaps", recursive = T)
  saveRDS(p_umap, 
          file = paste0("results/umaps/Batch_", row$batch_regime, "_effect_", row$effect_regime, "_scale_", row$scale_regime, ".RDS"))
})
