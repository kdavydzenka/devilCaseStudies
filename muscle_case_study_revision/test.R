
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ComplexHeatmap", "magick")
sapply(pkgs, require, character.only = TRUE)
library(ggplot2)
library(patchwork)
library(grid)

cell_group_colors = c(
  "old" = "darkorange",
  "young" = "steelblue"
)

# input HEATMAP ####
source("utils/utils.R")
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_counts_rna.RDS"

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

clusters = input_data$metadata$sample

old_fit = readRDS("results/MuscleRNA/fits/devil_interaction.RDS")

old_res = readRDS("results/MuscleRNA/full/devil_interaction.RDS")

# new_res = devil::test_de(old_fit, contrast = c(0,0,0,1), clusters = clusters)

devil:::deprecated_test_de

parallel.cores = 1
devil.fit = old_fit
contrast = c(0,0,0,1)

# Detect cores to use
max.cores <- parallel::detectCores()
if (is.null(parallel.cores)) {
  n.cores <- max.cores
} else {
  if (parallel.cores > max.cores) {
    message("Requested ", parallel.cores, " cores, but only ", max.cores, " available.")
  }
  n.cores <- min(max.cores, parallel.cores)
}

# Extract necessary information
ngenes <- nrow(devil.fit$input_matrix)
nsamples <- nrow(devil.fit$design_matrix)
contrast <- as.array(contrast)

# Calculate log fold changes
lfcs <- (devil.fit$beta %*% contrast) %>%
  unlist() %>%
  unname() %>%
  c()

if (!is.null(clusters) & !is.numeric(clusters)) {
  message("Converting clusters to numeric factors")
  clusters <- as.numeric(as.factor(clusters))
}

# Calculate p-values in parallel
gene_idx = 1
p_values <- parallel::mclapply(seq_len(nrow(devil.fit$input_matrix)), function(gene_idx) {
  mu_test <- lfcs[gene_idx]
  
  b = devil:::compute_hessian(devil.fit$beta[gene_idx,], 
                              devil.fit$overdispersion[gene_idx], 
                              devil.fit$input_matrix[gene_idx,], 
                              devil.fit$design_matrix, 
                              devil.fit$size_factors)
  
  msimple = devil:::compute_meat(
    devil.fit$design_matrix, 
    devil.fit$input_matrix[gene_idx,], 
    devil.fit$beta[gene_idx,], 
    devil.fit$overdispersion[gene_idx], 
    devil.fit$size_factors
  )
  
  m = devil:::compute_clustered_meat(devil.fit$design_matrix, 
                                     devil.fit$input_matrix[gene_idx,], 
                                     devil.fit$beta[gene_idx,], 
                                     devil.fit$overdispersion[gene_idx], 
                                     devil.fit$size_factors, 
                                     clusters)
  
  S = (b %*% m %*% b) * dim(devil.fit$input_matrix)[2]
  H = (b %*% msimple %*% b) * dim(devil.fit$input_matrix)[2]
  
  total_variance <- t(contrast) %*% H %*% contrast
  p <- 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = FALSE)
  
  total_variance <- t(contrast) %*% S %*% contrast
  pnull <- 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = FALSE)
  
  max(p, pnull)
}, mc.cores = n.cores) %>% unlist()

# Create tibble with results
result_df <- dplyr::tibble(
  name = rownames(devil.fit$beta),
  pval = p_values,
  adj_pval = stats::p.adjust(p_values, method = pval_adjust_method),
  lfc = lfcs / log(2)
)

# Filter results based on max_lfc
result_df <- result_df %>%
  dplyr::mutate(lfc = ifelse(.data$lfc >= max_lfc, max_lfc, .data$lfc)) %>%
  dplyr::mutate(lfc = ifelse(.data$lfc <= -max_lfc, -max_lfc, .data$lfc))