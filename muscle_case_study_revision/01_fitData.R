
rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 01_fitData.R METHOD DESIGN_TEST SUBSAMPLE(TRUE/FALSE)")
}

m  <- args[1]             # "devil", "glmGamPoi", "nebula"
dt <- args[2]             # "age_only", "interaction"
sub <- as.logical(args[3]) # TRUE / FALSE

if (dt == "age_only") {
  contrasts = c("age_only")
} else {
  contrasts = c("age_type1", "age_type2", "interaction")
}

pkgs <- c("ggplot2", "dplyr", "tidyr", "tibble", "reshape2", "Seurat", "glmGamPoi", "devil", "nebula", "Matrix")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils.R")

set.seed(12345)

## Input data
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_counts_rna.RDS"

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name), recursive = T)
}

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

# # for test let's work small
# gene.idxs = sample(1:nrow(input_data$counts), 1000)
# cell.idxs = sample(1:ncol(input_data$counts), 10000)
# 
# input_data$counts = input_data$counts[gene.idxs, cell.idxs]
# input_data$metadata = input_data$metadata[cell.idxs, ]

### Fit a single DE test ###

message("Fit coefficients for method: ", m, 
        " | design_test: ", dt, 
        " | subsampling: ", sub)

design = model.matrix(~ age_cluster * cell_type, input_data$metadata)

input_data$metadata %>% dplyr::mutate(idx = row_number()) %>% 
  dplyr::filter(cell_type == "Type II", age_pop != "old_pop")

if (sub) {
  balanced_input <- subsample_balanced_cells(input_data)
  fit <- fit_de(balanced_input, method = m, design_type = dt)
  save_path <- paste0("results/", dataset_name, "/fits_sub/")
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  saveRDS(fit,     paste0(save_path, m, "_", dt, ".RDS"))
  
  for (cont in contrasts) {
    de_res <- de_test(
      balanced_input,
      fit,
      method = m,
      design_test = cont
    )
    
    save_path <- paste0("results/", dataset_name, "/subsampled/")
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    saveRDS(de_res, paste0(save_path, m, "_", cont, ".RDS"))  
  }
  
} else {
  fit <- fit_de(input_data, method = m, design_type = dt)
  save_path <- paste0("results/", dataset_name, "/fits/")
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  saveRDS(fit,     paste0(save_path, m, "_", dt, ".RDS"))
  
  for (cont in contrasts) {
    de_res <- de_test(
      input_data,
      fit,
      method = m,
      design_test = cont
    )
    
    save_path <- paste0("results/", dataset_name, "/full/")
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    saveRDS(de_res, paste0(save_path, m, "_", cont, ".RDS"))  
  }
}
