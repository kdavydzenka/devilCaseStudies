rm(list=ls())
pkgs <- c("ggplot2", "dplyr", "tidyr", "tibble", "reshape2", "Seurat", "glmGamPoi", "devil", "nebula", "Matrix")
sapply(pkgs, require, character.only = TRUE)

setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/sc_devil/devilCaseStudies/muscle_case_study/")
source("utils/utils.R")

set.seed(12345)

## Input data
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_counts_rna.RDS"


if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

# Load and prepare input
input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

# Define parameters
methods <- c("nebula")
design_tests <- c("age_type1", "age_type2", "interaction")
path_fit_res <- "/orfeo/cephfs/scratch/cdslab/kdavydzenka/sc_devil/results/muscle/new/"

# Run DE tests
for (m in methods) {
  for (dt in design_tests) {
    
    message("Running DE test for method: ", m, " | design_test: ", dt)
    
    # Load model fit
    fit <- readRDS(paste0(path_fit_res, m, "_fit_interaction",  ".RDS"))
    
    # Run DE test
    de_res <- de_test(
      input_data,
      fit,
      method = m,
      design_test = dt
    )
    
    # Save results
    save_path <- paste0("results/", dataset_name, "/full/")
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    
    saveRDS(de_res, paste0(save_path, m, "_", dt, ".RDS"))
  }
}
