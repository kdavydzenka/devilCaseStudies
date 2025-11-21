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

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

# Define parameters
methods <- c("devil")
design_tests <- c("age_type1", "age_type2", "interaction")


### Fit an DE test ###

#time <- dplyr::tibble()

for (m in methods) {
  for (dt in design_tests) {
    #s <- Sys.time()	  
    message("Fit coefficients for method: ", m, " | design_test: ", dt)

    balanced_input <- subsample_balanced_cells(input_data)

    fit <- fit_de(balanced_input, method = m, design_type = "interaction")

    # Run DE test
    de_res <- de_test(
      balanced_input,
      fit,
      method = m,
      design_test = dt
    )

    # Save results
    save_path <- paste0("results/", dataset_name, "/subsampled/")
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

    saveRDS(de_res, paste0(save_path, m, "_", dt, ".RDS"))

    #e <- Sys.time()
    #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
    #print(time}
  }
}



