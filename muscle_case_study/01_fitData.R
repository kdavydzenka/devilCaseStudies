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


### RNA analysis ###

#time <- dplyr::tibble()
m <- 'devil'
for (m in c("devil", "glmGamPoi", "nebula")) {
 # s <- Sys.time()

  #balanced_input <- subsample_balanced_cells(input_data)

  # age effect across cells
  #res_fit_age <- fit_de(input_data, method = m, design_type = "age_only")

  # Interaction - test age effect controlling for cell type
  res_fit_interaction <- fit_de(input_data, method = m, design_type = "interaction")
  
  path_fit_res <- '/orfeo/cephfs/scratch/cdslab/kdavydzenka/sc_devil/results/muscle/new/'
  
  saveRDS(res_fit_interaction, paste0(path_fit_res, m, '_fit_interaction', '.RDS'))
  
  #saveRDS(res_fit_interaction, paste0('results/', dataset_name, '/subsampled/', m, '_age_cellType_rna', '.RDS'))

  #e <- Sys.time()
  #saveRDS(de_res, paste0('results/', dataset_name, '/subsampled/', m, '_age_cellType_rna', '.RDS'))
  #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  #print(time)
}



