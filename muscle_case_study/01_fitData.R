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
  de_res <- perform_analysis_rna(input_data, method = m, design_type = "age_only")

  # age effect + cell_type
  #de_res <- perform_analysis_rna(input_data, method = m, design_type = "age_plus_celltype")


  # Test age effect controlling for cell type
  #de_res <- perform_analysis_rna(input_data, method = m, design_type = "interaction", cell_type_of_interest = "Type II")

  #e <- Sys.time()
  saveRDS(de_res, paste0('results/', dataset_name, '/full/', m, '_age_rna', '.RDS'))
  #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  #print(time)
}



