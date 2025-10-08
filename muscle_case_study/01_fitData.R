rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula")
sapply(pkgs, require, character.only = TRUE)

setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/sc_devil/devilCaseStudies/muscle_case_study/")
source("utils/utils.R")

set.seed(12345)

## Input data
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_rna_adj.RDS"


if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

### RNA analysis 

#time <- dplyr::tibble()
m <- 'devil'
for (m in c("devil", "nebula", 'glmGamPoi')) {
 # s <- Sys.time()
  de_res <- perform_analysis_rna(input_data, method = m)
  #e <- Sys.time()
  saveRDS(de_res, paste0('results/', dataset_name, '/new/', m, '_rna', '.RDS'))
  #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  #print(time)
}



