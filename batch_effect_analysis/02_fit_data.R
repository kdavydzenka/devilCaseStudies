
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
})

set.seed(1)
source("utils_methods.R")
source("utils_sim.R")

param_grid = readRDS("data/param_grid.rds")
idx  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
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

truth <- sim$truth_treatment %>% 
  dplyr::select(gene, is_de_treatment) %>% distinct() %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(is_de_treatment = sum(is_de_treatment) > 0)

bad_genes = rowMeans(Y) <= .1
Y = Y[!bad_genes,]

truth = truth %>% dplyr::filter(gene %in% names(bad_genes)[!bad_genes])
Y = as.matrix(Y)

INIT_OV = NULL
OV_TYPE = "MOM"

d1 = run_devil_grid(Y, meta, include_batch = c(FALSE), size_factor_methods = c(NA), 
               cluster_by = c("none"), parallel.cores = 1, 
               init_overdispersion = INIT_OV, overdispersion = OV_TYPE,
               verbose_fit = TRUE) %>% na.omit()

d2 = run_devil_grid(Y, meta, include_batch = c(FALSE), size_factor_methods = c(NA), 
                    cluster_by = c("patient"), parallel.cores = 1,
                    init_overdispersion = INIT_OV, overdispersion = OV_TYPE,
                    verbose_fit = TRUE) %>% na.omit()

d3 = run_devil_grid(Y, meta, include_batch = c(TRUE), size_factor_methods = c(NA), 
                    cluster_by = c("patient"), parallel.cores = 1,
                    init_overdispersion = INIT_OV, overdispersion = OV_TYPE,
                    verbose_fit = TRUE) %>% na.omit()

res_all = dplyr::bind_rows(d1, d2, d3)

res = res_all %>% 
  dplyr::left_join(truth, by = "gene")

saveRDS(res, paste0("results/sim_", idx, ".rds"))

res %>%
  dplyr::rename(is_de = is_de_treatment, name = method) %>%
  dplyr::group_by(name) %>%
  mutate(
    predicted = p_adj <= 0.05,
    TP = as.numeric(is_de & predicted),    # True Positive
    TN = as.numeric(!is_de & !predicted),  # True Negative
    FP = as.numeric(!is_de & predicted),   # False Positive
    FN = as.numeric(is_de & !predicted),    # False Negative
  ) %>%
  summarise(
    TP = sum(TP),
    TN = sum(TN),
    FP = sum(FP),
    FN = sum(FN),
    numerator = (TP * TN) - (FP * FN),
    denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    MCC = ifelse(denominator == 0, 0, numerator / denominator),
    TPR = TP / (TP + FN),
    FDR = FP / (TP + FP),
    FPR = FP / sum(!is_de)
  ) %>%
  dplyr::select(name, MCC, TPR, FDR, FPR) %>%
  dplyr::arrange(-MCC)

res %>%
  dplyr::rename(is_de = is_de_treatment, name = method) %>%
  dplyr::group_by(name) %>%
  mutate(
    predicted = p_adj <= 0.05,
    TP = as.numeric(is_de & predicted),    # True Positive
    TN = as.numeric(!is_de & !predicted),  # True Negative
    FP = as.numeric(!is_de & predicted),   # False Positive
    FN = as.numeric(is_de & !predicted),    # False Negative
  ) %>%
  summarise(
    TP = sum(TP),
    TN = sum(TN),
    FP = sum(FP),
    FN = sum(FN),
    numerator = (TP * TN) - (FP * FN),
    denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    MCC = ifelse(denominator == 0, 0, numerator / denominator),
    TPR = TP / (TP + FN),
    FDR = FP / (TP + FP),
    FPR = FP / sum(!is_de)
  ) %>%
  dplyr::select(name, MCC, TPR, FDR, FPR) %>%
  dplyr::arrange(-MCC)

