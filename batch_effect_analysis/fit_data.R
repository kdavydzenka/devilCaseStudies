
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

# source("utils_methods.R")
# source("utils_sim.R")
# 
# # --- load grid & select row (SLURM-friendly) ---
# param_grid <- readRDS("data/param_grid.rds")   # FIX: was .csv
# idx        <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
# stopifnot(idx >= 1L, idx <= nrow(param_grid))
# row        <- param_grid[idx, ]
# 
# # Optional: set the seed per-row for reproducibility of this task
# set.seed(row$seed)
# 
# # --- handle optional treated_batches (may be NA/NULL) ---
# treated_batches_arg <- NULL
# if ("treated_batches" %in% names(row)) {
#   tb <- row$treated_batches[[1]]
#   if (!is.null(tb) && !is.na(tb)) treated_batches_arg <- tb
# }
# 
# # --- simulate ---
# sim <- sim_treatment_clusters(
#   ngenes                     = row$ngenes,
#   cells_per_patient          = row$cells_per_patient,
#   n_batches                  = row$n_batches,
#   patients_per_batch         = row$patients_per_batch,
#   # batch nuisance
#   batch_lfc_sd               = row$batch_lfc_sd,
#   libsize_cv                 = row$libsize_cv,
#   disp_mult_sd               = row$disp_mult_sd,
#   # assignment & confounding
#   assignment                 = row$assignment,
#   confounding_mode           = row$confounding_mode,
#   treated_batches            = treated_batches_arg,
#   # cluster/treatment linkage
#   link_cluster_to_treatment  = row$link_cluster_to_treatment,
#   # treatment DE
#   de_prop_treat              = row$de_prop_treat,
#   lfc_treat_loc              = row$lfc_treat_loc,
#   lfc_treat_sd               = row$lfc_treat_sd,
#   lfc_signed                 = row$lfc_signed,
#   # optional cluster markers
#   de_prop_cluster            = row$de_prop_cluster,
#   # dispersion regimes
#   phi0                       = row$phi0,
#   gene_batch_dispersion      = row$gene_batch_dispersion,
#   # seed (already set above, but pass for clarity)
#   seed                       = row$seed
# )

Y    <- sim$counts
meta <- sim$meta
truth <- sim$truth_treatment %>% 
  dplyr::select(gene, is_de_treatment) %>% distinct()

bad_genes = rowSums(Y) == 0
Y = Y[!bad_genes,]
truth = truth %>% dplyr::filter(gene %in% names(bad_genes)[!bad_genes])

# res_all <- run_devil_grid(
#   Y, meta,
#   include_batch       = c(TRUE, FALSE),
#   size_factor_methods = c(NA, "psinorm"),  # NA -> none
#   cluster_by          = c("none", "patient"),
#   parallel.cores      = 1,   # tweak for your machine
#   verbose_fit         = TRUE, 
#   overdispersion      = "new"
# )

d1 = run_devil_grid(Y, meta, include_batch = c(FALSE), size_factor_methods = c(NA), 
               cluster_by = c("none"), parallel.cores = 1,
               verbose_fit = TRUE)

d2 = run_devil_grid(Y, meta, include_batch = c(FALSE), size_factor_methods = c(NA), 
                    cluster_by = c("patient"), parallel.cores = 1,
                    verbose_fit = TRUE)

d3 = run_devil_grid(Y, meta, include_batch = c(FALSE), size_factor_methods = c("psinorm"), 
                    cluster_by = c("patient"), parallel.cores = 1,
                    verbose_fit = TRUE)

d4 = run_devil_grid(Y, meta, include_batch = c(TRUE), size_factor_methods = c("psinorm"), 
                    cluster_by = c("patient"), parallel.cores = 1,
                    verbose_fit = TRUE)

res_all = dplyr::bind_rows(d1, d2, d3, d4)

# Fit other methods
#edgeR_batch_res = run_edgeR_with_batch(Y, meta)
#edgeR_batch_res_no = run_edgeR_wo_batch(Y, meta)
#nebula_batch_res = run_nebula_w_batch(Y, meta)
#nebula_batch_res_no = run_nebula_wo_batch(Y, meta)
#res_all = res_all %>% dplyr::filter(grepl("devil", method))

# funcs = list(run_edgeR_with_batch, run_edgeR_wo_batch, run_nebula_w_batch, run_nebula_wo_batch)
# 
funcs = list(run_edgeR_with_batch)
additional_res = lapply(1:length(funcs), function(i){
	funcs[[i]](Y, meta)
}) %>% do.call(bind_rows, .)


#additional_res = dplyr::bind_rows(
#  edgeR_batch_res,
#  edgeR_batch_res_no#,
#  #nebula_batch_res,
#  #nebula_batch_res_no
#)

res_all = dplyr::bind_rows(res_all, additional_res)

res = res_all %>% 
  dplyr::left_join(truth, by = "gene")

saveRDS(res, paste0("results/sim_", idx, ".rds"))

res %>%
  dplyr::rename(is_de = is_de_treatment, name = method) %>% 
  dplyr::group_by(name) %>% 
  mutate(
    predicted = p_adj <= 0.05 & abs(lfc) > .5,
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

