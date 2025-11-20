
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

param_grid = readRDS("data/param_grid.csv")
idx  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
row = param_grid[idx,]

sim = sim_treatment_clusters(
  ngenes                   = row$ngenes,             # fixed 1000 here
  cells_per_patient        = row$cells_per_patient,
  n_batches                = row$n_batches,
  patients_per_batch       = row$patients_per_batch,
  batch_lfc_sd             = row$batch_lfc_sd,
  libsize_cv               = row$libsize_cv,
  disp_mult_sd             = row$disp_mult_sd,
  assignment               = row$assignment,
  link_cluster_to_treatment= row$link_cluster_to_treatment,
  de_prop_treat            = row$de_prop_treat,      # fixed 0.10 here
  lfc_treat_loc            = row$lfc_treat_loc,
  lfc_treat_sd             = row$lfc_treat_sd,
  de_prop_cluster          = row$de_prop_cluster,
  seed                     = row$seed
)


Y    <- sim$counts
meta <- sim$meta
truth <- sim$truth_treatment %>% 
  dplyr::select(gene, is_de_treatment) %>% distinct()

bad_genes = rowSums(Y) == 0
Y = Y[!bad_genes,]
truth = truth %>% dplyr::filter(gene %in% names(bad_genes)[!bad_genes])

res_all <- run_devil_grid(
  Y, meta,
  include_batch       = c(TRUE, FALSE),
  size_factor_methods = c(NA),  # NA -> none
  cluster_by          = c("none","patient"),
  parallel.cores      = 1,   # tweak for your machine
  verbose_fit         = TRUE
)

# Fit other methods
edgeR_batch_res = run_edgeR_with_batch(Y, meta)
#edgeR_batch_res_no = run_edgeR_wo_batch(Y, meta)
nebula_batch_res = run_nebula_w_batch(Y, meta)
nebula_batch_res_no = run_nebula_wo_batch(Y, meta)

additional_res = dplyr::bind_rows(
  edgeR_batch_res,
  #edgeR_batch_res_no,
  nebula_batch_res,
  nebula_batch_res_no
)

res_all = dplyr::bind_rows(res_all, additional_res)

res = res_all %>% 
  dplyr::left_join(truth, by = "gene")

saveRDS(res, paste0("results/sim_", idx, ".rds"))
