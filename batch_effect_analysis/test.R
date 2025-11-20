
# -------------------- setup --------------------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(glmGamPoi)
  library(limma)
  library(nebula)
})

set.seed(1)
source("utils_methods.R")
source("utils_sim.R")

# One scenario
sim = sim_treatment_clusters(
  ngenes = 1000, 
  cells_per_patient = 500, 
  n_batches = 4, 
  patients_per_batch = 3, 
  batch_lfc_sd = 0.1, 
  libsize_cv = 0.4, 
  disp_mult_sd = 0.2, 
  assignment = "by_patient", 
  link_cluster_to_treatment = T,
  de_prop_treat = .1,
  lfc_treat_loc = log(.5),   # on log scale (ln FC)
  lfc_treat_sd  = 0.2,
  de_prop_cluster = 0.001,
  seed = 1
)


Y    <- sim$counts
meta <- sim$meta
truth <- sim$truth_treatment %>% 
  dplyr::select(gene, is_de_treatment) %>% distinct()

res_list <- list(
  # edgeR = run_edgeR_glmQLF(Y, meta, cl),
  # dseq2 = run_DESeq2(Y, meta, cl),
  ggp   = run_glmGamPoi(Y, meta, cl),
  devil_BS = run_DEVIL_batch_and_sandwich(Y, meta, cl), 
  devil_sP = run_DEVIL_nobatch_sandwichbatches(Y, meta, cl), 
  devil_sB = run_DEVIL_nobatch_sandwichpatients(Y, meta, cl)
)

res = dplyr::bind_rows(res_list) %>% 
  dplyr::left_join(truth, by = "gene")


res %>% 
  ggplot(mapping = aes(x = lfc, y = -log10(p_adj), col = is_de_treatment)) +
  geom_point() +
  facet_wrap(~method)

res %>% 
  dplyr::group_by(method) %>% 
  dplyr::mutate(is_pred = p_adj <= .05 & abs(lfc) > .5) %>% 
  dplyr::summarise(TP = sum(is_pred & is_de_treatment), FP = sum(is_pred & !is_de_treatment))


