library(tidyverse)

batch_regimes <- tribble(
  ~batch_regime, ~libsize_cv, ~batch_lfc_sd, ~disp_mult_sd,
  "low",         0.10,         0.10,          0.10,
  "high",        0.40,         0.40,          0.50  # slightly stronger
)

effect_regimes <- tribble(
  ~effect_regime, ~lfc_treat_loc, ~lfc_treat_sd, ~de_prop_treat,
  "weak",          log(0.5),       0.20,          0.10,
  "strong",        log(1.2),       0.40,          0.20
)

scale_regimes <- tribble(
  ~scale_regime, ~n_batches, ~patients_per_batch, ~cells_per_patient, ~ngenes,
  "small",        2,          3,                   500,               1000,
  "large",        4,          5,                   500,               1000
)

heterogeneity_regimes <- tribble(
  ~prop_clusters_respond,
  1.0,
  0.5
)

batch_treat_patterns <- tribble(
  ~batch_treat_pattern,
  "balanced",
  "confounded"
)

design_full <- tidyr::crossing(
  batch_regimes,
  effect_regimes,
  scale_regimes,
  heterogeneity_regimes,
  batch_treat_patterns,
  assignment = c("by_patient"),
  link_cluster_to_treatment = FALSE,  # important for Regime A
  de_prop_cluster = 0.05,
  iters = 1:10
) %>%
  mutate(seed = row_number())

saveRDS(design_full, "data/param_grid.rds")
