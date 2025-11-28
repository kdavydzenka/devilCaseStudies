# 
# library(tidyverse)
# 
# batch_regimes <- tribble(
#   ~batch_regime, ~libsize_cv, ~batch_lfc_sd, ~disp_mult_sd,
#   "low",         1e-6,          1e-6,           1e-6,
#   "high",        1e-2,          1e-2,           1e-2
# )
# 
# effect_regimes <- tribble(
#   ~effect_regime, ~lfc_treat_loc, ~lfc_treat_sd, ~de_prop_treat,
#   "weak",          5e-5,           5e-5,          0.050,
#   "strong",        1e-3,           1e-3,          0.050
# )
# 
# scale_regimes <- tribble(
#   ~scale_regime, ~n_batches, ~patients_per_batch, ~cells_per_patient, ~ngenes,
#   "small",        2,          3,                   500,               1000,
#   "large",        5,          4,                   1000,               1000
# )
# 
# heterogeneity_regimes <- tribble(
#   ~prop_clusters_respond,
#   1.0
# )
# 
# batch_treat_patterns <- tribble(
#   ~batch_treat_pattern,
#   "balanced",
#   "confounded"
# )
# 
# design_full <- tidyr::crossing(
#   batch_regimes,
#   effect_regimes,
#   scale_regimes,
#   heterogeneity_regimes,
#   batch_treat_patterns,
#   assignment = c("by_patient"),
#   link_cluster_to_treatment = TRUE,
#   de_prop_cluster = 0.05,
#   iters = 1:10
# ) %>%
#   mutate(seed = row_number())
# 
# saveRDS(design_full, "data/param_grid.rds")


# Use following version

library(tidyverse)

batch_regimes <- tribble(
  ~batch_regime, ~libsize_cv, ~batch_lfc_sd, ~disp_mult_sd,
  "low",         1e-6,          1e-6,           1e-6,
  "high",        1e-1,          1e-1,           1e-1
)

effect_regimes <- tribble(
  ~effect_regime, ~lfc_treat_loc, ~lfc_treat_sd, ~de_prop_treat,
  "weak",          4e-5,           4e-5,          0.050,
  "strong",        2e-1,           1e-1,          0.050
)

scale_regimes <- tribble(
  ~scale_regime, ~n_batches, ~patients_per_batch, ~cells_per_patient, ~ngenes,
  "small",        2,          3,                   500,               1000,
  "large",        5,          4,                   1000,               1000
)

heterogeneity_regimes <- tribble(
  ~prop_clusters_respond,
  1.0
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
  link_cluster_to_treatment = TRUE,
  de_prop_cluster = 0.05,
  iters = 1:10
) %>%
  mutate(seed = row_number())

saveRDS(design_full, "data/param_grid.rds")
