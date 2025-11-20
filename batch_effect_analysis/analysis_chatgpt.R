
rm(list = ls())
library(tidyverse)

# -------------------------------------------------------------------
# 1. Load parameter grid (updated version with new columns)
# -------------------------------------------------------------------
grid <- readRDS("data/param_grid.rds") %>%
  mutate(sim_id = row_number())   # keep an explicit ID

# Quick check
print(head(grid))

# -------------------------------------------------------------------
# 2. Helper: summarise one simulation setting (one row of grid)
# -------------------------------------------------------------------
summarise_one_sim <- function(i, grid, results_dir = "results") {
  message("Processing sim ", i, " / ", nrow(grid))
  this <- grid[i, ]
  
  res_path <- file.path(results_dir, paste0("sim_", i, ".rds"))
  if (!file.exists(res_path)) {
    message("  -> results file not found, skipping: ", res_path)
    return(tibble())  # empty tibble, will be dropped in bind_rows()
  }
  
  res_raw <- readRDS(res_path)
  
  # Expect columns: gene, method, p_adj, is_de_treatment (or similar)
  # Adjust names here if needed.
  res <- res_raw %>%
    # defensively handle NAs
    mutate(
      p_adj = ifelse(is.na(p_adj), 1, p_adj),
      lfc   = ifelse(is.na(lfc), 0, lfc),
      is_de = is_de_treatment,
      name  = method
    ) %>%
    group_by(name) %>%
    mutate(
      predicted = p_adj <= 0.05,
      TP = as.numeric(is_de & predicted),
      TN = as.numeric(!is_de & !predicted),
      FP = as.numeric(!is_de & predicted),
      FN = as.numeric(is_de & !predicted)
    ) %>%
    summarise(
      TP = sum(TP),
      TN = sum(TN),
      FP = sum(FP),
      FN = sum(FN),
      # Matthews Correlation Coefficient
      numerator   = (TP * TN) - (FP * FN),
      denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
      MCC = ifelse(denominator == 0, 0, numerator / denominator),
      # Sensitivity (recall)
      TPR = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
      # FDR: among calls, how many are false
      FDR = ifelse((TP + FP) > 0, FP / (TP + FP), NA_real_),
      # FPR: among true nulls, how many are falsely called
      FPR = ifelse(sum(!is_de) > 0, FP / sum(!is_de), NA_real_),
      .groups = "drop"
    )
  
  # Bind metrics to the simulation parameters for this i
  bind_cols(
    res,
    this %>% 
      # put sim_id and key scenario descriptors next to the metrics
      dplyr::select(
        sim_id,
        iters,
        batch_regime,
        effect_regime,
        scale_regime,
        prop_clusters_respond,
        batch_treat_pattern,
        assignment,
        link_cluster_to_treatment,
        libsize_cv,
        batch_lfc_sd,
        disp_mult_sd,
        lfc_treat_loc,
        lfc_treat_sd,
        de_prop_treat,
        de_prop_cluster,
        n_batches,
        patients_per_batch,
        cells_per_patient,
        ngenes
      )
  )
}

# -------------------------------------------------------------------
# 3. Run over all grid rows
# -------------------------------------------------------------------
whole_res <- map_dfr(seq_len(nrow(grid)), summarise_one_sim, grid = grid)
whole_res = whole_res %>% na.omit()

# Optional: save summary for later use
saveRDS(whole_res, "results/sim_summary.rds")

# Quick look
whole_res %>% 
  arrange(batch_regime, effect_regime, scale_regime, name) %>%
  head()


whole_res %>% 
  ggplot(aes(x = name, y = MCC, colour = name)) +
  facet_grid(
    prop_clusters_respond ~ effect_regime + batch_treat_pattern,
    labeller = label_both
  ) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Method",
    y = "Matthews Correlation Coefficient",
    title = "MCC across simulation regimes (cluster heterogeneity × batch–treatment patterns)"
  )

whole_res %>% 
  ggplot(aes(x = name, y = FDR, colour = name)) +
  facet_grid(
    prop_clusters_respond ~ effect_regime,
    labeller = label_both
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Method",
    y = "FDR",
    title = "FDR across regimes (target 5%)"
  )


whole_res$name %>% unique()
core_methods <- c("devil+sf:none+se:none", "devil+sf:psinorm+se:patient", "devil+batch+sf:psinorm+se:patient", "edgeR+batch")

whole_res %>% 
  filter(name %in% core_methods) %>%
  ggplot(aes(x = name, y = MCC, colour = name)) +
  facet_grid(
    prop_clusters_respond ~ batch_treat_pattern,
    labeller = label_both
  ) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "MCC",
    title = "DEVIL vs pseudo-bulk vs naïve across cluster-specific and confounded regimes"
  )

