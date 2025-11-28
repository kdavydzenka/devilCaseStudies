
rm(list = ls())
library(tidyverse)
library(patchwork)

NAME_MAPPING = c(
  "devil+batch+sf:none+se:patient" = "devil (Sandwich+Batch)",
  "devil+sf:none+se:patient" = "devil (Sandwich)",
  "devil+sf:none+se:none" = "devil (plain)"
)

PALETTE = c(
  "devil (Sandwich+Batch)" = "#6e016b",
  "devil (Sandwich)" = "#8c6bb1", 
  "devil (plain)" = "#9ebcda"
)

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
whole_res$name = NAME_MAPPING[whole_res$name]
whole_res$name = factor(whole_res$name, levels = NAME_MAPPING)

# Optional: save summary for later use
saveRDS(whole_res, "results/sim_summary.rds")
whole_res = readRDS("results/sim_summary.rds")

plot_evaluation = function(scale_regimes) {
  #if (!is.null(batch_regimes)) whole_res = whole_res %>% dplyr::filter(batch_regime %in% batch_regimes)
  if (!is.null(scale_regimes)) whole_res = whole_res %>% dplyr::filter(scale_regime %in% scale_regimes)
  #if (!is.null(effect_regimes)) whole_res = whole_res %>% dplyr::filter(effect_regime %in% effect_regimes)
  
  mcc_plot = whole_res %>% 
    ggplot(aes(x = name, y = MCC, colour = name)) +
    ggh4x::facet_nested(paste0("Scale (", scale_regime, ")") ~ paste0("Batch noise (",batch_regime, ")") + paste0("DE signal (", effect_regime, ")")) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(), axis.ticks.x = element_blank()
    ) +
    labs(
      x = "Method",
      y = "MCC"
    ) +
    scale_color_manual(values = PALETTE) +
    labs(col = "Method")
  
  fdr_plot = whole_res %>% 
    ggplot(aes(x = name, y = FDR, colour = name))  +
    ggh4x::facet_nested(paste0("Scale (", scale_regime, ")") ~ paste0("Batch noise (",batch_regime, ")") + paste0("DE signal (", effect_regime, ")")) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(), axis.ticks.x = element_blank()
    ) +
    labs(
      x = "Method",
      y = "FDR"
    ) +
    scale_color_manual(values = PALETTE) +
    labs(col = "Method")
  
  # Load Umaps
  umap_names = list.files("results/umaps/")
  umap_names = umap_names[grepl(scale_regimes,umap_names)]
  
  
  n = umap_names[1]
  umaps = lapply(umap_names, function(n){
    
    u = readRDS(file.path("results/umaps/", n))
    
    if (scale_regimes == "small") u@layers[[1]]$aes_params$size = 0.8
    if (scale_regimes == "large") u@layers[[1]]$aes_params$size = 0.6
  
    u +
      labs(col = "DE Group", shape = "Batch") +
      scale_color_manual(values = c("firebrick3", "darkblue")) +
      scale_shape_manual(values = c(0, 1, 2, 5, 6)) +
      guides(
        colour = guide_legend(override.aes = list(size = 2)),
        shape  = guide_legend(override.aes = list(size = 2))
      )
  })
  names(umaps) = stringr::str_replace_all(umap_names, ".RDS", "")
  
  des = "ABCD\nEEEE\nFFFF"
  umaps[[1]] + umaps[[2]] + umaps[[3]] + umaps[[4]] + 
    mcc_plot + fdr_plot + 
    plot_annotation(tag_levels = list(c("A","","","","B","C"))) +
    plot_layout(design = des,  guides = "collect") & 
    theme(legend.direction = "vertical", plot.tag = element_text(face = 'bold'))
}

p_small = plot_evaluation(scale_regimes = "small")
p_large = plot_evaluation(scale_regimes = "large")

ggsave("figures/bench_batch_small.pdf", plot = p_small, width = 10, height = 8, units = "in")
ggsave("figures/bench_batch_large.pdf", plot = p_large, width = 10, height = 8, units = "in")
