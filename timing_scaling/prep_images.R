
rm(list = ls())
require(tidyverse)
library(patchwork)
source("utils_img.R")

model_levels = c("NEBULA - cpu", "glmGamPoi - cpu", "devil - cpu", "devil - a100", "devil - h100")
#model_levels = c("glmGamPoi - cpu", "devil - cpu", "devil - a100 (1)", "devil - h100 (1)", "devil - a100 (8)", "devil - h100 (8)")

dataset_name = "MacaqueBrain"
add_predicttions = TRUE
for (dataset_name in c("MacaqueBrain", "baronPancreas")) {
  results_folder <- paste0("results/", dataset_name)
  fits_folder = paste0("results/", dataset_name, "/fits/")
  
  time_results = get_time_results(results_folder)
  mem_results = get_memory_results(results_folder)
  
  if (add_predicttions) {
    time_results = plyr::rbind.fill(time_results, predict_time_results(time_results))
    mem_results = plyr::rbind.fill(mem_results, predict_memory_results(mem_results))
    time_results$Measure[is.na(time_results$Measure)] = "observed"
    mem_results$Measure[is.na(mem_results$Measure)] = "observed"
  } else {
    time_results$Measure = "observed"
    mem_results$Measure = "observed"
  }
  
  #time_results$model_name = ifelse(!is.na(time_results$n_gpus), paste0(time_results$model_name, " (", time_results$n_gpus ,")"), time_results$model_name)
  time_results$model_name = factor(time_results$model_name, levels = model_levels)
  mem_results$model_name = factor(mem_results$model_name, levels = model_levels)

  # Comparisons ####
  # pA <- time_comparison(time_results)
  # pB <- memory_comparison(mem_results)
  # pC <- time_comparison(time_results, ratio = "glmGamPoi - cpu")
  # pD <- memory_comparison(mem_results, ratio = "devil - cpu")
  
  time_comp_v2(time_results%>% dplyr::filter(n_cells == 1e06)) + memory_comp_v2(mem_results %>% dplyr::filter(n_cells == 1e06))
  
  pA <- time_comp_v2(time_results)
  pB <- memory_comp_v2(mem_results)
  pC <- time_comp_v2(time_results, ratio = "NEBULA - cpu") + scale_y_sqrt()
  pD <- memory_comp_v2(mem_results, ratio = "devil - cpu") + scale_y_sqrt()

  # Correlations ####
  corr_plots <- plot_correlations_by_model(fits_folder)
  # plot_correlations_single(fits_folder)
  # (fits_folder)

  # Plot UpSet ####
  # upset <- plot_upset(fits_folder, lfc_cut = 1, pval_cut = .05)
  # upset = ggplotify::as.ggplot(upset)
  
  venn_plot = plot_venn(fits_folder, lfc_cut = 1, pval_cut = .05)

  # Plot large test
  dir.create(file.path("img/RDS/", dataset_name), recursive = T, showWarnings = F)
  if (dataset_name == "MacaqueBrain") {
    pLarge = plot_large_test(dataset_name, add_competitors = T)
    pLarge = pLarge + scale_x_log10() + scale_y_sqrt() + geom_vline(xintercept = 1e6, linetype = "dotted", linewidth = .5, color = "darkslategray") +
      annotate("text", x = 1e6, y = 1000, label = "1 million cells", angle = 90, vjust = 1.4)
    saveRDS(pLarge, file.path("img/RDS/", dataset_name, "large.RDS"))
  }

  # Overdispersion
  # p_disp_runtime = plot_overdispersion_comparison(dataset_name, ratio = FALSE)
  # p_disp_speedup = plot_overdispersion_comparison(dataset_name, ratio = TRUE)

  saveRDS(pA, file.path("img/RDS/", dataset_name, "runtime.RDS"))
  saveRDS(pB, file.path("img/RDS/", dataset_name, "memory.RDS"))
  saveRDS(pC, file.path("img/RDS/", dataset_name, "speedup.RDS"))
  saveRDS(pD, file.path("img/RDS/", dataset_name, "memory_ratio.RDS"))
  saveRDS(corr_plots, file.path("img/RDS/", dataset_name, "correlation.RDS"))
  saveRDS(venn_plot, file.path("img/RDS/", dataset_name, "venn_plot.RDS"))
  #saveRDS(upset, file.path("img/RDS/", dataset_name, "upset.RDS"))
  # saveRDS(p_disp_runtime, file.path("img/RDS/", dataset_name, "disp_runtime.RDS"))
  # saveRDS(p_disp_speedup, file.path("img/RDS/", dataset_name, "disp_speedup.RDS"))
}
