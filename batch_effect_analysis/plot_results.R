
rm(list = ls())
library(tidyverse)

grid = readRDS("data/param_grid.csv")

i = 21
whole_res = lapply(1:nrow(grid), function(i) {
  print(i)
  this = grid[i,]
  
  if (!file.exists(paste0("results/sim_",i,".rds"))) {
    print(paste0("Skipping ", i))
    return(dplyr::tibble())
  }
  
  # readRDS(paste0("results/sim_",i,".rds"))$method %>% unique()
  # readRDS(paste0("results/sim_",i,".rds")) %>% 
  #   dplyr::filter(method == "NEBULA+batch") %>% view()
  
  res = readRDS(paste0("results/sim_",i,".rds")) %>% 
    # dplyr::mutate(p_adj = ifelse(is.na(p_adj), 1, p_adj)) %>% 
    # dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc)) %>% 
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
  
  cbind(res, this)
})

whole_res = dplyr::bind_rows(whole_res)

whole_res %>% 
  ggplot(mapping = aes(x = name, y = MCC, col = name, linetype = assignment)) +
  facet_grid(scale_regime~effect_regime+batch_regime) +
  geom_boxplot() +
  theme(legend.position = "bottom")

whole_res %>% 
  ggplot(mapping = aes(x = name, y = FDR, col = name)) +
  facet_grid(scale_regime~effect_regime) +
  geom_boxplot() +
  theme(legend.position = "bottom")
  
