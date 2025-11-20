
rm(list = ls())
library(tidyverse)

param_grid = readRDS("data/param_grid.rds")
param_grid$i = 1:nrow(param_grid)

res = lapply(param_grid$i, function(idx) {
  print(idx)
  if (!file.exists(paste0("results/sim_", idx, ".rds"))) {
    print(paste0("Skipping ", idx))
    return(dplyr::tibble())
  }
  readRDS(paste0("results/sim_", idx, ".rds")) %>% 
    dplyr::mutate(idx = idx)
}) %>% do.call("bind_rows", .)

dir.create("summarized_results")
saveRDS(res, "summarized_results/results.rds")
