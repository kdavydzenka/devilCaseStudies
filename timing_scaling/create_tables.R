
rm(list = ls())
require(tidyverse)
library(googlesheets4)
library(patchwork)
source("utils_img.R")
gs4_auth()

# MacaqueBrain ####
results_folder <- "results/MacaqueBrain/"
fits_folder = "results/MacaqueBrain/fits/"
time_results = get_time_results(results_folder) %>%
  dplyr::mutate(measure = "time") %>% dplyr::rename(value = time)
mem_results = get_memory_results(results_folder) %>%
  dplyr::mutate(measure = "memory") %>%
  dplyr::mutate(memory = as.numeric(memory)) %>%
  dplyr::rename(value = memory)
results = dplyr::bind_rows(time_results, mem_results)

table <- results %>%
  dplyr::group_by(n_genes, n_cells, model_name, measure) %>%
  dplyr::summarise(mean = mean(value) / 1) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(values_from = mean, names_from = measure) %>%
  dplyr::group_by(n_genes, n_cells) %>%
  dplyr::mutate(time_ratio = time / time[model_name == 'devil - a100']) %>%
  dplyr::mutate(memory_ratio = memory / memory[model_name == 'devil - a100']) %>%
  dplyr::select(model_name, n_genes, time, time_ratio, memory, memory_ratio) %>%
  tidyr::pivot_wider(names_from = model_name, values_from = c(time, time_ratio, memory, memory_ratio))
table <- table[,!grepl("ratio_devil (GPU)", colnames(table), fixed = T)]
table <- table[,!grepl("ratio", colnames(table), fixed = T)]
write.csv(table, "results/macaque_table.csv")

sheet_write(table, ss = "https://docs.google.com/spreadsheets/d/1xRKVD-H5w2YOISFDRALWH7uYbvrLAgIfhZpRrGQ95gU", sheet = "scaling_Chiou_et_al")

# xtable::xtable(table, caption = "Macauqe Brain", label = "tab:MacaqueBrain", digits = 2)

# baronPancreas ####
results_folder <- "results/baronPancreas/"
fits_folder = "results/baronPancreas/fits/"
time_results = get_time_results(results_folder) %>%
  dplyr::mutate(measure = "time") %>% dplyr::rename(value = time)
mem_results = get_memory_results(results_folder) %>%
  dplyr::mutate(measure = "memory") %>%
  dplyr::mutate(memory = as.numeric(memory)) %>%
  dplyr::rename(value = memory)
results = dplyr::bind_rows(time_results, mem_results)

table <- results %>%
  dplyr::group_by(n_genes, n_cells, model_name, measure) %>%
  dplyr::summarise(mean = mean(value) / 1) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(values_from = mean, names_from = measure) %>%
  dplyr::group_by(n_genes, n_cells) %>%
  dplyr::mutate(time_ratio = time / time[model_name == 'devil - a100']) %>%
  dplyr::mutate(memory_ratio = memory / memory[model_name == 'devil - a100']) %>%
  dplyr::select(model_name, n_genes, time, time_ratio, memory, memory_ratio) %>%
  tidyr::pivot_wider(names_from = model_name, values_from = c(time, time_ratio, memory, memory_ratio))
table <- table[,!grepl("ratio_devil (GPU)", colnames(table), fixed = T)]
table <- table[,!grepl("ratio", colnames(table), fixed = T)]

write.csv(table, "results/pancreas_table.csv")
sheet_write(table, ss = "https://docs.google.com/spreadsheets/d/1xRKVD-H5w2YOISFDRALWH7uYbvrLAgIfhZpRrGQ95gU", sheet = "scaling_Baron_et_al")
