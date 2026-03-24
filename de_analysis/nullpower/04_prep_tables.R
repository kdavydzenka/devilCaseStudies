
rm(list = ls())
require(tidyverse)
library(googlesheets4)
gs4_auth()

methods =  c("edgeR (Pb)", "MAST (cell)", "Seurat (cell)", "limma (Pb)", "limma (cell)",
             "glmGamPoi (cell)", "NEBULA", "limmaDupCorr (cell)", "devil (MLE-NB)", "devil (MOM-NB)")

res <- readRDS("final_res/results.rds")

## Cell-wise ####
cw_table <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% methods) %>%
  dplyr::group_by(name, author, patients, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(name, author, patients, MCC) %>% 
  dplyr::arrange(name) %>% 
  tidyr::pivot_wider(names_from = name, values_from = MCC)

## PatientWise ####
pw_table <- res %>%
  dplyr::filter(is.pb == TRUE) %>%
  dplyr::filter(name %in% methods) %>%
  dplyr::group_by(name, author, patients, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(name, author, patients, MCC) %>% 
  dplyr::arrange(name) %>% 
  tidyr::pivot_wider(names_from = name, values_from = MCC)

# Save tables ####

write.csv(cw_table, "final_res/cw_table.csv")
write.csv(pw_table, "final_res/pw_table.csv")

sheet_write(pw_table, ss = "https://docs.google.com/spreadsheets/d/1xRKVD-H5w2YOISFDRALWH7uYbvrLAgIfhZpRrGQ95gU", sheet = "Patient-wise MCC")
sheet_write(cw_table, ss = "https://docs.google.com/spreadsheets/d/1xRKVD-H5w2YOISFDRALWH7uYbvrLAgIfhZpRrGQ95gU", sheet = "Cell-wise MCC")

# Prep timing tables
time_df = res %>% 
  dplyr::filter(name %in% methods) %>%
  dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
  dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>% 
  dplyr::select(name, Time, cell_order, author) %>% 
  dplyr::group_by(name, cell_order, author) %>% 
  dplyr::summarise(
    across(
      c(Time), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  )

time_df = time_df %>% 
  tidyr::pivot_wider(values_from = Time, names_from = name)

write.csv(time_df, "final_res/timing.csv")
sheet_write(time_df, ss = "https://docs.google.com/spreadsheets/d/1xRKVD-H5w2YOISFDRALWH7uYbvrLAgIfhZpRrGQ95gU", sheet = "Patient and Cell-wise timing")
