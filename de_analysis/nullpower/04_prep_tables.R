
rm(list = ls())
require(tidyverse)

methods_main = c("edgeR (Pb)", "MAST (cell)", "limma (Pb)",
                 "glmGamPoi (cell)", "NEBULA", "devil (MOM-NB)")

methods_supp = c("edgeR (Pb)", "MAST (cell)", "Seurat (cell)", "limma (Pb)", "limma (cell)",
            "glmGamPoi (cell)", "NEBULA", "limmaDupCorr (cell)", "devil (MLE-NB)", "devil (MOM-NB)")

res <- readRDS("final_res/results.rds")

## Cell-wise ####
cw_table_main <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% methods_main) %>%
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

cw_table_supp <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% methods_supp) %>%
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
pw_table_main <- res %>%
  dplyr::filter(is.pb == TRUE) %>%
  dplyr::filter(name %in% methods_main) %>%
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

pw_table_supp <- res %>%
  dplyr::filter(is.pb == TRUE) %>%
  dplyr::filter(name %in% methods_supp) %>%
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
write.csv(cw_table_main, "final_res/cw_table_main.csv")
write.csv(pw_table_main, "final_res/pw_table_main.csv")
write.csv(cw_table_supp, "final_res/cw_table_supp.csv")
write.csv(pw_table_supp, "final_res/pw_table_supp.csv")


# Prep timing tables
folder_path = "timing_results/"
lf <- list.files(folder_path)
df <- lapply(1:length(lf), function(i) {
  author = strsplit(lf[i], "_")[[1]][1]
  readRDS(paste0(folder_path, lf[i]))   %>% 
    dplyr::filter(name %in% methods_supp) %>% 
    dplyr::mutate(is_main = name %in% methods_main) %>% 
    dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
    dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>% 
    dplyr::select(name, Time, , cell_order, is_main) %>% 
    dplyr::mutate(author = author)
}) %>% do.call("bind_rows", .)

df %>% 
  dplyr::group_by(name, cell_order, author) %>% 
  dplyr::summarise(m = median(Time), s = sd(Time)) %>% 
  dplyr::mutate(time = sprintf("%.3f ± %.3f", m, s)) %>% 
  #dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  #dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::select(name, cell_order, time, author) %>% 
  tidyr::pivot_wider(names_from = name, values_from = time) %>% 
  xtable::xtable()




