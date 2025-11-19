
rm(list = ls())
library(tidyverse)

MY_PALETTE = c(
  "Devil (base)" = "#099668",
  "Devil (mixed)" = "#099668",
  "Devil" = "#099668",
  "devil" = "#099668",
  "Nebula" = "steelblue",
  "NEBULA" = "steelblue",
  "edgeR" = "#7D629E",
  "edgeR (Pb)" = "#7D629E",
  "limma" = "#B96461",
  "limma (Pb)" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "glmGamPoi" = "#EAB578",
  "limmaDupCorr (cell)" = "#8B0000",
  "limmaDupCorr" = "#8B0000",
  "Seurat (cell)" = "#708090",
  "Seurat (Wilcox)" = "#708090",
  "Seurat" = "#708090",
  "MAST (cell)" = "#D8BFD8",
  "MAST" = "#D8BFD8"
)

res = readRDS("summarized_results/results.rds")

res %>% 
  na.omit() %>% 
  tidyr::pivot_longer(!c(ncells, n_patients, n_covariates, n_genes, seed, continuous, idx)) %>% 
  ggplot(mapping = aes(x = n_patients, y = value, col = name, alpha = as.factor(n_covariates))) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  ggh4x::facet_nested(~"N covariates"+n_covariates) +
  scale_color_manual(values = MY_PALETTE) +
  theme_bw()

res %>% 
  na.omit() %>% 
  tidyr::pivot_longer(!c(ncells, n_patients, n_covariates, n_genes, seed, continuous)) %>% 
  dplyr::filter(name != "idx") %>% 
  ggplot(mapping = aes(x = n_covariates, y = value, col = name)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  ggh4x::facet_nested(~"N patients"+n_patients) +
  scale_color_manual(values = MY_PALETTE) +
  theme_bw()


BASE = "NEBULA"
res %>% 
  na.omit() %>% 
  tidyr::pivot_longer(!c(ncells, n_patients, n_covariates, n_genes, seed, continuous, idx)) %>% 
  dplyr::group_by(idx) %>% 
  dplyr::mutate(value = value[name == BASE] / value) %>% 
  ggplot(mapping = aes(x = n_patients, y = value, col = name, linetype = as.factor(n_covariates))) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  #ggh4x::facet_nested(~"N covariates"+n_covariates) +
  scale_color_manual(values = MY_PALETTE) +
  theme_bw() +
  labs(x = "N patients", y = "Speedup (w.r.t NEBULA)")

res %>% 
  na.omit() %>% 
  tidyr::pivot_longer(!c(ncells, n_patients, n_covariates, n_genes, seed, continuous)) %>% 
  ggplot(mapping = aes(x = n_covariates, y = value, col = name)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  ggh4x::facet_nested(~"N patients"+n_patients) +
  scale_color_manual(values = MY_PALETTE) +
  theme_bw()


