
rm(list = ls())
library(patchwork)
source("../figure_theme.R")

pA = readRDS("figures/plot_runtime_covariates.RDS") %>% ggplotify::as.ggplot()
pB = readRDS("figures/plot_speedup_covariates.RDS") %>% ggplotify::as.ggplot()
pC = readRDS("figures/plot_runtime_patients.RDS") %>% ggplotify::as.ggplot()
pD = readRDS("figures/plot_speedup_patients.RDS") %>% ggplotify::as.ggplot()

des = "AB\nCD"


p = patchwork::free(pA) + patchwork::free(pB) + patchwork::free(pC) + patchwork::free(pD) +
  plot_layout(design = des, guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  MY_THEME

saveRDS(p, "figures/plot.rds")
ggsave("figures/plot_scaling_cov_and_features.pdf", plot = p, width = 10, height = 7, units = "in")
