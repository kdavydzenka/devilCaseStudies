
rm(list = ls())
library(patchwork)

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

pA = readRDS("figures/plot_runtime_covariates.RDS") %>% ggplotify::as.ggplot()
pB = readRDS("figures/plot_speedup_covariates.RDS") %>% ggplotify::as.ggplot()
pC = readRDS("figures/plot_runtime_patients.RDS") %>% ggplotify::as.ggplot()
pD = readRDS("figures/plot_speedup_patients.RDS") %>% ggplotify::as.ggplot()

des = "AB\nCD"


p = patchwork::free(pA) + patchwork::free(pB) + patchwork::free(pC) + patchwork::free(pD) + 
  plot_layout(design = des, guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

saveRDS(p, "figures/plot.rds")
ggsave("figures/plot_scaling_cov_and_features.pdf", plot = p, width = 8, height = 6, units = "in")
