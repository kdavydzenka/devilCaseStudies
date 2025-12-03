
rm(list = ls())

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

pA = readRDS("figures/plot_runtime_covariates.RDS")
pB = readRDS("figures/plot_speedup_covariates.RDS")
pC = readRDS("figures/plot_runtime_patients.RDS")
pD = readRDS("figures/plot_speedup_patients.RDS")

des = "AB\nCD"

free(pA) + free(pB) + free(pC) + free(pD) + 
  plot_layout(design = des) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
