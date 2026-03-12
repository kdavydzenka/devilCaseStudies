
rm(list = ls())
require(tidyverse)
require(patchwork)
library(ggplot2)

dir.create("all_figures/scaling_and_sim/", recursive = T)

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  #legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank(),
  text = element_text(size = 12), 
  plot.tag = element_text(face = "bold")
)

IMG_FOLDER = "de_analysis/nullpower/figures/RDS/devil_comparison/"

pA = readRDS("de_analysis/nullpower/figures/RDS/devil_comparison/MCC_box.rds") + theme(legend.position = "none")
pB = readRDS("de_analysis/nullpower/figures/RDS/devil_comparison/ptiming_ratio.rds") + theme(legend.position = "bottom")

des = "
AAAAA
AAAAA
AAAAA
#BBB#
#BBB#
"

final_plot = pA + pB + 
  plot_layout(design = des) + plot_annotation(tag_levels = "A") & MY_THEME
final_plot

ggsave(paste0("all_figures/scaling_and_sim/supp_devil_comparison.png"), final_plot, width = 8, height = 8, dpi = 600, units = "in")
ggsave(paste0("all_figures/scaling_and_sim/supp_devil_comparison.pdf"), final_plot, width = 8, height = 8, units = "in")

