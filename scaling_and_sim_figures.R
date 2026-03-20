
rm(list = ls())
require(tidyverse)
require(patchwork)
library(ggplot2)

dir.create("all_figures/scaling_and_sim/", recursive = T)

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

# MAIN ####
# Use MacaqueBrain and HSC
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/runtime.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pAB = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect")

pC = readRDS("de_analysis/nullpower/figures/RDS/main/qq_plot.rds") +
  MY_THEME +
  theme(legend.position = "right")

pD = readRDS("de_analysis/nullpower/figures/RDS/main/power_curve.rds") +
  MY_THEME +
  theme(legend.position = "none")

pCD = pC + pD + plot_layout(guides = "collect")

pE = readRDS("de_analysis/nullpower/figures/RDS/main/MCC_boxplot.rds") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  MY_THEME +
  theme(legend.position = "none")

# design = "
# AAAAA
# AAAAA
# AAAAA
# AAAAA
# ##CCC
# ##CCC
# ##DDD
# ##DDD
# EEEEE
# EEEEE
# EEEEE
# EEEEE"

design = "
AAA
AAA
#BB
#BB
#BB
DDD
DDD
"

library(patchwork)
final_plot = free(pAB) +
  free(pCD) +
  free(pE) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A", "", "C", "D", "E"))) &
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(filename = "all_figures/scaling_and_sim/main_2_v0.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 11.7, units = "in")
rm(pA, pAB, pC, pD, pE, pB, pCD, final_plot, design)

# EXTENDED ####
# Use MacaqueBrain and HSC
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/correlation.RDS") + theme(legend.position = "bottom")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/venn_plot.RDS")
pC = readRDS("timing_scaling/img/RDS/MacaqueBrain/large.RDS") +
  theme(legend.direction='horizontal',
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.box.margin = margin(0, 0, 0, 0)) +
  guides(color = guide_legend(ncol = 2))

pFG = readRDS("de_analysis/nullpower/figures/RDS/main/ecfd_ks_plot.rds") +
  MY_THEME +
  theme(legend.position = "bottom")
pH = readRDS("de_analysis/nullpower/figures/RDS/main/ptiming_ratio.rds")

design = "
AAAAABBB
AAAAABBB
AAAAACCC
AAAAACCC
AAAAACCC
FFFFFHHH
FFFFFHHH
FFFFFHHH"

final_plot = free(pA) + free(pB) + free(pC) +
  free(pFG) + free(pH) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
# final_plot
ggsave(filename = "all_figures/scaling_and_sim/ext_2.pdf", plot = final_plot, dpi = 600, width = 13.7, height = 12, units = "in")
rm(pAB, pC, pFG, pH, final_plot, design)

# SUPP SCALING ####
## Times and Memory ####
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/speedup.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory_ratio.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")

final_plot = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
final_plot
ggsave(filename = "all_figures/scaling_and_sim/supp_scaling_times_and_memory.pdf", plot = final_plot, dpi = 600, width = 11, height = 5, units = "in")
rm(pAB, pCD, final_plot, pA, pB, pC, pD)


## Small dataset ####
dataset_name = "baronPancreas"

path_to_rds = file.path("timing_scaling/img/RDS/", dataset_name)
pA = readRDS(file.path(path_to_rds, "runtime.RDS"))
pB = readRDS(file.path(path_to_rds, "speedup.RDS"))
pAB = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect")
pC = readRDS(file.path(path_to_rds, "memory.RDS"))
pD = readRDS(file.path(path_to_rds, "memory_ratio.RDS"))
pCD = (pC + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pD + theme(legend.position = "none")) +
  plot_layout(guides = "collect")
pEF = readRDS(file.path(path_to_rds, "correlation.RDS"))
pG = readRDS(file.path(path_to_rds, "venn_plot.RDS"))

design = "
AAAAAA
AAAAAA
CCCCCC
CCCCCC
EEEEGG
EEEEGG
EEEEGG"

supp_fig = free(pAB) + free(pCD) +
  free(pEF) + free(pG) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "D"))) &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(1, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
ggsave(filename = paste0("all_figures/scaling_and_sim/supp_scaling_",dataset_name,".pdf"), plot = supp_fig, dpi = 600, width = 12.7, height = 12.7, units = "in")

# # SUPP DE ANALYSIS ####
# ## All methods ####
# all_models_plots = readRDS("de_analysis/img/RDS/all_models.RDS")
# pA = all_models_plots$MCC + theme(legend.position = "bottom")
# pB = all_models_plots$timing + theme(legend.position = "bottom")
# pC = all_models_plots$failure_rate + theme(legend.position = "bottom")
#
# des = "
# AAAAA
# AAAAA
# BBBCC
# BBBCC"
#
# all_models_plot = free(pA) + free(pB) + free(pC) +
#   plot_layout(design = des) +
#   plot_annotation(tag_levels = c("A")) &
#   theme(plot.tag = element_text(face = 'bold'))
# ggsave("all_figures/scaling_and_sim/supp_sim_allmodels.pdf", all_models_plot, width = 10, height = 10, dpi = 600, units = "in")
# rm(all_models_plots, pA, pB, pC, all_models_plot, des)


## All datasets ####
author = "yazar"
for (author in c("hsc", "kumar", "yazar", "bca")) {

  qq20 = readRDS(file.path("de_analysis/nullpower/figures/RDS",author,"qq20.rds")) +
    facet_wrap(~is.pb) + ggtitle("20 patients")
  qq4 = readRDS(file.path("de_analysis/nullpower/figures/RDS",author,"qq4.rds")) +
    facet_wrap(~is.pb)+ ggtitle("4 patients")

  qq20@layers[[2]]$aes_params$linewidth = .8
  qq4@layers[[2]]$aes_params$linewidth = .8

  power_curve_20 = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/power_curve_20.rds")) +
    facet_wrap(~splitval, ncol = 2) + ggtitle("20 patients")
  power_curve_4 = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/power_curve_4.rds")) +
    facet_wrap(~splitval, ncol = 2) + ggtitle("4 patients")

  power_curve_20@layers[[1]]$aes_params$linewidth = .8
  power_curve_4@layers[[1]]$aes_params$linewidth = .8

  MMC_box = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/MCC_box.rds"))
  ks_plot = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/ecfd_ks_plot.rds")) + theme(legend.position = "right")
  ks_plot@layers[[1]]$aes_params$linewidth = .8

  timing_ratio = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/ptiming_ratio.rds"))
  timing = readRDS(file.path("de_analysis/nullpower/figures/RDS/",author,"/ptiming.rds"))

  # pA = readRDS(file.path("de_analysis/img/RDS/", author, "pvalues.RDS"))
  # pA = (pA$null_pvalue + theme(legend.direction='vertical', legend.position = "none", legend.box = "horizontal")) +
  #   (pA$de_pvalue + theme(legend.position = "none"))
  #
  # pBC = readRDS(file.path("de_analysis/img/RDS/", author, "MCCs.RDS"))
  # pB = pBC$cellwise + ggtitle("Cell-wise")
  # pC = pBC$patientwise + ggtitle("Patient-wise")
  # pBC = (pB + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  #   (pC + theme(legend.position = "none")) +
  #   plot_layout(guides = "collect")
  #
  # pDE = readRDS(file.path("de_analysis/img/RDS/", author, "ks_test.RDS"))
  # pD = pDE$cellwise + ggtitle("Cell-wise")
  # pE = pDE$patientwise + ggtitle("Patient-wise")
  # pDE = (pD + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  #   (pE + theme(legend.position = "none")) +
  #   plot_layout(guides = "collect")
  #
  # pF = readRDS(file.path("de_analysis/img/RDS/", author, "timing.RDS"))

  timing_ratio = timing_ratio + scale_x_discrete(limits = rev)
  timing = timing + scale_x_discrete(limits = rev)

  qq = qq4 + qq20 +
    plot_layout(guides = "collect")

  power = power_curve_4 + power_curve_20 +
    plot_layout(guides = "collect")

  MMC_box = MMC_box + scale_x_discrete(limits = rev)

  timings = timing + timing_ratio +
    plot_layout(guides = "collect")

  design = "
AAAA
BBBB
CCCC
DDDD
EEEE
  "

  final_plot = free(qq) / free(power) / free(MMC_box) / free(ks_plot) / free(timings) +
    plot_layout(design = design) +
    plot_annotation(tag_levels = c("A", "B", "C", "D", "E")) &
    theme(plot.tag = element_text(face = 'bold'))

  # final_plot = qq / power / MMC_box / ks_plot / timings +
  #   plot_layout(design = design, guides = "collect") +
  #   plot_annotation(tag_levels = c("A", "B", "C", "D", "E")) &
  #   theme(plot.tag = element_text(face = 'bold'))

  # final_plot = pA / pBC / pDE / pF +
  #   plot_layout(design = design) +
  #   plot_annotation(tag_levels = c("A", "B", "C", "D")) &
  #   theme(plot.tag = element_text(face = 'bold'))

  ggsave(paste0("all_figures/scaling_and_sim/supp_sim_",author,".png"), final_plot, width = 8 * 1.5, height = 11 * 1.5, dpi = 600, units = "in")
  ggsave(paste0("all_figures/scaling_and_sim/supp_sim_",author,".pdf"), final_plot, width = 8 * 1.5, height = 11 * 1.5, dpi = 600, units = "in")
}

# Clean supplementary for QQ plots and power
