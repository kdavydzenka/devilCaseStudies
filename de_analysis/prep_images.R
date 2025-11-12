
rm(list = ls())
require(tidyverse)
require(patchwork)
source("utils_plots.R")
source("utils_supp_plots.R")

method_cellwise = c("glmGamPoi (cell)", "Devil (base)", "Nebula", "Seurat (Cell)", "MAST (cell)", "edgeR (Pb)", "limma (Pb)")
method_patientwise = c("glmGamPoi (cell)", "Devil (mixed)", "Nebula", "Seurat (Cell)", "MAST (cell)", "edgeR (Pb)", "limma (Pb)")
method_levels = c("devil", "DevilSF (mixed)", "Devil (base)","glmGamPoi", "glmGamPoi (cell)", "NEBULA", "Nebula", "Seurat (Cell)", "MAST (cell)", "edgeR (Pb)", "limma (Pb)")
#method_levels <- c("limma", "glmGamPoi", "glmGamPoi (cell)", "Nebula", "NEBULA", "Devil (mixed)", "Devil (base)", "Devil", "devil")
readRDS("nullpower/final_res/results.rds")$name %>% unique()

author = "kumar"
for (author in c("hsc", "kumar", "yazar", "bca")) {
  MCCs_boxplots = plot_MCCs_boxplots(author)
  pvalues_plot = plot_pvalues(author, method_cellwise, method_patientwise)
  MCCs_plot = plot_MCCs(author, method_cellwise, method_patientwise)
  ks_plots = plot_ks(author, method_cellwise, method_patientwise)
  timing_plot = plot_timing(author)

  dir.create(file.path("img/RDS/", author), recursive = TRUE, showWarnings = F)

  saveRDS(MCCs_boxplots, file.path("img/RDS/", author, "MCCs_boxplots.RDS"))
  saveRDS(pvalues_plot, file.path("img/RDS/", author, "pvalues.RDS"))
  saveRDS(MCCs_plot, file.path("img/RDS/", author, "MCCs.RDS"))
  saveRDS(ks_plots, file.path("img/RDS/", author, "ks_test.RDS"))
  saveRDS(timing_plot, file.path("img/RDS/", author, "timing.RDS"))
}

# All methods all together
all_models_plots = plot_all_models()
saveRDS(all_models_plots, file.path("img/RDS/", "all_models.RDS"))

des = "A\nB\nC"
all_models_plot = free(all_models_plots$MCC) + free(all_models_plots$timing) + free(all_models_plots$failure_rate) +
  plot_layout(design = des) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
ggsave("../figures/supp_03_a.pdf", all_models_plot, width = 10, height = 10, units = "in", dpi = 600)

# Supplemetary
cwise_pvalues_plot = get_cwise_pvalues_plot()
cwise_boxplots_plot = get_cwise_boxplots()
pwise_pvalues_plot = get_pwise_pvalues_plot()
pwise_boxplots_plot = get_pwise_boxplots()
timings_plot = get_timings_plot()


ggsave("../figures/supp_03_b.pdf", cwise_pvalues_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_c.pdf", cwise_boxplots_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_d.pdf", pwise_pvalues_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_e.pdf", pwise_boxplots_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_f.pdf", timings_plot, width = 14, height = 12, units = "in", dpi = 600)

# saveRDS(cwise_pvalues_plot, "img/RDS/all_cwise_pvalue.RDS")
# saveRDS(cwise_boxplots_plot, "img/RDS/all_cwise_MCC.RDS")
# saveRDS(pwise_pvalues_plot, "img/RDS/all_pwise_pvalue.RDS")
# saveRDS(pwise_boxplots_plot, "img/RDS/all_pwise_MCC.RDS")
# saveRDS(timings_plot, "img/RDS/timings_plot.RDS")
