
rm(list = ls())
library(tidyverse)
library(patchwork)

dir.create("figures/supplementary_figures", recursive = T)

# Volcanos 
volcano = readRDS("figures/volcano_all.RDS")
ggsave("figures/supplementary_figures/volcanos_all.png", width = 13, height = 10, plot = volcano, dpi = 400)

# GSEA sub v full
bar_plot_gsea_typeI = readRDS("results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/age_type1/venn_bar_plot_gsea.RDS")
bar_plot_gsea_typeII = readRDS("results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/age_type2/venn_bar_plot_gsea.RDS")
bar_plot_gsea = bar_plot_gsea_typeI + bar_plot_gsea_typeII + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/supplementary_figures/sub_bar_plot_gsea.pdf", width = 13, height = 5, plot = bar_plot_gsea)

# UMAPS
bar_venn_plot_type1 = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/venn_bar_plot.RDS")
umaps_type1 = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/umaps.RDS")
bar_venn_plot_type1 = bar_venn_plot_type1 +
  guides(fill = guide_legend(nrow = 4))

umap_glm_private_1 = umaps_type1$`glmGamPoi private` +
  scale_color_manual(values = c("Old - Type I" = "goldenrod3", "Young - Type I" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

umap_glm_devil_shared_1 = umaps_type1$`glmGamPoi and devil` +
  scale_color_manual(values = c("Old - Type I" = "goldenrod3", "Young - Type I" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))


bar_venn_plot_type2 = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type2/venn_bar_plot.RDS")
umaps_type2 = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type2/umaps.RDS")
bar_venn_plot_type2 = bar_venn_plot_type2 +
  guides(fill = guide_legend(nrow = 4))

umap_glm_private_2 = umaps_type2$`glmGamPoi private` +
  scale_color_manual(values = c("Old - Type II" = "goldenrod3", "Young - Type II" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

umap_glm_devil_shared_2 = umaps_type2$`glmGamPoi and devil` +
  scale_color_manual(values = c("Old - Type II" = "goldenrod3", "Young - Type II" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

des = "
AA#BBBBCCCC
DD#EEEEFFFF"


p = patchwork::free(bar_venn_plot_type1) + patchwork::free(umap_glm_private_1) + patchwork::free(umap_glm_devil_shared_1) +
  patchwork::free(bar_venn_plot_type2) + patchwork::free(umap_glm_private_2) + patchwork::free(umap_glm_devil_shared_2) + 
  plot_layout(design = des) +
  plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(face = "bold")) & theme(legend.position = "bottom")

ggsave("figures/supplementary_figures/umaps_devil_v_glmGamPoi.png", width = 10, height = 10, plot = p, dpi = 400)


# GO plots
go_plot_glm = readRDS("figures/go_plot_glmGamPoi.rds") +
  theme(legend.position = "bottom", legend.title.position = "top", legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2), size = guide_legend(nrow = 2))

go_plot_nebula = readRDS("figures/go_plot_NEBULA.rds") +
  theme(legend.position = "bottom", legend.title.position = "top", legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2), size = guide_legend(nrow = 2))

go_plots = patchwork::free(go_plot_glm) + patchwork::free(go_plot_nebula) +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/supplementary_figures/go_plots_glm_and_nebula.pdf", width = 10, height = 12, plot = go_plots)

# Semantic plot
semantic_plot = readRDS("figures/semantic_plot_grid.rds") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

ggsave("figures/supplementary_figures/semantic_plot_grid.pdf", width = 13, height = 15, plot = semantic_plot)
