
rm(list = ls())
library(patchwork)
library(ggplot2)

hm = readRDS("figures/hm_interaction.rds")
hm = ggplotify::as.ggplot(hm)

# Example of type I results

bar_venn_plot = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/venn_bar_plot.RDS")
tissue_plot = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/gsea_tissue_plot.RDS")
umaps = readRDS("results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/umaps.RDS")

bar_venn_plot = bar_venn_plot +
  guides(fill = guide_legend(nrow = 4))

umap_glm_private = umaps$`glmGamPoi private` +
  scale_color_manual(values = c("Old - Type I" = "goldenrod3", "Young - Type I" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

umap_glm_devil_shared = umaps$`glmGamPoi and devil` +
  scale_color_manual(values = c("Old - Type I" = "goldenrod3", "Young - Type I" = "#483D8B")) +
  theme(legend.position = "right", legend.title.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

bar_plot_gsea = readRDS("results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/age_type1/venn_bar_plot_gsea.RDS")

des = "
#AAAA#
#AAAA#
#AAAA#
BBBCCC
BBBCCC
DDEEFF
DDEEFF
"

des = "
AAAABBB
AAAABBB
AAAACCC
DEEEEEE
DEEEEEE
"

des = "
AAAABB
AAAABB
AAAACC
D#EEEE
D#EEEE
"

umaps = umap_glm_private + umap_glm_devil_shared + 
  plot_layout(ncol = 2, guides = "collect")


p = patchwork::free(hm) + 
  patchwork::free(tissue_plot) + patchwork::free(bar_plot_gsea) +
  #patchwork::free(bar_venn_plot) + patchwork::free(umap_glm_private) + patchwork::free(umap_glm_devil_shared) +
  patchwork::free(bar_venn_plot) + umaps +
  plot_layout(design = des) +
  plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 12), plot.tag = element_text(face = "bold"))
p
ggsave("figures/extended.pdf", plot = p, width = 16, height = 12, units = "in")
