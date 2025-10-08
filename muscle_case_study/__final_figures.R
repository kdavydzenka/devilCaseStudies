
rm(list = ls())
require(ggplot2)


# MAIN ####
pA = readRDS("plot/umap_both.rds")
pA = ggplotify::as.ggplot(pA)
pB = readRDS("plot/volcanos.rds")
pC = readRDS("plot/final_venn_plot.rds")
pD = readRDS("plot/final_GO_plot.rds")
pE = readRDS("plot/final_tissue_specific_dist_plot.rds") + coord_flip()

design <- "
AAABBBB
AAABBBB
DDEEEEE
DDEEEEE
FFFFFFF
FFFFFFF
FFFFFFF
FFFFFFF
"

final_plot <- free(pA) + free(pB) +
  free(pC) + free(pE) +
  free(pD) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

l = 15
x_ratio = 1
y_ratio = 1.2
ggsave("plot/final_main.png", width = x_ratio * l, height = y_ratio * l, units = "in", dpi = 300, plot = final_plot)


# EXTENDED ####

pA = readRDS("plot/hm.rds")
pA = ggplotify::as.ggplot(pA)
#pB = readRDS("plot/final_p_values_dist_plot.rds")
pC = readRDS("plot/final_umap_glm_private_plot.rds")
pD = readRDS("plot/final_umap_glm_devil_plot.rds")
pE = readRDS("plot/final_gene_group_impact_plot.rds")

design <- "
AAAAAA
AAAAAA
BBCCDD
BBCCDD
"

final_plot <- free(pA) + 
  free(pC) + 
  free(pD) +
  free(pE) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
l = 11
x_ratio = 1
y_ratio = .8
ggsave("plot/final_extended.png", width = x_ratio * l, height = y_ratio * l, units = "in", dpi = 300, plot = final_plot)
