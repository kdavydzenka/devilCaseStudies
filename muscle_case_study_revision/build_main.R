
rm(list = ls())
library(patchwork)
library(ggplot2)

umap = readRDS("figures/umap_all.rds") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  coord_fixed()
umap_zoom = readRDS("figures/umap_zoom.rds") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  coord_fixed()
# hm = readRDS("figures/hm_interaction.rds")
# hm = ggplotify::as.ggplot(hm)
volcano = readRDS("figures/volcano_devil.RDS") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  facet_wrap(~test_type, nrow = 1)

go_plot = readRDS("figures/go_plot_devil.rds") +
  theme(legend.position = "bottom", legend.title.position = "top", legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2), size = guide_legend(nrow = 2))
semantic_plot = readRDS("figures/semantic_plot_overlay.rds") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

des = "
ABCCC
ABCCC
DDDDD
DDDDD
DDDDD
EEEEE
EEEEE
EEEEE
"

p = umap + umap_zoom + volcano + 
  patchwork::free(go_plot) + 
  patchwork::free(semantic_plot) +
  plot_layout(design = des) +
  plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 12), plot.tag = element_text(face = "bold"))

ggsave("figures/main.pdf", plot = p, width = 16, height = 18, units = "in")
