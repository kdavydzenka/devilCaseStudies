
rm(list = ls())
library(patchwork)
library(ggrastr)
library(ggplot2)

umap = readRDS("figures/umap_all.rds") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  coord_fixed()
umap$layers[[1]]$aes_params$stroke <- 0
umap$layers[[1]]$aes_params$size <- .75
umap$layers[[1]]$aes_params$alpha <- .8
umap <- rasterise(umap, layers = "Point", dpi = 600)

# ggsave("~/Downloads/umap.pdf", plot = umap, width = 5, height = 5, units = "in")
# ggsave("~/Downloads/umap_rast.pdf", plot = umap_rast, width = 5, height = 5, units = "in")
# file.info(c("~/Downloads/umap.pdf", "~/Downloads/umap_rast.pdf"))$size

umap_zoom = readRDS("figures/umap_zoom.rds") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  coord_fixed()
umap_zoom$layers[[1]]$aes_params$stroke <- 0
umap_zoom$layers[[1]]$aes_params$size <- 1
umap_zoom$layers[[1]]$aes_params$alpha <- 0.8
umap_zoom <- rasterise(umap_zoom, layers = "Point", dpi = 600)
# hm = readRDS("figures/hm_interaction.rds")
# hm = ggplotify::as.ggplot(hm)
volcano = readRDS("figures/volcano_devil.RDS") +
  theme(legend.position = "bottom", legend.title.position = "top") +
  facet_wrap(~test_type, nrow = 1)
volcano$layers[[1]]$aes_params$stroke <- 0
volcano$layers[[1]]$aes_params$size <- 1.5
volcano$layers[[1]]$aes_params$alpha <- .8
volcano <- rasterise(volcano, layers = "Point", dpi = 600)

go_plot = readRDS("figures/go_plot_devil.rds") +
  theme(legend.position = "bottom", legend.title.position = "top", legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2), size = guide_legend(nrow = 2))
go_plot$layers[[1]]$aes_params$stroke <- 0

semantic_plot = readRDS("figures/semantic_plot_overlay.rds") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

semantic_plot$layers[[2]]$aes_params$stroke <- 0
semantic_plot$layers[[2]]$aes_params$alpha <- .8


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
  plot_annotation(tag_levels = "a") & 
  theme(text = element_text(size = 12), 
        plot.tag = element_text(face = "bold"))

p = umap + umap_zoom + volcano + 
  patchwork::free(go_plot) + 
  patchwork::free(semantic_plot) +
  plot_layout(design = des) +
  plot_annotation(tag_levels = "a") & 
  theme(text = element_text(size = 12), 
        plot.tag = element_text(face = "bold"))

ratio = 1.2
# ggsave("~/Downloads/Figure5_new.pdf", plot = p, width = 16, height = 18, units = "in")
ggsave("figures/main.pdf", plot = p, width = 16, height = 18, units = "in")
#ggsave("figures/main.pdf", plot = p, width = 16, height = 18, units = "in")
