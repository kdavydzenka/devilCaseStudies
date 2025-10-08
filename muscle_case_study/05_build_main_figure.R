
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ComplexHeatmap")
sapply(pkgs, require, character.only = TRUE)

venn <- readRDS("plot/venn_plot_v2.rds")
volcanos <- readRDS("plot/volcanos.rds")
go_plots <- readRDS("plot/enrichment_dotplot.RDS")
simp_plot <- readRDS("plot/simp_plot.rds")

design <- "
BBBC
BBBC
#DDD
#DDD
#DDD"

design <- "
BBBBBBCC
BBBBBBCC
DDEEEEEE
DDEEEEEE
DDEEEEEE
DDEEEEEE
DDEEEEEE"

design <- "
BBEEE
BBEEE
BBEEE
CDEEE
"

final_plot <- free(volcanos) +
  free(venn) +
  free(simp_plot + theme(legend.position = "bottom")) +
  free(go_plots) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("C", 'D', "E", "F"))) &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.position = "bottom"
  )
final_plot

x = 12
ggsave("plot/main_figure_bottom.pdf", width = 1.43 * x, height = x, dpi = 600, units = "in", plot = final_plot)
