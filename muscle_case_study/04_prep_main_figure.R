
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ComplexHeatmap", "magick")
sapply(pkgs, require, character.only = TRUE)
library(ggplot2)
library(patchwork)
library(grid)

cell_group_colors = c(
  "old" = "darkorange",
  "young" = "steelblue"
)

# input HEATMAP ####
source("utils.R")
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/multiomics/rna/seurat_counts_rna.RDS"
input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

# get markers 
de_res <- readRDS("results/MuscleRNA/devil_rna.RDS")
#de_res <- de_res %>% 
#  dplyr::filter(abs(lfc) >= 1, adj_pval <= .05)

gene_markers <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                  "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                  "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

de_res_top = de_res %>% 
  dplyr::filter(name %in% gene_markers)
  # dplyr::arrange(-abs(lfc)) %>% 
  # dplyr::slice_head(n = 50) %>% 
  # dplyr::pull(name)


N_subsample <- 10000
sample_idx = sample(1:ncol(input_data$counts), N_subsample, replace = FALSE)

rownames(input_data$counts)
mat <- input_data$counts[de_res_top,sample_idx] %>% as.matrix()
meta <- input_data$metadata[sample_idx,] 

meta = meta %>% 
  dplyr::mutate(Annotation = ifelse(Annotation == "Type II", 'Myonuclei TII', 'Myonuclei TI')) %>% 
  dplyr::mutate(age_pop = ifelse(age_pop == 'old_pop', "Old", "Young"))

# reorder by samples
ordered_indices = order(meta$sample)

mat <- mat[,ordered_indices]
meta = meta[ordered_indices,]

mat.scaled = t(apply(mat, 1, scale))

age_cluster = meta$age_pop
cell_type_cluster = meta$Annotation
patient_cluster = meta$sample

stats::quantile(mat.scaled, c(0.01, 0.95))

col_fun = circlize::colorRamp2(c(-.5, 0, .5), c("#FF00FF", "gray20", "#FFFF00"))

## Young pop ####
anno = as.data.frame(meta$age_pop)
colnames(anno) = "Age"
anno$Sample = meta$sample
anno$`Cell Type` = meta$Annotation

ha = HeatmapAnnotation(
  df = anno[,c(1:3)],
  annotation_height = unit(c(0.5,0.5, 0.5), "cm"),
  show_annotation_name = TRUE,
  annotation_name_offset = unit(2, "mm"),
  annotation_name_rot = c(0, 0, 0),
  col = list(
    Age = c("Old"="goldenrod3", "Young"="#483D8B"),
    `Cell Type` = c("Myonuclei TII"="#008080", "Myonuclei TI"="maroon"),
    Sample = c(
      "OM1" = '#1f77b4',
      'OM2' = '#aec7e8',
      'OM3' = '#ff7f0e',
      'OM5' = '#ffbb78',
      'OM6' = '#2ca02c',
      'OM7' = '#98df8a',
      'OM8' = '#d62728',
      'OM9' = '#ff9896',
      'P17' = '#9467bd',
      'P21' = '#c5b0d5',
      'P23' = '#8c564b',
      'P27' = '#c49c94',
      'P29' = '#e377c2',
      'P3' = '#f7b6d2',
      'P5' = '#7f7f7f',
      'YM1' = '#c7c7c7',
      'YM2' = '#bcbd22',
      'YM3' = '#dbdb8d',
      'YM4' = '#17becf'
    )
  )
)


hm <- Heatmap(
  mat.scaled,
  name = "Z-score", 
  km = 1,
  column_split = factor(meta$age_pop),
  #row_split = factor(de_res_top$lfc >= 0),
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = NULL,
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 5),
  column_gap = unit(0.5, "mm"),
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 8),
  column_title_rot = 90,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 10,
  top_annotation = ha
)
hm

saveRDS(hm, "plot/hm.rds")
saveRDS(mat.scaled, "plot/mat.scaled.rds")
saveRDS(meta, "plot/meta.rds")

pdf("plot/hm_complexHeatmp.pdf", width = 10, height = 6)
draw(hm)
dev.off()

hm <- ggplotify::as.ggplot(hm)
hm

ggsave("plot/hm.pdf", width = 10, height = 6, units = "in", dpi = 700)


# input UMAPs ####
load("results/data/metadata_rna_umap.Rdata")

# metadata_rna$cell_type %>% unique()
# 
# x_lims <- c(-5.5, 1)
# y_lims <- c(0, 11)
# 
# box = dplyr::tibble(
#   xmin=x_lims[1],
#   xmax=x_lims[2],
#   ymin=y_lims[1],
#   ymax=y_lims[2],
# )
# 
# 
# 
# umap_cell_types <- metadata_rna %>%
#   dplyr::sample_n(10001) %>% 
#   #`rownames<-`(1:nrow(metadata_atac)) %>%
#   dplyr::mutate(cell_type = if_else(cell_type %in% c("Myonuclei TI", "Myonuclei TII"), cell_type, "Other")) %>%
#   #dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
#   #dplyr::sample_n(5000) %>%
#   dplyr::select(umap_1, umap_2, cell_type, age_pop) %>%
#   dplyr::mutate(idx = row_number()) %>%
#   #dplyr::filter(umap_1 >= -1, umap_2 >= -6, umap_2 <= 7.5) %>%
#   ggplot(mapping = aes(x=umap_1, y=umap_2, col=cell_type)) +
#   geom_point(alpha = .5, size = .5) +
#   geom_rect(mapping = aes(xmin=x_lims[1], ymin=y_lims[1], xmax=x_lims[2], ymax=y_lims[2]), col="black", fill = alpha("white", 0)) +
#   scale_color_manual(values = c("Myonuclei TI" = "maroon", "Myonuclei TII"="#008080", "Other" = "gray90")) +
#   theme_bw() +
#   labs(x = "UMAP 1", y = "UMAP 2", col="Cell type")
# 
# umap_cell_types <- umap_cell_types +
#   theme_void() +
#   theme(legend.position = "none")
# 
# umap_ages <- metadata_rna %>%
#   #`rownames<-`(1:nrow(metadata_atac)) %>%
#   #dplyr::mutate(cell_type = if_else(cell_type %in% c("Myonuclei TI", "Myonuclei TII"), cell_type, "Other")) %>%
#   dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
#   dplyr::sample_n(10000) %>% 
#   dplyr::select(umap_1, umap_2, cell_type, age_pop) %>%
#   dplyr::mutate(idx = row_number()) %>%
#   #dplyr::filter(umap_1 >= -1, umap_2 >= -6, umap_2 <= 7.5) %>%
#   ggplot(mapping = aes(x=umap_1, y=umap_2, col=age_pop)) +
#   geom_point(alpha = .8, size = .5) +
#   theme_bw() +
#   scale_color_manual(values = c("old_pop" = "goldenrod3", "young_pop" = "#483D8B")) +
#   lims(x = x_lims) +
#   labs(x = "UMAP 1", y = "UMAP 2", col="Cell type") +
#   theme_void() +
#   theme(legend.position = "none")
# umap_ages
# 
# diff(x_lims)
# diff(y_lims)

# ggsave("plot/umap_cell_types.png", plot = umap_cell_types, width = 10, height = 10, units = "in", dpi = 1200)
# ggsave("plot/umap_ages.png", plot = umap_ages, width = diff(x_lims) / 2, height = diff(y_lims) / 2, units = "in", dpi = 1200)












# Define zoom area
df = metadata_rna %>%
  dplyr::sample_n(10000) %>% 
  dplyr::mutate(age_pop = ifelse(grepl("young", age_pop), "Young", "Old")) %>% 
  dplyr::mutate(cell_type = if_else(cell_type %in% c("Myonuclei TI", "Myonuclei TII"), cell_type, "Other")) %>%
  dplyr::select(umap_1, umap_2, cell_type, age_pop) %>%
  dplyr::mutate(idx = row_number())

x_lims <- c(-5.5, 1)
y_lims <- c(0, 11)
df_zoom <- df %>% 
  dplyr::filter(umap_1 >= x_lims[1], umap_1 <= x_lims[2], umap_2 >= y_lims[1], umap_2 <= y_lims[2])



# ----- Main plot with rectangle -----
p1 <- ggplot(df, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(
    name = "Cell Type",
    values = c(
      "Myonuclei TI" = "maroon",
      "Myonuclei TII" = "#008080",
      "Other" = "gray50"
    )
  ) +
  coord_fixed() +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))
  annotate("rect", xmin = x_lims[1], xmax = x_lims[2],
           ymin = y_lims[1], ymax = y_lims[2],
           color = "black", fill = NA, linewidth = 0.6, linetype = "dashed")

# ----- Zoomed-in plot (with panel border and no axes) -----
p2 <- ggplot(df_zoom, aes(x = umap_1, y = umap_2, color = age_pop)) +
  geom_point(alpha = 0.8, size = 0.5) +
  scale_color_manual(
    name = "Age Population",
    values = c(
      "Old" = "goldenrod3",
      "Young" = "#483D8B"
    )
  ) +
  coord_cartesian(xlim = x_lims, ylim = y_lims) +
  coord_fixed() +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom"  # Show this legend too
  ) +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

# ----- Patchwork combo -----
combined_plot <- p1 + 
  inset_element(p2, left = 0.6, bottom = 0.4, right = 1, top = 1, align_to = 'panel')

# Show both legends together below
final <- combined_plot + plot_layout(guides = "collect") & theme(legend.position = "right")

final

ggsave("plot/umap_both.png", plot = final, width = 8, height = 8, units = "in", dpi = 600)
saveRDS(final, "plot/umap_both.rds")
