
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
source("utils/utils.R")
dataset_name <- "MuscleRNA"
data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_counts_rna.RDS"

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)

# get markers 
de_res <- readRDS("results/MuscleRNA/full/devil_age_only.RDS")

gene_markers = de_res %>% 
  dplyr::filter(adj_pval <= .05) %>% 
  dplyr::arrange(-abs(lfc)) %>% 
  dplyr::slice_head(n = 100) %>% 
  dplyr::pull(name)
  
gene_markers_paper <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                  "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                  "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

gene_markers = c(gene_markers, gene_markers_paper)

de_res_top = de_res %>% 
  dplyr::filter(name %in% gene_markers)
de_res_top %>% ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point()

set.seed(1234)
N_subsample <- 10000
sample_idx = sample(1:ncol(input_data$counts), N_subsample, replace = FALSE)

gene_markers = gene_markers[gene_markers %in% rownames(input_data$counts)]
mat <- input_data$counts[gene_markers,sample_idx] %>% as.matrix()
meta <- input_data$metadata[sample_idx,] 

dim(input_data$counts)

input_data$metadata$sample %>% unique() %>% length()


dim(devil_interaction$input_matrix)

meta = meta %>% 
  dplyr::mutate(`Cell type` = ifelse(cell_type == "Type II", 'Myonuclei TII', 'Myonuclei TI')) %>% 
  dplyr::mutate("Age" = ifelse(age_pop == 'old_pop', "Old", "Young"))

# reorder by samples
ordered_indices = order(meta$sample)

mat <- mat[,ordered_indices]
meta = meta[ordered_indices,]

#mat.scaled = t(apply(mat, 1, scale))

mat.log <- log1p(mat)
mat.scaled <- t(scale(t(mat.log)))


age_cluster = meta$Age
cell_type_cluster = meta$`Cell type`
patient_cluster = meta$sample

q <- quantile(mat.scaled, c(0.05, 0.5, 0.95))  # 2%, 50%, 98%
lim <- max(abs(q[1]), abs(q[3]))

hm_colors = c("#2166AC", "gray95", "#B2182B")
#hm_colors = c("#FF00FF", "gray20", "#FFFF00")
#hm_colors = c("firebrick", "gray20", "forestgreen")

col_fun <- circlize::colorRamp2(
  breaks = c(-lim, 0, lim),
  colors = hm_colors
)


#col_fun = circlize::colorRamp2(c(-.5, 0, 1.7), c("steelblue", "white", "firebrick"))

## Young pop ####
anno = as.data.frame(meta$Age)
colnames(anno) = "Age"
anno$Sample = meta$sample
anno$`Cell Type` = meta$`Cell type`

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


# Build row labels: show name only if it's in the paper list
row_lab <- rownames(mat.scaled)
row_lab[!row_lab %in% gene_markers_paper] <- ""

genes_to_mark <- gene_markers_paper
mark_indices <- which(row_lab %in% genes_to_mark)

row_ha = rowAnnotation(
  mark = anno_mark(
    at = mark_indices,
    labels = row_lab[mark_indices],
    labels_gp = gpar(fontsize = 8)
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
  row_title = NULL,
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 5),
  column_gap = unit(0.5, "mm"),
  show_row_dend = FALSE,
  col = col_fun,
  cluster_rows = T,
  row_labels = row_lab,
  row_names_gp = gpar(fontsize = 5),
  column_title_rot = 90,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 10,
  top_annotation = ha, 
  show_row_names = FALSE,
  right_annotation = row_ha
)

hm

saveRDS(hm, "figures/hm.rds")
saveRDS(mat.scaled, "figures/mat.scaled.rds")
saveRDS(meta, "figures/meta.rds")

pdf("figures/hm_complexHeatmp.pdf", width = 8, height = 6)
draw(hm)
dev.off()

hm <- ggplotify::as.ggplot(hm)
hm

ggsave("figures/hm.pdf", width = 8, height = 6, units = "in", dpi = 700)

# Heatmap interaction ####
gene_markers_type1 = readRDS("results/MuscleRNA/full/devil_age_type1.RDS") %>% 
  dplyr::filter(adj_pval <= .05) %>% 
  dplyr::arrange(-abs(lfc)) %>% 
  dplyr::slice_head(n = 100) %>% 
  dplyr::pull(name)

gene_markers_type2 = readRDS("results/MuscleRNA/full/devil_age_type2.RDS") %>% 
  dplyr::filter(adj_pval <= .05) %>% 
  dplyr::arrange(-abs(lfc)) %>% 
  dplyr::slice_head(n = 100) %>% 
  dplyr::pull(name)

gene_markers_paper <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                        "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                        "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

gene_markers = c(gene_markers_type1, gene_markers_type2, gene_markers_paper)

N_subsample <- 10000
sample_idx = sample(1:ncol(input_data$counts), N_subsample, replace = FALSE)

gene_markers = gene_markers[gene_markers %in% rownames(input_data$counts)]
mat <- input_data$counts[gene_markers,sample_idx] %>% as.matrix()
meta <- input_data$metadata[sample_idx,] 

meta = meta %>% 
  dplyr::mutate(`Cell type` = ifelse(cell_type == "Type II", 'Myonuclei TII', 'Myonuclei TI')) %>% 
  dplyr::mutate("Age" = ifelse(age_pop == 'old_pop', "Old", "Young"))

# reorder by samples
ordered_indices = order(meta$sample)

mat <- mat[,ordered_indices]
meta = meta[ordered_indices,]

#mat.scaled = t(apply(mat, 1, scale))

mat.log <- log1p(mat)
mat.scaled <- t(scale(t(mat.log)))


age_cluster = meta$Age
cell_type_cluster = meta$`Cell type`
patient_cluster = meta$sample

q <- quantile(mat.scaled, c(0.05, 0.5, 0.95))  # 2%, 50%, 98%
lim <- max(abs(q[1]), abs(q[3]))

#hm_colors = c("#2166AC", "gray95", "#B2182B")
#hm_colors = c("#FF00FF", "gray20", "#FFFF00")
#hm_colors = c("firebrick", "gray20", "forestgreen")

col_fun <- circlize::colorRamp2(
  breaks = c(-lim, 0, lim),
  colors = hm_colors
)


#col_fun = circlize::colorRamp2(c(-.5, 0, 1.7), c("steelblue", "white", "firebrick"))

## Young pop ####
anno = as.data.frame(meta$Age)
colnames(anno) = "Age"
anno$Sample = meta$sample
anno$`Cell Type` = meta$`Cell type`

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


# Build row labels: show name only if it's in the paper list
row_lab <- rownames(mat.scaled)
row_lab[!row_lab %in% gene_markers_paper] <- ""

genes_to_mark <- gene_markers_paper
mark_indices <- which(row_lab %in% genes_to_mark)

row_ha = rowAnnotation(
  mark = anno_mark(
    at = mark_indices,
    labels = row_lab[mark_indices],
    labels_gp = gpar(fontsize = 8)
  )
)

hm <- Heatmap(
  mat.scaled,
  name = "Z-score", 
  km = 1,
  #column_split = factor(meta$age_pop),
  column_split = factor(paste0(meta$age_pop, " - ", meta$cell_type)),
  #row_split = factor(de_res_top$lfc >= 0),
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = NULL,
  row_title = NULL,
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 5),
  column_gap = unit(0.5, "mm"),
  show_row_dend = FALSE,
  col = col_fun,
  cluster_rows = T,
  row_labels = row_lab,
  row_names_gp = gpar(fontsize = 5),
  column_title_rot = 90,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 10,
  top_annotation = ha, 
  show_row_names = FALSE,
  right_annotation = row_ha
)

hm

saveRDS(hm, "figures/hm_interaction.rds")
saveRDS(mat.scaled, "figures/mat.scaled_interaction.rds")
saveRDS(meta, "figures/meta_interaction.rds")

pdf("figures/hm_complexHeatmp_interaction.pdf", width = 8, height = 6)
draw(hm)
dev.off()

hm <- ggplotify::as.ggplot(hm)
hm

ggsave("figures/hm_interaction.pdf", width = 8, height = 6, units = "in", dpi = 700)

# UMAPs ####
load("results/metadata_rna_umap.Rdata")

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

# Main plot
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
  # coord_fixed() +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

# Zoomed-in plot
p2 <- ggplot(df_zoom, aes(x = umap_1, y = umap_2, color = age_pop)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(
    name = "Age",
    values = c(
      "Old" = "goldenrod3",
      "Young" = "#483D8B"
    )
  ) +
  coord_cartesian(xlim = x_lims, ylim = y_lims) +
  # coord_fixed() +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))

p1
p2

# Patchwork combo

# Show both legends together below
# final <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "right")
# final
# 
# ggsave("figures//umap_both.png", plot = final, width = 8, height = 8, units = "in", dpi = 600)
saveRDS(p1, "figures/umap_all.rds")
saveRDS(p2, "figures/umap_zoom.rds")
