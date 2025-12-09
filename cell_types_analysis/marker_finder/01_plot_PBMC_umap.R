rm(list = ls())

library(readr)
library(tidyverse)
library(Seurat)
library(ggplot2)

Liu_datasets <- read_csv("data/Datasets/Liu/Liu_datasets.csv")
Labels       <- read_csv("data/Datasets/Liu/Labels.csv")

# ----- build count matrix -----
cnt_mat  = as.matrix(Liu_datasets)
cnt_mat  = cnt_mat[, 2:ncol(cnt_mat)]          # drop first column if it's an ID
gene_names = colnames(Liu_datasets)[-1]        # gene names are all but first column
cnt_mat  = t(cnt_mat)

# N = 1000
# gene_idxs = sample(1:nrow(cnt_mat), size = N)
gene_idxs = 1:nrow(cnt_mat)
cnt_mat = cnt_mat[gene_idxs, ]

cnt_mat <- as.matrix(
  apply(cnt_mat, 2, function(x) as.numeric(as.character(x)))
)

rownames(cnt_mat) = gene_names[gene_idxs]
cnt_mat = cnt_mat[rowMeans(cnt_mat) > .05, ]
gene_idxs = 1:nrow(cnt_mat)

# make sure column names of count matrix are character
colnames(cnt_mat) <- paste0("Cell", 1:ncol(cnt_mat))
Labels$cell_id = colnames(cnt_mat)

meta <- Labels %>% 
  filter(cell_id %in% colnames(cnt_mat)) %>%
  column_to_rownames("cell_id")

# reorder metadata to match columns of cnt_mat
meta <- meta[colnames(cnt_mat), , drop = FALSE]

# ----- create Seurat object and run UMAP -----
seu <- CreateSeuratObject(counts = cnt_mat, meta.data = meta)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 20)
seu <- RunUMAP(seu, dims = 1:10)

# ----- UMAP coloured by cell-type -----

my_large_palette <- c(
  "steelblue4",
  "#D98880",
  "goldenrod",
  "indianred3",
  "mediumpurple",
  "brown",
  "plum",
  "tan",
  "darkseagreen",
  'cyan3',
  "lightsteelblue",
  "peru",
  "olivedrab",
  "palevioletred",
  "firebrick"
)

seu$cell_type = seu$x

f = .2
idxs = sample(1:ncol(seu), as.integer(f * ncol(seu)))
seu_sub = seu[,idxs]

umap_plot = DimPlot(
  seu_sub,
  reduction = "umap",
  group.by  = "cell_type",
  pt.size   = 0.5
) + theme_bw() +
  scale_color_manual(values = my_large_palette) +
  ggtitle("") +
  labs(x = "UMAP 1", y = "UMAP 2")
umap_plot

umap_plot = ggplotify::as.ggplot(umap_plot)

saveRDS(umap_plot, "img/umap_plot.RDS")
