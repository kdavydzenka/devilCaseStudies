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

umap_data = dplyr::tibble(
  UMAP_1 = seu@reductions$umap@cell.embeddings[,1],
  UMAP_2 = seu@reductions$umap@cell.embeddings[,2],
  cell_type = seu$cell_type
)

saveRDS(umap_data, "results/umap_data.rds")
  
  

f = .5
umap_plot = umap_data %>% 
  dplyr::slice_sample(prop = f) %>% 
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2, col = cell_type)) +
  geom_point( size=.2) +
  scale_color_manual(values = my_large_palette) +
  theme_bw() +
  labs(x = "UMAP 1", y = "UMAP 2", col="Cluster") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme(text=element_text(size=12)) +
  labs(col = "")

saveRDS(umap_plot, "img/umap_plot.RDS")
