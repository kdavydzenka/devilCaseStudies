
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


# Extract glm private genes
# Loading data #
rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

# Gene selection based on LFC & pvalue cutoff #
lfc_cut <- 1.0
pval_cut <- .05

rna_deg_devil <- rna_devil %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "devil")

rna_deg_glm <- rna_glm %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "glmGamPoi")

rna_deg_nebula <- rna_nebula %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "nebula")

glm_private = rna_deg_glm$geneID[!rna_deg_glm$geneID %in% rna_deg_nebula$geneID & !rna_deg_glm$geneID %in% rna_deg_devil$geneID]
glm_devil_shared = intersect(rna_deg_glm$geneID, rna_deg_devil$geneID)
all_shared = intersect(rna_deg_nebula$geneID, glm_devil_shared)

gene_list = list(
  "glmGamPoi private" = glm_private,
  "glmGamPoi and devil" = glm_devil_shared,
  "shared" = all_shared
)

N_subsample <- 10000
sample_idx = sample(1:ncol(input_data$counts), N_subsample, replace = FALSE)
for (i in 1:length(gene_list)) {
  list_name = names(gene_list)[i]
  intersting_genes = gene_list[[i]]
  
  mat <- input_data$counts[intersting_genes,sample_idx] %>% as.matrix()
  meta <- input_data$metadata[sample_idx,] 
  
  meta = meta %>% 
    dplyr::mutate(Annotation = ifelse(Annotation == "Type II", 'Myonuclei TII', 'Myonuclei TI')) %>% 
    dplyr::mutate(age_pop = ifelse(age_pop == 'old_pop', "Old", "Young"))
  
  # reorder by samples
  ordered_indices = order(meta$sample)
  
  mat <- mat[,ordered_indices]
  meta = meta[ordered_indices,]
  
  mat.scaled = t(apply(mat, 1, scale))
  
  # Check UMAP
  # 1. Create Seurat object
  # Ensure `mat` has gene names as rows and sample/cell names as columns
  seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta)
  
  # 2. Normalize and scale
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  
  # 3. PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # 4. UMAP
  seurat_obj = RunUMAP(seurat_obj, dims = 1:10)
  
  umap_df = as_tibble(seurat_obj@reductions$umap@cell.embeddings)
  umap_df$age = seurat_obj$age
  umap_df$age_pop = seurat_obj$age_pop
  umap_df$sample = seurat_obj$sample
  umap_df$Annotation = seurat_obj$Annotation
  
  saveRDS(object = umap_df, file = paste0("test_analysis/", list_name, ".RDS"))
  
  # 5. Plot UMAP — color by sample (or change to 'Age', 'Cell Type', etc.)
  p = DimPlot(seurat_obj, reduction = "umap", group.by = "age_pop", pt.size = .05) +
    ggtitle(list_name)
  #png(paste0("test_analysis/", list_name, ".png"), width = 8, height = 8, units="in", res=300)
  #print(p)
  #dev.off()
}
