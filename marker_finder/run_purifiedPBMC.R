
library(readr)
library(tidyverse)

Liu_datasets <- read_csv("data/Datasets/Liu/Liu_datasets.csv")
Labels <- read_csv("data/Datasets/Liu/Labels.csv")

cnt_mat = as.matrix(Liu_datasets)
cnt_mat = cnt_mat[,2:ncol(cnt_mat)]
gene_names = colnames(Liu_datasets)
cnt_mat = t(cnt_mat)

# N = 1000
# gene_idxs = sample(1:nrow(cnt_mat), size = N)
gene_idxs = 1:nrow(cnt_mat)
cnt_mat = cnt_mat[gene_idxs,]

cnt_mat <- as.matrix(
  apply(cnt_mat, 2, function(x) as.numeric(as.character(x)))
)

rownames(cnt_mat) = gene_names[gene_idxs]
cnt_mat = cnt_mat[rowMeans(cnt_mat) > .05,]
gene_idxs = 1:nrow(cnt_mat)

# UTILS

devil.fit = function(ct, cnt_mat, design_matrix, contrast, sf) {
  fit = devil::fit_devil(cnt_mat, design_matrix, overdispersion = T, size_factors = sf, 
                         max_iter = 500, parallel.cores = 1, verbose = T)
  test.res = devil::test_de(fit, c(0,1), max_lfc = 20, rep(1, nrow(design_matrix)))
  test.res %>% colnames()
  colnames(test.res) = c("gene", "pval", "padj", "lfc")
  test.res %>% dplyr::mutate(cell_type = ct)
}

glm.fit = function(ct, cnt_mat, design_matrix, contrast) {
  fit = glmGamPoi::glm_gp(cnt_mat, design = design_matrix)
  test.res = glmGamPoi::test_de(fit, contrast = contrast)
  test.res = test.res %>% dplyr::select(name, pval, adj_pval, lfc)
  colnames(test.res) = c("gene", "pval", "padj", "lfc")
  test.res %>% dplyr::mutate(cell_type = ct)
}

nebula.fit = function(ct, cnt_mat, design_matrix, contrast) {
  fit = nebula::nebula(cnt_mat, id = rep(1, nrow(design_matrix)), 
                       pred = design_matrix, mincp = 0, cpc = 0, ncore = 1)
  lfc = fit$summary$logFC_is_cellTRUE
  pval = fit$summary$p_is_cellTRUE
  gene = fit$summary$gene
  padj = p.adjust(pval, method = "BH")
  dplyr::tibble(gene=gene, pval=pval, padj=padj, lfc=lfc, cell_type=ct)
}


cell_types = unique(Labels$x)
ct = cell_types[1]

for (ct in cell_types) {
  Labels$is_cell = Labels$x == ct
  design_matrix = model.matrix(~is_cell, Labels)
  contrast = c(0,1)
  
  # Fit and Save
  dir.create("results/fits/", recursive = T)
  f = devil.fit(ct, cnt_mat, design_matrix, contrast, sf = "psinorm")
  saveRDS(f, paste0("results/fits/devil_", ct, ".rds"))
  f = glm.fit(ct, cnt_mat, design_matrix, contrast)
  saveRDS(f, paste0("results/fits/glmGamPoi_", ct, ".rds"))
  f = nebula.fit(ct, cnt_mat, design_matrix, contrast)
  saveRDS(f, paste0("results/fits/nebula_", ct, ".rds"))  
}