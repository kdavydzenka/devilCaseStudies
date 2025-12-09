
rm(list = ls())
library(dplyr)
library(readr)
library(ggpubr)
library(stringr)
library(ggplot2)
library(tidyr)
library(pROC)

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "nebula" =  "steelblue",
  "NEBULA" =  "steelblue",
  "devil" = "#099668"
)

syn_map <- tribble(
  ~panglao,                      ~target,
  
  # --- B cells ---
  "B cells",                     "Naive_B",
  "B cells naive",               "Naive_B",
  "Naive B cells",               "Naive_B",
  "B cells memory",              "Memory_B",
  "Memory B cells",              "Memory_B",
  "Plasma cells",                "Memory_B",    
  
  # --- Monocytes / myeloid ---
  "Monocytes",                   "Monocytes",
  "Macrophages",                 "Monocytes",
  "Kupffer cells",               "Monocytes",
  "Alveolar macrophages",        "Monocytes",
  "Dendritic cells",             "Monocytes",
  "Plasmacytoid dendritic cells","Monocytes",
  
  # --- Granulocytes ---
  "Neutrophils",                 "Granulocyte",
  "Eosinophils",                 "Granulocyte",
  "Basophils",                   "Granulocyte",
  "Mast cells",                  "Granulocyte",
  
  # --- NK / innate lymphoid ---
  "NK cells",                    "NK",
  "Natural killer T cells",      "NK",
  "Nuocytes",                    "NK",
  
  # --- T cells: general / helper (CD4 lineage) ---
  "T cells",                     "Naive_CD4_T",
  # "T cells naive",               "Naive_CD4_T",
  "T helper cells",              "Memory_CD4_T",
  "T memory cells",              "Memory_CD4_T",
  "T follicular helper cells",   "Memory_CD4_T",
  "T regulatory cells",          "Memory_CD4_T",
  
  # --- T cells: cytotoxic (CD8 lineage) ---
  "T cytotoxic cells",           "Memory_CD8_T",
  "Cytotoxic T cells",           "Memory_CD8_T",
  "T cells naive",               "Naive_CD8_T",
  "CD8 T cells",                 "Naive_CD8_T",
  "Naive CD8 T cells",           "Naive_CD8_T",
  "Effector CD8 T cells",        "Memory_CD8_T",
  "Memory CD8 T cells",          "Memory_CD8_T",
  
  # --- Gamma-delta T ---
  "Gamma delta T cells",         "Gamma_delta_T"
)

evaluate_markers <- function(df, ref, top_n = 50) {
  ct = unique(df$cell_type)[1]
  results <- lapply(unique(df$cell_type), function(ct) {
    print(ct)
    de_ct <- df %>% filter(cell_type == ct)
    ref_ct <- ref %>% filter(cell_type == ct)
    
    if (nrow(ref_ct) == 0) return(NULL)
    
    # 1. Rank DE genes by adjusted p-value (or abs(lfc))
    de_ct <- de_ct %>% 
      dplyr::filter(lfc > 0) %>% 
      dplyr::arrange(padj, -lfc) %>%
      dplyr::mutate(rank = row_number(),
             is_marker = gene %in% ref_ct$gene)
    
    # 2. Overlap metrics
    top_genes <- de_ct %>% slice_head(n = top_n)
    n_ref <- length(unique(ref_ct$gene))
    n_top <- nrow(top_genes)
    overlap <- sum(top_genes$is_marker)
    precision <- overlap / n_top
    recall <- overlap / n_ref
    
    # 3. Enrichment Fisher test (top N)
    tbl <- matrix(c(overlap,
                    n_top - overlap,
                    n_ref - overlap,
                    nrow(de_ct) - n_ref - n_top + overlap),
                  nrow = 2)
    fisher_p <- fisher.test(tbl, alternative="greater")$p.value
    
    # 4. AUC: do DE rank vs marker binary label
    # de_ct$score <- -log10(de_ct$padj + 1e-10)
    de_ct$score <- -log10(de_ct$padj + 1e-10) * de_ct$lfc
    auc <- tryCatch({
      roc(response = de_ct$is_marker, predictor = de_ct$score, quiet = TRUE)$auc
    }, error = function(e) NA_real_)
    auc = as.numeric(auc)
    
    tibble(cell_type = ct,
           top_n = top_n,
           overlap = overlap,
           precision = precision,
           recall = recall,
           fisher_p = fisher_p,
           auc = auc)
  })
  bind_rows(results)
}

# Example: PanglaoDB from decoupler or local file
pang <- decoupleR::get_resource(name = "PanglaoDB", organism = "human", license = "academic")
pang <- pang %>% transmute(source="PanglaoDB", gene=genesymbol, cell_type=cell_type)

pang_mapped <- pang %>%
  left_join(syn_map, by = c("cell_type" = "panglao")) %>%
  mutate(target = coalesce(target, cell_type)) %>%  # keep unmapped original names
  filter(target %in% c(
    "Gamma_delta_T", "Granulocyte", "Memory_B",
    "Memory_CD4_T", "Memory_CD8_T", "Monocytes",
    "Naive_B", "Naive_CD4_T", "Naive_CD8_T", "NK"
  )) %>%
  select(gene, cell_type = target) %>%
  distinct()
rm(pang)
# Read df

models = c("devil", "glmGamPoi", "nebula")
m = "nebula"
df_res = lapply(models, function(m) {
  file.paths = list.files("results/fits/", full.names = T)[grepl(m, list.files("results/fits/"))]
  df = lapply(file.paths, function(x) readRDS(x)) %>% do.call("bind_rows", .)
  
  # Keep only the cell types relevant to your dataset
  target_labels <- unique(df$cell_type)
  pang_mapped <- pang_mapped %>% filter(cell_type %in% target_labels)
  evaluate_markers(df, ref = pang_mapped, top_n = 50) %>% dplyr::mutate(model = m)
}) %>% do.call("bind_rows", .)
df_res = df_res %>% dplyr::select(cell_type, auc, model)

auc_bar_plot = df_res %>% 
  ggplot2::ggplot(mapping = aes(x = cell_type, y = auc, fill = model)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = method_colors) +
  labs(y = "AUC", x = "Cell type", fill = "Method") +
  coord_flip()
auc_bar_plot
saveRDS(auc_bar_plot, "img/auc_bar_plot.RDS")

methods <- unique(df_res$model)
comparisons <- combn(methods, 2, simplify = FALSE)
auc_boxplot = df_res %>% 
  ggplot2::ggplot(mapping = aes(x = model, y = auc, col = model)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Method", y = "AUC") +
  scale_color_manual(values = method_colors) +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  )
auc_boxplot
saveRDS(auc_boxplot, "img/auc_boxplot.RDS")

# Save summary
df_summary_mean = df_res %>% 
  dplyr::group_by(model) %>% 
  dplyr::summarise(mean_auc = mean(auc))
write_csv(df_summary_mean, "results/summary_mean_auc.csv")

df_summary_auc = df_res %>% 
  tidyr::pivot_wider(names_from = model, values_from = auc)
write_csv(df_summary_mean, "results/summary_auc.csv")


# Prepare PDF
umap = readRDS("img/umap_plot.RDS")
auc_box = readRDS("img/auc_boxplot.RDS")
auc_barplot = readRDS("img/auc_bar_plot.RDS")

des = "
AABBC
"

library(patchwork)
p = free(umap) + free(auc_barplot) + free(auc_box) + 
  patchwork::plot_layout(design = des) + 
  patchwork::plot_annotation(tag_levels = c("A"))

saveRDS(p, "img/10_pbmc_plot.RDS")
ggsave("img/10_pbmc_plot.pdf", p, width = 12, height = 5)  
