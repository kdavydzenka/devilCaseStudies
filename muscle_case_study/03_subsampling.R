rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("tidyverse")
sapply(pkgs, require, character.only = TRUE)
#source("utils/utils_analysis.R")

###------------Load results---------------###

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

read_rna <- function(method, condition, set = c("subsampled", "full")) {
  set <- match.arg(set)
  path <- glue::glue("results/MuscleRNA/{set}/{method}_{condition}.RDS")
  readRDS(path)
}

fix_gene_column <- function(df) {
  if ("name" %in% names(df)) {
    df <- df %>% rename(geneID = name)
  }
  df
}

all_results <- cross_df(list(
  method = methods,
  condition = conditions,
  set = c("full", "subsampled")
)) %>%
  mutate(
    data = pmap(list(method, condition, set), read_rna),
    data = map(data, fix_gene_column)
  )

results_long <- all_results %>%
  unnest(data)

###----------Compare lfc & adj_pvalue----------###

comparison <- results_long %>%
  dplyr::select(geneID, method, condition, set, lfc, adj_pval) %>%
  pivot_wider(
    names_from = set,
    values_from = c(lfc, adj_pval),
    names_sep = "_"
  )

comparison <- comparison %>%
  dplyr::mutate(
    delta_lfc = lfc_subsampled - lfc_full,
    delta_adj_pval = adj_pval_subsampled - adj_pval_full
  )

summary_stats <- comparison %>%
  group_by(method, condition) %>%
  summarise(
    cor_lfc = cor(lfc_full, lfc_subsampled, use = "pairwise"),
    cor_adj_pval = cor(adj_pval_full, adj_pval_subsampled, use = "pairwise"),
    mean_abs_lfc_diff = mean(abs(lfc_full - lfc_subsampled), na.rm = TRUE),
    mean_abs_padj_diff = mean(abs(adj_pval_full - adj_pval_subsampled), na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  )

saveRDS(summary_stats, file = "results/subsampled_test/summary_stats.RDS")
saveRDS(comparison, file = "results/subsampled_test/joint_results.RDS")


### Scatter plot per method × condition ###

plot_lfc <- comparison %>%
  ggplot(aes(lfc_full, lfc_subsampled)) +
  geom_point(alpha = 0.3) +
  facet_grid(method ~ condition) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_bw() +
  labs(
    title = "Full vs Subsampled lfc comparison",
    x = "log2FC (full)",
    y = "log2FC (subsampled)"
  )
plot_lfc

ggsave("plot/revision/scatter_comparison_lfc.png", dpi = 400, width = 12.0, height = 15.0, plot = plot_lfc)

plot_padj <- comparison %>%
  ggplot(aes(-log10(adj_pval_full), -log10(adj_pval_subsampled))) +
  geom_point(alpha = 0.3) +
  #facet_grid(method ~ condition, scales = "free") +
  ggh4x::facet_nested(factor(method, levels = c("devil", "glmGamPoi", "nebula"))
                      ~factor(condition, levels = c("age_only", "age_type1", "age_type2", "interaction")), 
                      scales ="free", independent = "x")+
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_bw()+
  labs(
    title = "Full vs Subsampled adj_pval comparison",
    x = "-log10(adj pvalue full)",
    y = "-log10(adj pvalue subsampled)"
  )
plot_padj

ggsave("plot/revision/scatter_comparison_padj.png", dpi = 400, width = 12.0, height = 15.0, plot = plot_padj)
