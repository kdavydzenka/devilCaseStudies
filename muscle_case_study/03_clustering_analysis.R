### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
#set.seed(1234)
source("utils_analysis.R")

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "devil" = "#099668"
)

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

rna_deg_devil$adj_pval[rna_deg_devil$adj_pval == 0] <- min(rna_deg_devil$adj_pval[rna_deg_devil$adj_pval != 0])
rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval == 0] <- min(rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval != 0])
rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

### Gene set Enrichement analysis ###

gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)

# clusterProfiler::simplify analysis ####
df_simp = get_simplified_GOterms()
saveRDS(df_simp, "results/gsea_GO/simplified_df.RDS")
df_simp = readRDS("results/gsea_GO/simplified_df.RDS")

auc_scores = lapply(unique(df_simp$model), function(m) {
  d = df_simp %>% dplyr::filter(model == m)
  dplyr::tibble(model = m, AUC = auc(d$c, d$f))
}) %>% do.call("bind_rows", .)
saveRDS(auc_scores, "results/gsea_GO/auc_scores.RDS")
auc_scores = readRDS("results/gsea_GO/auc_scores.RDS")

simp_plot <- df_simp %>%
  dplyr::select(model, n_simplified, f, c) %>%
  tidyr::pivot_longer(c(f, n_simplified)) %>%
  dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  dplyr::filter(name == "Fraction simplified") %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  scale_color_manual(values = method_colors) +
  labs(y = "Value", x="Clustering cutoff", col="")
simp_plot
saveRDS(simp_plot, file = "plot/simp_plot.rds")

s_cutoff = 0.5
gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

saveRDS(gseGO_devil_s, "results/gsea_GO/gseGO_devil_s.RDS")
saveRDS(gseGO_devil, "results/gsea_GO/gseGO_devil.RDS")

saveRDS(gseGO_glm_s, "results/gsea_GO/gseGO_glm_s.RDS")
saveRDS(gseGO_glm, "results/gsea_GO/gseGO_glm.RDS")

saveRDS(gseGO_nebula_s, "results/gsea_GO/gseGO_nebula_s.RDS")
saveRDS(gseGO_nebula, "results/gsea_GO/gseGO_nebula.RDS")

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot
saveRDS(GO_plot, "plot/enrichment_dotplot.RDS")
