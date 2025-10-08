### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
#set.seed(1234)
source("utils_analysis.R")

#setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")

s_cutoff = 0.75

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


# TEST 1 ####
enrichmentGO <- function(rna_deg_data) {
  rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] = min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0])
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * rna_deg_data$lfc
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID

  gseGO <- clusterProfiler::gseGO(
    genes,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    maxGSSize = 350,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 0,
    nPermSimple = 10000
    #by ="fgsea",
    #nPerm = 10000,
    #seed = 123
  )
  return(gseGO)
  #return(gseGO@result %>% as.data.frame())
}

gseGO_devil <- enrichmentGO(rna_devil)
gseGO_glm <- enrichmentGO(rna_glm)
gseGO_nebula <- enrichmentGO(rna_nebula)

# clusterProfiler::simplify analysis ####
df_simp = get_simplified_GOterms()
auc_scores = lapply(unique(df_simp$model), function(m) {
  d = df_simp %>% dplyr::filter(model == m)
  dplyr::tibble(model = m, AUC = auc(d$c, d$f))
}) %>% do.call("bind_rows", .)

simp_plot <- df_simp %>%
  dplyr::select(model, n_simplified, f, c) %>%
  tidyr::pivot_longer(c(f, n_simplified)) %>%
  dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  labs(y = "Value", x="Clustering cutoff", col="")
simp_plot


gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot

require(patchwork)
p = simp_plot | GO_plot
ggsave("plot/test_1.png", plot = p, width = 15,height = 10, dpi = 400, units = "in")

# TEST 2 ####
gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)

# clusterProfiler::simplify analysis ####
df_simp = get_simplified_GOterms()
auc_scores = lapply(unique(df_simp$model), function(m) {
  d = df_simp %>% dplyr::filter(model == m)
  dplyr::tibble(model = m, AUC = auc(d$c, d$f))
}) %>% do.call("bind_rows", .)

simp_plot <- df_simp %>%
  dplyr::select(model, n_simplified, f, c) %>%
  tidyr::pivot_longer(c(f, n_simplified)) %>%
  dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  labs(y = "Value", x="Clustering cutoff", col="")
simp_plot


gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot

require(patchwork)
p = simp_plot | GO_plot
ggsave("plot/test_2.png", plot = p, width = 15,height = 10, dpi = 400, units = "in")


# TEST 3 ####
enrichmentGO <- function(rna_deg_data) {
  rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] = min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0])
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * rna_deg_data$lfc
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID

  gseGO <- clusterProfiler::gseGO(
    genes,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    maxGSSize = 350,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 0,
    nPermSimple = 100000
    #by ="fgsea",
    #nPerm = 10000,
    #seed = 123
  )
  return(gseGO)
  #return(gseGO@result %>% as.data.frame())
}

gseGO_devil <- enrichmentGO(rna_devil)
gseGO_glm <- enrichmentGO(rna_glm)
gseGO_nebula <- enrichmentGO(rna_nebula)

# clusterProfiler::simplify analysis ####
df_simp = get_simplified_GOterms()
auc_scores = lapply(unique(df_simp$model), function(m) {
  d = df_simp %>% dplyr::filter(model == m)
  dplyr::tibble(model = m, AUC = auc(d$c, d$f))
}) %>% do.call("bind_rows", .)

simp_plot <- df_simp %>%
  dplyr::select(model, n_simplified, f, c) %>%
  tidyr::pivot_longer(c(f, n_simplified)) %>%
  dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  labs(y = "Value", x="Clustering cutoff", col="")
simp_plot


gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot

require(patchwork)
p = simp_plot | GO_plot
ggsave("plot/test_3.png", plot = p, width = 15,height = 10, dpi = 400, units = "in")

# TEST 4 ####
gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)

# clusterProfiler::simplify analysis ####
df_simp = get_simplified_GOterms()
auc_scores = lapply(unique(df_simp$model), function(m) {
  d = df_simp %>% dplyr::filter(model == m)
  dplyr::tibble(model = m, AUC = auc(d$c, d$f))
}) %>% do.call("bind_rows", .)

simp_plot <- df_simp %>%
  dplyr::select(model, n_simplified, f, c) %>%
  tidyr::pivot_longer(c(f, n_simplified)) %>%
  dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  labs(y = "Value", x="Clustering cutoff", col="")
simp_plot


gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot

require(patchwork)
p = simp_plot | GO_plot
ggsave("plot/test_4.png", plot = p, width = 15,height = 10, dpi = 400, units = "in")
