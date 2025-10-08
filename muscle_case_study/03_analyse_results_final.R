### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db", "GOSemSim")
sapply(pkgs, require, character.only = TRUE)

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue',
  "devil" = "#099668"
)

source("utils_analysis.R")

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")

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


# Vienn diagram ####

x <- list(
  devil = rna_deg_devil$geneID,
  glmGamPoi = rna_deg_glm$geneID,
  nebula = rna_deg_nebula$geneID
)

venn_plot <- ggVennDiagram::ggVennDiagram(x, color = 1, lwd = 0.8) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
venn_plot

saveRDS(venn_plot, "plot/venn_plot_v2.rds")


# Volcano plot
rna_deg_devil$adj_pval[rna_deg_devil$adj_pval == 0] <- min(rna_deg_devil$adj_pval[rna_deg_devil$adj_pval != 0])
rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval == 0] <- min(rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval != 0])
rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

rna_join <- rbind(rna_deg_devil, rna_deg_glm, rna_deg_nebula)
rna_join <- rna_join %>%
  dplyr::mutate(
    isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
    DEtype = if_else(!isDE, "n.s.", if_else(lfc > 0, "Up-reg", "Down-reg")))

de_colors <- c("Down-reg" = "steelblue", "Up-reg" = "indianred", "n.s." = "grey")

gene_markers <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                  "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                  "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

saveRDS(gene_markers, file = "results/gsea_GO/gene_markers_BRelevant.RDS")

p_volcanos <- rna_join %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(col = DEtype), size = 2.0, alpha = 0.2) +
  scale_color_manual(values = de_colors) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  geom_label_repel(
    data = rna_join %>% filter(geneID %in% gene_markers),
    aes(label = geneID),
    size = 4.0,
    #fontface = "bold",
    color = "black",
    fill = "white",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.5,
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.4, "lines"),
    min.segment.length = 1
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)),
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
  labs(x = expression(Log[2] ~ FC),
       y = expression(-log[10] ~ Pvalue),
       col = "DE type") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p_volcanos

ggsave("plot/volcanos_v2.png", dpi = 400, width = 16.0, height = 5.0, plot = p_volcanos)
saveRDS(p_volcanos, "plot/volcanos.rds")


### Gene set Enrichement analysis ###

gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)


## clusterProfiler::simplify analysis ##
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
  dplyr::mutate(name = ifelse(name=="f", "Fraction of Terms Retained", "Reduced Terms Count")) %>%
  ggplot(mapping = aes(x=c, y=value, col=model)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~name, scales = "free", ncol = 1,strip.position = "left") +
  labs(y = "", x="Clustering cutoff", col="") +
  scale_color_manual(values = method_colors)
simp_plot
saveRDS(simp_plot, file = "plot/simp_plot.rds")

s_cutoff = 0.6
gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=s_cutoff)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=s_cutoff)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=s_cutoff)

devil_s <- gseGO_devil_s@result
glm_s <- gseGO_glm_s@result
nebula_s <- gseGO_nebula_s@result

saveRDS(gseGO_devil_s, "results/gsea_GO/gseGO_devil_s.RDS")
saveRDS(gseGO_devil, "results/gsea_GO/gseGO_devil.RDS")

saveRDS(gseGO_glm_s, "results/gsea_GO/gseGO_glm_s.RDS")
saveRDS(gseGO_glm, "results/gsea_GO/gseGO_glm.RDS")

saveRDS(gseGO_nebula_s, "results/gsea_GO/gseGO_nebula_s.RDS")
saveRDS(gseGO_nebula, "results/gsea_GO/gseGO_nebula.RDS")

gseGO_devil_s <- readRDS("results/gsea_GO/gseGO_devil_s.RDS")
gseGO_glm_s <- readRDS("results/gsea_GO/gseGO_glm_s.RDS")
gseGO_nebula_s <- readRDS("results/gsea_GO/gseGO_nebula_s.RDS")

GO_plot = plot_dotplot_GO(gseGO_devil_s@result, gseGO_glm_s@result, gseGO_nebula_s@result)
GO_plot

saveRDS(GO_plot, "plot/enrichment_dotplot.RDS")
ggsave("plot/enrichment_dotplot.png", dpi = 400, width = 10.0, height = 9.0, plot = GO_plot)



### glmGamPoi private genes/ shared genes enrichment - to test for potential overcalling / permissiveness ###

run_enrichment <- function(df, gene_ids) {
  filtered <- df %>% dplyr::filter(geneID %in% gene_ids)
  enrichmentGO(filtered)
}

rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

glm_private_g <- setdiff(rna_deg_glm$geneID, union(rna_deg_devil$geneID, rna_deg_nebula$geneID))
shared_glm_devil <- intersect(rna_deg_glm$geneID, rna_deg_devil$geneID)

gseGO_glm_private <- run_enrichment(rna_deg_glm, glm_private_g)
gseGO_shared_glm  <- run_enrichment(rna_deg_glm, shared_glm_devil)
gseGO_shared_devil <- run_enrichment(rna_deg_devil, shared_glm_devil)

#gseGO_glm_private_s = clusterProfiler::simplify(gseGO_glm_private, cutoff=0.6)

saveRDS(gseGO_glm_private, "results/gsea_GO/gseGO_glm_private.RDS")


# Test only shared pathways
go_glm <- gseGO_shared_glm@result$ID
go_devil <- gseGO_shared_devil@result$ID

shared_terms <- intersect(go_glm, go_devil)
private_glm  <- setdiff(go_glm, go_devil)
private_devil <- setdiff(go_devil, go_glm)

shared_df <- gseGO_shared_glm@result[gseGO_shared_glm@result$ID %in% shared_terms, ]
shared_df$Group <- "Shared"

private_glm_df <- gseGO_shared_glm@result[gseGO_shared_glm@result$ID %in% private_glm, ]
private_glm_df$Group <- "Private glmGamPoi"

private_devil_df <- gseGO_shared_devil@result[gseGO_shared_devil@result$ID %in% private_devil, ]
private_devil_df$Group <- "Private devil"

go_terms_combined <- rbind(shared_df, private_glm_df, private_devil_df)

# Calculate semantic specificity 

hsGO_BP <- godata('org.Hs.eg.db', ont = "BP")

calculate_specificity_scores <- function(gsea_result, semData) {
  go_terms <- gsea_result@result$ID
  sim_matrix <- GOSemSim::mgoSim(go_terms, go_terms, semData = semData, measure = "Wang", combine = NULL)
  specificity_scores <- rowMeans(sim_matrix, na.rm = TRUE)
  gsea_result@result$specificity <- specificity_scores
  return(gsea_result)
}

gseGO_glm_private <- calculate_specificity_scores(gseGO_glm_private, hsGO_BP)
gseGO_shared_glm  <- calculate_specificity_scores(gseGO_shared_glm, hsGO_BP)
gseGO_shared_devil <- calculate_specificity_scores(gseGO_shared_devil, hsGO_BP)

shared_combined <- rbind(gseGO_shared_glm@result, gseGO_shared_devil@result)
shared_combined <- shared_combined[!duplicated(shared_combined$ID), ]

# Specificity comparison plot

specificity_comparison <- data.frame(
  GeneSet = rep(c("Private glmGamPoisson", "Shared glmGamPoisson/devil"),
                c(nrow(gseGO_glm_private@result), nrow(shared_combined))),
  Specificity = c(gseGO_glm_private@result$specificity, shared_combined$specificity)
)

specificity_boxplot <- ggplot(specificity_comparison, aes(x = GeneSet, y = Specificity, fill = GeneSet)) +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = list(c("Private glmGamPoisson", "Shared glmGamPoisson/devil")),
    label = "p.format"
  ) +
  scale_fill_manual(values = c("Private glmGamPoisson" = "#f99379", 
                               "Shared glmGamPoisson/devil" = "#099668")) +
  labs(title = "Comparison of Specificity Scores",
       x = "Gene set",
       y = "Specificity Score",
       fill = "Gene set group") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

specificity_boxplot

ggsave("plot/specificity_score_comparison_onlyShared.png", dpi = 400, width = 6.0, height = 4.0, plot = specificity_boxplot)
ggsave("plot/specificity_score_comparison.png", dpi = 400, width = 7.0, height = 5.0, plot = specificity_boxplot)
saveRDS(specificity_boxplot, "plot/specificity_boxplot.RDS")


# GO terms overlap #

terms_glm <- gseGO_shared_glm@result$ID
terms_devil <- gseGO_shared_devil@result$ID

venn_data <- list(glmGamPoisson = terms_glm, devil = terms_devil)

venn <- ggVennDiagram(venn_data, color = 1, lwd = 0.8,
                      set_label = list(size = 2))+
  scale_fill_gradient(low = "#F9E0B0", high = "#BF8B4B")+
  theme(legend.position = "none")
venn

saveRDS(venn, "plot/venn_GOterms_overplap.rds")
ggsave("plot/venn_GOterms_overplap.png", venn, width = 6, height = 6)

# ReactomePA enrichment
#gseRe_devil <- enrichmentReactomePA(rna_deg_devil)
#gseRe_glm <- enrichmentReactomePA(rna_deg_glm)
#gseRe_nebula <- enrichmentReactomePA(rna_deg_nebula)

#saveRDS(gseRe_devil, "results/gse_Reactome/gseRE_devil.RDS")
#saveRDS(gseRe_glm, "results/gse_Reactome/gseRE_glm.RDS")
#saveRDS(gseRe_nebula, "results/gse_Reactome/gseRE_nebula.RDS")

# Plot
#gseRe_devil <- readRDS("results/gse_Reactome/gseRE_devil.RDS")
#gseRe_glm <- readRDS("results/gse_Reactome/gseRE_glm.RDS")
#gseRe_nebula <- readRDS("results/gse_Reactome/gseRE_nebula.RDS")

#RE_plot = plot_dotplot_RE(gseRe_devil, gseRe_glm, gseRe_nebula)
#RE_plot
#saveRDS(RE_plot, "plot/enrichment_dotplot_RE.RDS")

#ggsave("plot/enrichment_dotplot_RE.png", dpi = 400, width = 10.0, height = 9.0, plot = RE_plot)

