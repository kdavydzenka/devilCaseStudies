### Results downstream analysis ###

rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr')
sapply(pkgs, require, character.only = TRUE)

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
#rna_deg_devil <- rna_deg_devil %>% dplyr::filter(adj_pval > 4.467475e-90 ,)
#rna_deg_nebula <- rna_deg_nebula %>% dplyr::filter(adj_pval > 4.787713e-22 ,)

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
  #ggplot2::theme(legend.position = 'right',
                 #legend.text = element_text(size = 16, color = "black"),
                 #legend.title = element_text(size = 18, color = "black"),  
                 #strip.text = element_text(size = 20, face = "plain", color = "black"),
                 #axis.text = element_text(size = 18, color = "black"),
                 #axis.title = element_text(size = 20, color = "black"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p_volcanos

ggsave("plot/volcanos_v2.png", dpi = 400, width = 16.0, height = 5.0, plot = p_volcanos)
saveRDS(p_volcanos, "plot/volcanos.rds")


### Gene set Enrichement analysis ###

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", 
          "clusterProfiler", "enrichplot", "ggpubr", "ReactomePA")
sapply(pkgs, require, character.only = TRUE)

enrichmentGO <- function(rna_deg_data) {
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID
  
  set.seed(1234)
  
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
    eps = 1e-10,
    by ="fgsea",
    nPerm = 10000,
    seed = 1234
  )
  return(gseGO)
  #return(gseGO@result %>% as.data.frame())
}

gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)

gseGO_list <- list(
  devil = gseGO_devil,
  glm = gseGO_glm,
  nebula = gseGO_nebula
)

gseGO_simplified <- lapply(gseGO_list, function(x) clusterProfiler::simplify(
  x, cutoff = 0.6, 
  by = "p.adjust", 
  select_fun = min, 
  measure = "Wang", 
  semData = NULL
))

gseGO_devil_s  <- gseGO_simplified$devil
gseGO_glm_s    <- gseGO_simplified$glm
gseGO_nebula_s <- gseGO_simplified$nebula


# ReactomePA

enrichmentReactomePA <- function(rna_deg_data) {
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID
  
  entrez_ids <- mapIds(org.Hs.eg.db, keys = names(genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- as.data.frame(entrez_ids)
  names(genes) <- entrez_ids$entrez_ids
  
  gseReactome <- ReactomePA::gsePathway(genes, 
                                        pvalueCutoff = 0.2,
                                        pAdjustMethod = "BH", 
                                        verbose = FALSE)
  
  return(gseReactome@result %>% as.data.frame())
}

gseRe_devil <- enrichmentReactomePA(rna_deg_devil)
gseRe_glm <- enrichmentReactomePA(rna_deg_glm)
gseRe_nebula <- enrichmentReactomePA(rna_deg_nebula)


# Plot Venn over results
ggVennDiagram::ggVennDiagram(
  list(
    glmGamPoi = gseGO_glm@result$Description,
    devil = gseGO_devil@result$Description,
    NEBULA = gseGO_nebula@result$Description
  ), color = 1, lwd = 0.8) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

gseGO_devil_s = clusterProfiler::simplify(gseGO_devil)
gseGO_glm_s = clusterProfiler::simplify(gseGO_glm)
gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula)

ggVennDiagram::ggVennDiagram(
  list(
    glmGamPoi = gseGO_glm_s@result$Description,
    devil = gseGO_devil_s@result$Description,
    NEBULA = gseGO_nebula_s@result$Description
  ), color = 1, lwd = 0.8) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

gseGO_devil_s@result$Description %>% unlist()
gseGO_glm_s@result$Description %>% unlist()
gseGO_nebula_s@result$Description %>% unlist()

nrow(gseGO_devil_s@result) / nrow(gseGO_devil@result)
nrow(gseGO_glm_s@result) / nrow(gseGO_glm@result)
nrow(gseGO_nebula_s@result) / nrow(gseGO_nebula@result)

require(rrvgo)

simplify_w_rrvgo = function(go_analysis, cut=.2) {
  simMatrix <- calculateSimMatrix(go_analysis$ID,
                                  orgdb=org.Hs.eg.db,
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=cut,
                                  orgdb=org.Hs.eg.db)  
  list(
    simMatrix=simMatrix,
    reducedTerms=reducedTerms
  )
}

rrvgo_devil = simplify_w_rrvgo(gseGO_devil)
rrvgo_glm = simplify_w_rrvgo(gseGO_glm)
rrvgo_nebula = simplify_w_rrvgo(gseGO_nebula)

heatmapPlot(rrvgo_devil$simMatrix, rrvgo_devil$reducedTerms, annotateParent = TRUE, annotationLabel = "parentTerm")
heatmapPlot(rrvgo_glm$simMatrix, rrvgo_glm$reducedTerms, annotateParent = TRUE, annotationLabel = "parentTerm")

scatterPlot(rrvgo_devil$simMatrix, rrvgo_devil$reducedTerms)
scatterPlot(rrvgo_devil$simMatrix, rrvgo_devil$reducedTerms)

unique(rrvgo_devil$reducedTerms$parentTerm)
unique(rrvgo_glm$reducedTerms$parentTerm)
rrvgo_glm$reducedTerms$cluster %>% unique()
# Enrichment results processing #

# Remove redundant terms 
source("utils.R")

filtered_devil <- remove_redundant_terms(gseGO_devil,
                                         enrichment_col = "enrichmentScore",
                                         core_col = "core_enrichment",
                                         desc_col = "Description",
                                         threshold = 0.5)

redundant_devil <- as.data.frame(filtered_devil[["redundant_data"]])
nonRed_devil <- as.data.frame(filtered_devil[["non_redundant_data"]])
nrow(redundant_devil)
nrow(nonRed_devil)

filtered_glm <- remove_redundant_terms(gseGO_glm,
                                       enrichment_col = "enrichmentScore",
                                       core_col = "core_enrichment",
                                       desc_col = "Description",
                                       threshold = 0.5)

redundant_glm <- as.data.frame(filtered_glm[["redundant_data"]])
nonRed_glm <- as.data.frame(filtered_glm[["non_redundant_data"]])
nrow(redundant_glm)
nrow(nonRed_glm)

filtered_nebula <- remove_redundant_terms(gseGO_nebula,
                                          enrichment_col = "enrichmentScore",
                                          core_col = "core_enrichment",
                                          desc_col = "Description",
                                          threshold = 0.5)

redundant_nebula <- as.data.frame(filtered_nebula[["redundant_data"]])
nonRed_nebula <- as.data.frame(filtered_nebula[["non_redundant_data"]])
nrow(redundant_nebula)
nrow(nonRed_nebula)

saveRDS(filtered_devil, file = "results/gsea_GO/gseGO_devil_list.RDS")
saveRDS(filtered_glm, file = "results/gsea_GO/gseGO_glmGamPoi_list.RDS")
saveRDS(filtered_nebula, file = "results/gsea_GO/gseGO_nebula_list.RDS")


# Enrichment results | Pathways selection #

filtered_devil <- readRDS("results/gsea_GO/gseGO_devil_list.RDS")
filtered_glm <- readRDS("results/gsea_GO/gseGO_glmGamPoi_list.RDS")
filtered_nebula <- readRDS("results/gsea_GO/gseGO_nebula_list.RDS")

filtered_devil_nonRed <- filtered_devil[["non_redundant_data"]]
filtered_glm_nonRed <- filtered_glm[["non_redundant_data"]]
filtered_nebula_nonRed <- filtered_nebula[["non_redundant_data"]]

filtered_devil_nonRed$DE_type <- ifelse(filtered_devil_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_glm_nonRed$DE_type <- ifelse(filtered_glm_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_nebula_nonRed$DE_type <- ifelse(filtered_nebula_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")

filtered_devil_nonRed <- filtered_devil_nonRed %>%
  dplyr::mutate(method = "devil")

filtered_glm_nonRed <- filtered_glm_nonRed %>%
  dplyr::mutate(method = "glmGamPoi")

filtered_nebula_nonRed <- filtered_nebula_nonRed %>%
  dplyr::mutate(method = "nebula")

combined_data <- bind_rows(filtered_devil_nonRed, 
                           filtered_glm_nonRed, 
                           filtered_nebula_nonRed)

down_terms <- combined_data %>%
  filter(DE_type == "Down-regulated")

up_terms <- combined_data %>%
  filter(DE_type == "Up-regulated") %>%
  group_by(method) %>%
  slice_max(order_by = abs(enrichmentScore), n = 10) %>% 
  ungroup()

filtered_terms <- bind_rows(down_terms, up_terms)

# Enrichment plot

plot_GO = filtered_terms %>%
  dplyr::mutate(
    Description = factor(
      Description, 
      levels = filtered_terms %>%
        arrange(factor(method, levels = c("devil", "glmGamPoi", "nebula")), enrichmentScore) %>%
        pull(Description) %>%
        unique()
    ),
    DE_type = factor(DE_type, levels = c("Up-regulated", "Down-regulated")) # Ensure DE_type is ordered
  ) %>%
  ggplot(aes(x = method, y = Description, size = setSize, color = p.adjust)) +
  geom_point() +
  #facet_grid(~DE_type, space = "free", scales = "free") +
  facet_wrap(~factor(DE_type, levels = c("Up-regulated", "Down-regulated")), scales = "free_y", ncol = 2, shrink = F) +
  scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p.adjust)") +
  theme_bw() +
  labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Count")+
  theme(
    strip.text = element_text(face = "plain", size = 16), 
    axis.text.x = element_text(size = 14, color = "black"), 
    axis.text.y = element_text(size = 14, color = "black"), 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12) 
  )
plot_GO

ggsave("plot/enrichment_dotplot.png", dpi = 400, width = 20.0, height = 9.0, plot = plot_GO)

saveRDS(plot_GO, "plot/enrichment_dotplot.RDS")



# Select not biologically specific pathways 

terms_notSpecific_devil <- c("embryo development")

terms_notSpecific_nebula <- c("response to tumor necrosis factor", 
                              "positive regulation of gene expression",
                              "embryonic morphogenesis", 
                              "tube development", 
                              "gene expression",
                              "regulation of biological process")


terms_notSpecific_glm <- c("response to gamma radiation",
                           "artery morphogenesis",
                           "positive regulation of tumor necrosis factor superfamily cytokine production",
                           "regulation of viral process",
                           "cognition",
                           "inner ear development",
                           "cell fate commitment",
                           "response to steroid hormone",
                           "viral process",
                           "response to tumor necrosis factor",
                           "reproductive structure development",
                           "kidney development",
                           "positive regulation of gene expression",
                           "negative regulation of gene expression",
                           "embryonic morphogenesis",
                           "nervous system development",
                           "central nervous system development",
                           "cardiac muscle hypertrophy")

not_biol_specific_terms <- list(terms_notSpecific_devil, terms_notSpecific_glm, terms_notSpecific_nebula)
names(not_biol_specific_terms) <- c("devil", "glmGamPoi", "nebula")

saveRDS(not_biol_specific_terms, file = "results/gsea_GO/notBiolSpecific_terms.RDS")


# Check glmGamPoi private genes enrichment #

glm_private_genes <- setdiff(rna_deg_glm$geneID, union(rna_deg_devil$geneID, rna_deg_nebula$geneID))

cat("Number of private genes in GLM:", length(glm_private_genes), "\n")

glm_private_data <- rna_deg_glm %>%
  dplyr::filter(geneID %in% glm_private_genes)

glm_private_data$adj_pval[glm_private_data$adj_pval == 0] <- min(glm_private_data$adj_pval[glm_private_data$adj_pval != 0])

gseGO_glm_private <- enrichmentGO(glm_private_data)

saveRDS(gseGO_glm_private, file = "results/gsea_GO/gseGO_glm_private.RDS")


### Summary of results across methods ###

gse_devil <- readRDS("results/gsea_GO/gseGO_devil_list.RDS")
gse_glm <- readRDS("results/gsea_GO/gseGO_glmGamPoi_list.RDS")
gse_nebula <- readRDS("results/gsea_GO/gseGO_nebula_list.RDS")

notSpecific_terms <- readRDS("results/gsea_GO/notBiolSpecific_terms.RDS")

calculate_counts <- function(redundant_data, non_redundant_data, non_specific_terms) {
  redundant_terms_count <- length(redundant_data)
  non_redundant_terms_count <- length(non_redundant_data)
  non_specific_terms_count <- length(non_specific_terms)
  
  return(c(redundant_terms_count, non_redundant_terms_count, non_specific_terms_count))
}

devil_counts <- calculate_counts(
  gse_devil[["non_redundant_data"]][["Description"]],
  gse_devil[["redundant_data"]][["Description"]],
  notSpecific_terms[["devil"]]
)

glm_counts <- calculate_counts(
  gse_glm[["non_redundant_data"]][["Description"]],
  gse_glm[["redundant_data"]][["Description"]],
  notSpecific_terms[["glmGamPoi"]]
)

nebula_counts <- calculate_counts(
  gse_nebula[["non_redundant_data"]][["Description"]],
  gse_nebula[["redundant_data"]][["Description"]],
  notSpecific_terms[["nebula"]]
)

summary_table <- data.frame(
  method = c("devil", "glmGamPoi", "nebula"),
  non_redundant_terms = c(devil_counts[1], glm_counts[1], nebula_counts[1]),
  redundant_terms = c(devil_counts[2], glm_counts[2], nebula_counts[2]),
  non_specific_terms = c(devil_counts[3], glm_counts[3], nebula_counts[3])
)


summary_table_long <- summary_table %>%
  pivot_longer(cols = -method,  
               names_to = "term_type",  
               values_to = "count")  


barplot <- ggplot(summary_table_long, aes(x = method, y = count, fill = term_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), # Position labels at the center of each section
            color = "white", # Label color
            size = 4) + # Font size
  scale_fill_manual(values = c("redundant_terms" = "#AD002AB2", 
                               "non_redundant_terms" = "#00468BB2", 
                               "non_specific_terms" = "#66CC99"),
                    labels = c("redundant_terms" = "Redundant terms", 
                               "non_redundant_terms" = "Non-Redundant terms", 
                               "non_specific_terms" = "Non-Specific terms")) +
  labs(title = "",
       x = "Method",
       y = "Count",
       fill = "Term type") +
  theme_classic()

barplot

ggsave("plot/barplot_enrichment_summary.png", dpi = 400, width = 6.0, height = 5.0, plot = barplot)
