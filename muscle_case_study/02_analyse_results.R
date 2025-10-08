### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram')
sapply(pkgs, require, character.only = TRUE)


### scaDA results analysis ###

# grange_scaDA_path <- "multiomics_analysis/results/grange_annot_scADA.RDS"
# grange_scaDA <- readRDS(grange_scaDA_path)
#
# grange <- "multiomics_analysis/results/grange_annot.RDS"
# grange <- readRDS(grange)

# atac_scaDA_path <- "multiomics_analysis/results/MuscleATAC/scADA_res.RDS"
# atac_scaDA <- readRDS(atac_scaDA_path)

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

# Gene selection based on LFC & pvalue cutoff #
#quantile <- .5
lfc_cut <- 0.5
pval_cut <- .01

# atac_deg <- atac_nodup %>%
#   #dplyr::filter(FDR < 0.05, abs(log2fc) > stats::quantile(abs(atac_nodup$log2fc), quantile))
#   dplyr::filter(FDR < pval_cut, abs(log2fc) > lfc_cut)

rna_deg_devil <- rna_devil %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_devil$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)

rna_deg_glm <- rna_glm %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_glm$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)

rna_deg_nebula <- rna_nebula %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_nebula$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)


# Vienn diagram ####

x <- list(
  devil = rna_deg_devil$geneID,
  glmGamPoi = rna_deg_glm$geneID,
  nebula = rna_deg_nebula$geneID
)

venn_plot <- ggVennDiagram::ggVennDiagram(x, color = 1, lwd = 0.7) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
venn_plot

saveRDS(venn_plot, "plot/venn_plot.rds")


pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

# Gene set Enrichement analysis devil ####
path <- "results/MuscleRNA/devil_rna.RDS"
de_res <- readRDS(path)
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$name

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "All",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 500,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

devil.dot.plot <- dotplot(gseGO, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("devilS - GSEA")
devil.dot.plot
saveRDS(devil.dot.plot, "plot/devil.dot.plot.rds")
ggsave("plot/devil.dot.plot.pdf", width = 12, height = 10, units = 'in', plot = devil.dot.plot, dpi = 600)


# Gene set Enrichement analysis glmGamPoi ####
path <- "results/MuscleRNA/glmGamPoi_rna.RDS"
de_res <- readRDS(path)
de_res$adj_pval[de_res$adj_pval == 0] <- min(de_res$adj_pval[de_res$adj_pval != 0])
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$name

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "All",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 500,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

glm.dot.plot <- dotplot(gseGO, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("glmGamPoi - GSEA")
glm.dot.plot
saveRDS(glm.dot.plot, "plot/glm.dot.plot.rds")
ggsave("plot/glmGamPoi.dot.plot.pdf", width = 12, height = 10, units = 'in', plot = glm.dot.plot, dpi = 600)

# Gene set Enrichement analysis nebula ####
path <- "results/MuscleRNA/nebula_rna.RDS"
de_res <- readRDS(path)
de_res$adj_pval[de_res$adj_pval == 0] <- min(de_res$adj_pval[de_res$adj_pval != 0])
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$name

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "All",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 500,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

nebula.dot.plot <- dotplot(gseGO, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("NEBULA - GSEA")
nebula.dot.plot
saveRDS(nebula.dot.plot, "plot/nebula.dot.plot.rds")
ggsave("plot/nebula.dot.plot.pdf", width = 12, height = 10, units = 'in', plot = nebula.dot.plot, dpi = 600)

