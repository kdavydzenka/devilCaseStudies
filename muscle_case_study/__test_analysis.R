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

# Loading data #
rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

rna_devil$adj_pval[rna_devil$adj_pval == 0] <- min(rna_devil$adj_pval[rna_devil$adj_pval != 0])
rna_nebula$adj_pval[rna_nebula$adj_pval == 0] <- min(rna_nebula$adj_pval[rna_nebula$adj_pval != 0])
rna_glm$adj_pval[rna_glm$adj_pval == 0] <- min(rna_glm$adj_pval[rna_glm$adj_pval != 0])

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
saveRDS(venn_plot, "plot/final_venn_plot.rds")

# P-values plot
p_values_dist_plot = dplyr::bind_rows(
  dplyr::tibble(method = "devil", pvalues = rna_devil$adj_pval),
  dplyr::tibble(method = "glmGamPoi", pvalues = rna_glm$adj_pval),
  dplyr::tibble(method = "NEBULA", pvalues = rna_nebula$adj_pval),
) %>% 
  ggplot(mapping = aes(x = method, y = pvalues, col = method)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  labs(x = "", y="Adjusted p-value", col="")
saveRDS(p_values_dist_plot, "plot/final_p_values_dist_plot.rds")

enrichmentGO = function(rna_deg_data) {
  rna_deg_data$pval[rna_deg_data$pval == 0] = min(rna_deg_data$pval[rna_deg_data$pval != 0])
  rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] = min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0])
  #rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] = min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0])
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data$RankMetric <- -log10(rna_deg_data$pval) * sign(rna_deg_data$lfc)
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * rna_deg_data$lfc
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID
  
  gseGO <- clusterProfiler::gseGO(
    genes,
    ont = "BP",  
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    #maxGSSize = 200,
    maxGSSize = 350,
    #maxGSSize = 500,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05, 
    seed = TRUE,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 0,
    nPermSimple = 50000
    #nPermSimple = 10000
    #by ="fgsea",
    #nPerm = 10000,
    #seed = 123
  )
  return(gseGO)
  #return(gseGO@result %>% as.data.frame())
}

gseGO_devil <- enrichmentGO(rna_devil %>% dplyr::filter(adj_pval <= .05))
gseGO_glm <- enrichmentGO(rna_glm %>% dplyr::filter(adj_pval <= .05))
gseGO_nebula <- enrichmentGO(rna_nebula %>% dplyr::filter(adj_pval <= .05))

saveRDS(gseGO_devil, "test_analysis_res/gseGO_devil.rds")
saveRDS(gseGO_glm, "test_analysis_res/gseGO_glm.rds")
saveRDS(gseGO_nebula, "test_analysis_res/gseGO_nebula.rds")

gseGO_devil@result %>% nrow()
gseGO_glm@result %>% nrow()
gseGO_nebula@result %>% nrow()

plot_dotplot = function(topK = 10) {
  go_terms <- c(
    "myofibril assembly",
    "striated muscle cell development",
    "actin-mediated cell contraction",
    "actin filament-based movement",
    "muscle contraction",
    "muscle system process",
    "muscle cell development",
    "cell-cell adhesion via plasma-membrane adhesion molecules",
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "acute-phase response",
    "regulation of cell division",
    "immune response-activating cell surface receptor signaling pathway",
    "defense response to other organism",
    "response to biotic stimulus",
    "defense response to symbiont",
    "response to external biotic stimulus",
    "response to other organism",
    "neurotransmitter reuptake",
    "cellular component assembly involved in morphogenesis",
    "amelogenesis",
    "actomyosin structure organization",
    "striated muscle cell differentiation",
    "regulation of extrinsic apoptotic signaling pathway via death domain receptors",
    "regulation of membrane lipid metabolic process",
    "negative regulation of protein localization to nucleus",
    "cellular response to fibroblast growth factor stimulus",
    "negative regulation of extrinsic apoptotic signaling pathway via death domain receptors",
    "response to fibroblast growth factor",
    "cardiac atrium development",
    "clathrin-dependent endocytosis",
    "regulation of extrinsic apoptotic signaling pathway",
    "negative regulation of extrinsic apoptotic signaling pathway",
    "innate immune response",
    "leukocyte migration",
    "nitrogen compound transport",
    "cytokine production",
    "regulation of cytokine production",
    "defense response",
    "tube development",
    "cellular response to cytokine stimulus",
    "positive regulation of gene expression",
    "response to cytokine"
  )
  
  clusters <- c(
    # Muscle development & contraction
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Cell Adhesion",
    "Cell Adhesion",
    "Immune and Defense Response",
    "Cell Cycle, Apoptosis, and Signaling",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Transport and Endocytosis",
    "Morphogenesis and Tissue Organization",
    "Organ Development (non-muscle)",
    "Muscle Development and Contraction",
    "Muscle Development and Contraction",
    "Cell Cycle, Apoptosis, and Signaling",
    "Other / General Regulation",
    "Other / General Regulation",
    "Response to signaling molecules",
    "Cell Cycle, Apoptosis, and Signaling",
    "Response to signaling molecules",
    "Organ Development (non-muscle)", 
    "Transport and Endocytosis",
    "Cell Cycle, Apoptosis, and Signaling",
    "Cell Cycle, Apoptosis, and Signaling",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Transport and Endocytosis",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Immune and Defense Response",
    "Organ Development (non-muscle)",
    "Response to signaling molecules",
    "Other / General Regulation",
    "Response to signaling molecules"
  )
  
  redundant_terms = list(
    "cell-cell adhesion via plasma-membrane adhesion molecules" = c("cell-cell adhesion via plasma-membrane adhesion molecules", "homophilic cell adhesion via plasma membrane adhesion molecules"),
    "regulation of extrinsic apoptotic signaling pathway" = c("regulation of extrinsic apoptotic signaling pathway", "regulation of extrinsic apoptotic signaling pathway via death domain receptors", "negative regulation of extrinsic apoptotic signaling pathway via death domain receptors", "negative regulation of extrinsic apoptotic signaling pathway"),
    "response to fibroblast growth factor" = c("regulation of extrinsic apoptotic signaling pathway", "cellular response to fibroblast growth factor stimulus"),
    "response to cytokine" = c("response to cytokine", "cellular response to cytokine stimulus"),
    "response to biotic stimulus" = c("response to external biotic stimulus","response to biotic stimulus"),
    "defense response" = c("defense response", "defense response to symbiont", "defense response to other organism")
  )
  
  go_cluster_map = c(clusters)
  names(go_cluster_map) = go_terms
  
  devil_res = gseGO_devil@result %>% 
    dplyr::group_by(ID) %>% 
    dplyr::mutate(GeneRatio = (sum(str_count(core_enrichment, "/")) + 1) / setSize) %>% 
    dplyr::mutate(sign = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>% 
    dplyr::arrange(-GeneRatio) %>% 
    dplyr::group_by(sign) %>% 
    dplyr::slice_head(n = topK)
  
  glm_res = gseGO_glm@result %>% 
    dplyr::group_by(ID) %>% 
    dplyr::mutate(GeneRatio = (sum(str_count(core_enrichment, "/")) + 1) / setSize) %>% 
    dplyr::mutate(sign = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>% 
    dplyr::arrange(-GeneRatio) %>% 
    dplyr::group_by(sign) %>% 
    dplyr::slice_head(n = topK)
  
  nebula_res = gseGO_nebula@result %>% 
    dplyr::group_by(ID) %>% 
    dplyr::mutate(GeneRatio = (sum(str_count(core_enrichment, "/")) + 1) / setSize) %>% 
    dplyr::mutate(sign = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>% 
    dplyr::arrange(-GeneRatio) %>% 
    dplyr::group_by(sign) %>% 
    dplyr::slice_head(n = topK)
  
  df = dplyr::bind_rows(
    devil_res %>% dplyr::mutate(name = "devil"),
    glm_res %>% dplyr::mutate(name = "glmGamPoi"),
    nebula_res %>% dplyr::mutate(name = "NEBULA")
  )
  
  df$Description = lapply(1:nrow(df), function(i) {
    r = df[i,]
    for (j in 1:length(redundant_terms)) {
      if (r$Description %in% redundant_terms[[j]]) {
        return(names(redundant_terms)[j])
      }
    }
    return(r$Description)  
  }) %>% unlist()
  
  df = df %>% dplyr::group_by(Description, name) %>% dplyr::filter(GeneRatio == max(GeneRatio))
  df$Biological_process = go_cluster_map[df$Description]
  
  ggplot(df, mapping = aes(x = name, y=Description, size = GeneRatio, col=pvalue )) +
    geom_point() +
    facet_grid(Biological_process~sign, space = "free", scales = "free") +
    scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p-value") +
    theme_bw() +
    labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Ratio") +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
}

GO_plot = plot_dotplot()

GO_plot
saveRDS(GO_plot, "plot/final_GO_plot.rds")

# Load required packages
library(clusterProfiler)
library(TissueEnrich)
library(org.Hs.eg.db)

get_tissue_specific_res = function(gse_result, method) {
  # Extract gene sets (significant ones)
  significant_terms <- subset(gse_result@result, qvalue < 0.05)
  
  # Get the genes for each GO term (pathway)
  # Map ENTREZ IDs to SYMBOLs if needed
  all_gene_sets <- lapply(significant_terms$ID, function(go_id) {
    genes <- DOSE::geneInCategory(gse_result)[[go_id]]
    genes
  })
  names(all_gene_sets) <- significant_terms$Description
  
  # Run TissueEnrich on each gene set
  results <- lapply(seq_along(all_gene_sets), function(i) {
    gene_set <- all_gene_sets[[i]]
    gs <- GeneSet(geneIds=gene_set,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
    if (length(gene_set) >= 5) { # Minimum size
      TissueEnrich::teEnrichment(gs)
    } else {
      NULL
    }
  })
  names(results) <- names(all_gene_sets)
  
  # Extract enrichment score (p-value or fold change) for the tissue of interest
  res_df = lapply(1:length(results), function(i) {
    path_name = names(results)[i]
    te_result = results[[i]]
    if (is.null(te_result)) {
      print(i)
      return(NULL)
    }
    seEnrichmentOutput<-te_result[[1]]
    enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
    enrichmentOutput$Tissue<-row.names(enrichmentOutput)
    enrichmentOutput %>% dplyr::mutate(path_name = path_name)
  }) %>% do.call("bind_rows", .)
  rownames(res_df) = NULL
  res_df %>% dplyr::mutate(method = method)
}

tissue = target_tissue <- "Skeletal Muscle"
tissue_gse_devil = get_tissue_specific_res(gseGO_devil, "devil")
tissue_gse_glm = get_tissue_specific_res(gseGO_glm, "glmGamPoi")
tissue_gse_nebula = get_tissue_specific_res(gseGO_nebula, "NEBULA")

saveRDS(tissue_gse_devil, "test_analysis_res/tissue_gse_devil.rds")
saveRDS(tissue_gse_glm, "test_analysis_res/tissue_gse_glm")
saveRDS(tissue_gse_nebula, "test_analysis_res/tissue_gse_nebula")

pval_cut = .05
best_df = lapply(list(tissue_gse_devil, tissue_gse_glm, tissue_gse_nebula), function(tissue_gse) {
  print(unique(tissue_gse$method))
  path = "response to other organism"
  lapply(unique(tissue_gse$path_name), function(path) {
    r = tissue_gse %>% 
      dplyr::filter(path_name == path) %>% 
      dplyr::filter(Log10PValue >= -log10(pval_cut))
    
    if (nrow(r) > 0) {
      r %>% 
        na.omit() %>% 
        dplyr::filter(fold.change == max(fold.change))
      if (nrow(r) > 1) {
        r = r %>% dplyr::filter(Tissue.Specific.Genes == max(Tissue.Specific.Genes))
        r = r[1,]
      }
    } else {
      r = tissue_gse %>% 
        dplyr::filter(path_name == path) %>% 
        dplyr::sample_n(1) %>% 
        dplyr::mutate(Tissue = "Generic")
    }
    r
  }) %>% do.call("bind_rows", .)
}) %>% do.call("bind_rows", .)

tissue_specific_dist_plot = best_df %>%   
  dplyr::group_by(method, Tissue) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::group_by(method) %>% 
  dplyr::mutate(f = n / sum(n)) %>% 
  dplyr::mutate(Tissue = factor(Tissue, levels = rev(c("Skeletal Muscle", "Generic", "Cerebral Cortex", "Liver", "Prostate", "Adipose Tissue")))) %>% 
  ggplot(mapping = aes(x = Tissue, y=f, fill=method)) +
  geom_col(position = "dodge") +
  theme(legend.position = "bottom") +
  coord_flip() +
  theme_bw() +
  labs(y = "Fraction", fill="") +
  scale_fill_manual(values = method_colors)
saveRDS(tissue_specific_dist_plot, "plot/final_tissue_specific_dist_plot.rds")

# Plot UMAPS
umap_glm_private = readRDS("test_analysis/glmGamPoi private.RDS")
umap_glm_and_devil = readRDS("test_analysis/glmGamPoi and devil.RDS")

umap_glm_private_plot = umap_glm_private %>% 
  ggplot(mapping = aes(x = umap_1, y = umap_2, col=age_pop)) +
  geom_point(size = .5) +
  theme_bw() +
  scale_color_manual(values = c("Old" = "goldenrod3", "Young" = "#483D8B")) +
  labs(x = "UMAP 1", y="UMAP 2", col = "Age group") +
  ggtitle("glmGamPoi private DEGs")
saveRDS(umap_glm_private_plot, "plot/final_umap_glm_private_plot.rds")


umap_glm_and_devil_plot = umap_glm_and_devil %>% 
  ggplot(mapping = aes(x = umap_1, y = umap_2, col=age_pop)) +
  geom_point(size = .5) +
  theme_bw() +
  scale_color_manual(values = c("Old" = "goldenrod3", "Young" = "#483D8B")) +
  labs(x = "UMAP 1", y="UMAP 2", col = "Age group") +
  ggtitle("glmGamPoi/devil shared DEGs")
saveRDS(umap_glm_and_devil_plot, "plot/final_umap_glm_devil_plot.rds")

# Compare gene impact
compare_gse = function(gseGO_res1, gseGO_res2, deg1, deg2, name1, name2) {
  shared_genes = intersect(deg1, deg2)
  private_1 = deg1[!deg1 %in% deg2]
  private_2 = deg2[!deg2 %in% deg1]
  
  private_paths_1 = gseGO_res1@result %>% dplyr::filter(!Description %in%  gseGO_res2@result$Description)
  private_paths_2 = gseGO_res2@result %>% dplyr::filter(!Description %in%  gseGO_res1@result$Description)
  shared_paths = gseGO_res1@result %>% dplyr::filter(Description %in% gseGO_res2@result$Description)
  
  dplyr::bind_rows(
    get_role_of_specific_genes(private_paths_1, private_1) %>% dplyr::mutate(group = name1),
    get_role_of_specific_genes(private_paths_2, private_2) %>% dplyr::mutate(group = name2),
    get_role_of_specific_genes(shared_paths, shared_genes) %>% dplyr::mutate(group = "Shared")
  ) %>% 
    ggplot(mapping = aes(x = group, y = gene_frac)) +
    geom_boxplot() +
    theme_bw() +
    labs(x = "Group", y = "Prevalence of gene group over pathway")
}

get_role_of_specific_genes = function(gseGO_res, genes_of_interest) {
  gene_role = lapply(1:nrow(gseGO_res), function(i) {
    ce = gseGO_res$core_enrichment[i]
    ce = unlist(strsplit(ce, "/"))
    sum(genes_of_interest %in% ce) / length(ce)  
  }) %>% unlist()
  
  dplyr::tibble(Description = gseGO_res$Description, gene_frac = gene_role, setSize = gseGO_res$setSize)
}

gene_group_impact_plot = compare_gse(
  deg1 = rna_deg_devil$geneID,
  deg2 = rna_deg_glm$geneID, 
  gseGO_res1 = gseGO_devil, 
  gseGO_res2 = gseGO_glm,
  name1 = "devil private", 
  name2 = "glmGamPoi private"
)
saveRDS(gene_group_impact_plot, "plot/final_gene_group_impact_plot.rds")



