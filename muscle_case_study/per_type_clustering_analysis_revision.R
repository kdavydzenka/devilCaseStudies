### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
#set.seed(1234)
source("utils/utils_analysis.R")

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "nebula" =  'steelblue', #"#B0C4DE",
  "devil" = "#099668"
)

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

c = conditions[1]

for (c in conditions) {
  res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/", c)
  dir.create(res.dir, recursive = TRUE)
  
  # Loading data #
  rna_devil_name <- paste0("results/MuscleRNA/full/devil_", c, ".RDS")
  rna_devil <- readRDS(rna_devil_name)
  
  rna_glm_name <- paste0("results/MuscleRNA/full/glmGamPoi_", c, ".RDS")
  rna_glm <- readRDS(rna_glm_name)
  
  rna_nebula_name <- paste0("results/MuscleRNA/full/nebula_", c, ".RDS")
  rna_nebula <- readRDS(rna_nebula_name) %>%
    dplyr::mutate(lfc = lfc / log(2))  # convert to log2 if needed
  
  # Cutoffs (used as hyperparameters in Bayesian-style prob)
  lfc_cut <- 1.0
  pval_cut <- 0.05
  
  # k = 8
  # # --- Bayesian-style DE probabilities --- #
  # rna_devil  <- bayes_de_prob(rna_devil,  lfc_cut = lfc_cut, alpha = pval_cut, k = k) %>%
  #   dplyr::mutate(method = "devil")
  # 
  # rna_glm    <- bayes_de_prob(rna_glm,    lfc_cut = lfc_cut, alpha = pval_cut, k = k) %>%
  #   dplyr::mutate(method = "glmGamPoi")
  # 
  # rna_nebula <- bayes_de_prob(rna_nebula, lfc_cut = lfc_cut, alpha = pval_cut, k = k) %>%
  #   dplyr::mutate(method = "nebula")
  # 
  # thr <- 0.5
  # rna_deg_devil <- rna_devil  %>% dplyr::filter(p_de >= thr)
  # rna_deg_glm   <- rna_glm    %>% dplyr::filter(p_de >= thr)
  # rna_deg_nebula <- rna_nebula %>% dplyr::filter(p_de >= thr)
  
  rna_deg_devil <- rna_devil %>%
    dplyr::filter(adj_pval <= pval_cut, abs(lfc) >= lfc_cut) %>%
    dplyr::mutate(method = "devil")

  rna_deg_glm <- rna_glm %>%
    dplyr::filter(adj_pval <= pval_cut, abs(lfc) >= lfc_cut) %>%
    dplyr::mutate(method = "glmGamPoi")

  rna_deg_nebula <- rna_nebula %>%
    dplyr::filter(adj_pval <= pval_cut, abs(lfc) >= lfc_cut) %>%
    dplyr::mutate(method = "nebula")
  
  rna_deg_devil$adj_pval[rna_deg_devil$adj_pval == 0] <- min(rna_deg_devil$adj_pval[rna_deg_devil$adj_pval != 0])
  rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval == 0] <- min(rna_deg_nebula$adj_pval[rna_deg_nebula$adj_pval != 0])
  rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])
  
  # Venn diagram ####
  x <- list(
    devil = rna_deg_devil$name,
    glmGamPoi = rna_deg_glm$name,
    nebula = rna_deg_nebula$name
  )
  
  venn_plot <- ggVennDiagram::ggVennDiagram(x, color = 1, lwd = 0.8) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
  venn_plot
  saveRDS(venn_plot, file.path(res.dir, "venn_plot.RDS"))
  
  # Volcano Plot ####
  rna_join <- rbind(rna_deg_devil, rna_deg_glm, rna_deg_nebula)
  rna_join <- rna_join %>%
    dplyr::mutate(
      isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
      DEtype = if_else(!isDE, "n.s.", if_else(lfc > 0, "Up-reg", "Down-reg")))
  
  de_colors <- c("Down-reg" = "steelblue", "Up-reg" = "indianred", "n.s." = "grey")
  
  p_volcanos <- rna_join %>%
    ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
    geom_point(aes(col = DEtype), size = 2.0, alpha = 0.2) +
    scale_color_manual(values = de_colors) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)),
                                    ceiling(max(rna_join$lfc)), by = 2)) +
    facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
    labs(x = expression(Log[2] ~ FC),
         y = expression(-log[10] ~ Pvalue),
         col = "DE type") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  p_volcanos
  saveRDS(p_volcanos, file.path(res.dir, "volcanos.RDS"))
  
  # (GSEA) Gene set Enrichement analysis ####
  ## Actual GSEA analysis ####
  gseGO_devil <- runGO(rna_devil, padj_col = "adj_pval", lfc_col = "lfc", method = "GSE")
  gseGO_glm <- runGO(rna_glm, padj_col = "adj_pval", lfc_col = "lfc", method = "GSE")
  gseGO_nebula <- runGO(rna_nebula, padj_col = "adj_pval", lfc_col = "lfc", method = "GSE")
  
  saveRDS(gseGO_devil, file.path(res.dir, "gseGO_devil.RDS"))
  saveRDS(gseGO_glm, file.path(res.dir, "gseGO_glm.RDS"))
  saveRDS(gseGO_nebula, file.path(res.dir, "gseGO_nebula.RDS"))
  
  ## Simplified GSEA analysis ####
  df_simp_gsea = get_simplified_GOterms(gseGO_devil, gseGO_glm, gseGO_nebula)
  
  auc_scores = lapply(unique(df_simp_gsea$model), function(m) {
    d = df_simp_gsea %>% dplyr::filter(model == m)
    dplyr::tibble(model = m, AUC = auc(d$c, d$f))
  }) %>% do.call("bind_rows", .)
  
  simp_plot <- df_simp_gsea %>%
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
  
  saveRDS(df_simp_gsea, file.path(res.dir, "gsea_df_simp.RDS"))
  saveRDS(auc_scores, file.path(res.dir, "gsea_auc_scores.RDS"))
  saveRDS(simp_plot, file.path(res.dir, "gsea_simp_plot.RDS"))
  
  # Get simplified GSEA
  cutoff = .5
  gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=cutoff)
  gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=cutoff)
  gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=cutoff)
  
  
  ## GSEA GO plot ####
  GO_list_df = dplyr::bind_rows(
    gseGO_devil@result %>% dplyr::mutate(method = "devil"),
    gseGO_glm@result %>% dplyr::mutate(method = "glmGamPoi"),
    gseGO_nebula@result %>% dplyr::mutate(method = "NEBULA")
  ) %>% 
    dplyr::group_by(method) %>% 
    # dplyr::arrange(p.adjust) %>% 
    # dplyr::slice_head(n = 10) %>% 
    dplyr::mutate(regulation = ifelse(enrichmentScore <= 0, "Down-regulated", "Up-regulated")) %>%
    dplyr::select(Description, method, p.adjust, enrichmentScore, regulation, setSize)
  
  GO_list_df_simp = dplyr::bind_rows(
    gseGO_devil_s@result %>% dplyr::mutate(method = "devil"),
    gseGO_glm_s@result %>% dplyr::mutate(method = "glmGamPoi"),
    gseGO_nebula_s@result %>% dplyr::mutate(method = "NEBULA")
  ) %>% 
    dplyr::group_by(method) %>% 
    # dplyr::arrange(p.adjust) %>% 
    # dplyr::slice_head(n = 10) %>% 
    dplyr::mutate(regulation = ifelse(enrichmentScore <= 0, "Down-regulated", "Up-regulated")) %>%
    dplyr::select(Description, method, p.adjust, enrichmentScore, regulation, setSize) %>% print(n = 100)
  
  # GO_plot = plot_dotplot_GO(devil_res = gseGO_devil@result, 
  #                           glm_res = gseGO_glm@result, 
  #                           nebula_res = gseGO_nebula@result, 
  #                           top_K = 10)
  # GO_plot
  #saveRDS(GO_plot, file.path(res.dir, "gsea_GO_plot.RDS"))
  saveRDS(GO_list_df, file.path(res.dir, "gsea_GO_list_df.RDS"))
  saveRDS(GO_list_df_simp, file.path(res.dir, "gsea_GO_list_df_simp.RDS"))
  
  ## TissueEnrich analysis ####
  library(TissueEnrich)
  tissue = target_tissue <- "Skeletal Muscle"
  tissue_gse_devil = get_tissue_specific_res(gseGO_devil, "devil")
  tissue_gse_glm = get_tissue_specific_res(gseGO_glm, "glmGamPoi")
  tissue_gse_nebula = get_tissue_specific_res(gseGO_nebula, "NEBULA")
  
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
    # na.omit() %>% 
    # dplyr::mutate(Tissue = factor(Tissue, levels = rev(c("Skeletal Muscle", "Generic", "Cerebral Cortex", "Liver", "Prostate", "Adipose Tissue")))) %>% 
    ggplot(mapping = aes(x = Tissue, y=f, fill=method)) +
    geom_col(position = "dodge") +
    theme(legend.position = "bottom") +
    coord_flip() +
    theme_bw() +
    labs(y = "Fraction", fill="") +
    scale_fill_manual(values = method_colors)
  
  saveRDS(tissue_gse_devil, file.path(res.dir, "tissue_gse_devil.RDS"))
  saveRDS(tissue_gse_glm, file.path(res.dir, "tissue_gse_glm.RDS"))
  saveRDS(tissue_gse_nebula, file.path(res.dir, "tissue_gse_nebula.RDS"))
  saveRDS(tissue_specific_dist_plot, file.path(res.dir, "gsea_tissue_plot.RDS"))
  
  # (ENRICH) ORA analysis ####
  ## Actual ORA analysis ####
  ORA_devil <- runGO(rna_devil, padj_col = "adj_pval", lfc_col = "lfc", method = "ENRICH", pval_cut = .05, lfc_cut = 1)
  ORA_glm <- runGO(rna_glm, padj_col = "adj_pval", lfc_col = "lfc", method = "ENRICH", pval_cut = .05, lfc_cut = 1)
  ORA_nebula <- runGO(rna_nebula, padj_col = "adj_pval", lfc_col = "lfc", method = "ENRICH", pval_cut = .05, lfc_cut = 1)
  
  ORA_devil_df = dplyr::bind_rows(ORA_devil$up_res_df, ORA_devil$down_res_df) %>% dplyr::mutate(method = "devil")
  ORA_glm_df = dplyr::bind_rows(ORA_glm$up_res_df, ORA_glm$down_res_df) %>% dplyr::mutate(method = "glmGamPoi")
  ORA_nebula_df = dplyr::bind_rows(ORA_nebula$up_res_df, ORA_nebula$down_res_df) %>% dplyr::mutate(method = "NEBULA")
  ORA_all = dplyr::bind_rows(ORA_devil_df, ORA_glm_df, ORA_nebula_df)
  
  ORA_list_df = ORA_all %>% 
    tidyr::separate(GeneRatio, c("n", "N"), sep = "/") %>% 
    dplyr::mutate(geneRatio = as.numeric(n) / as.numeric(N)) %>% 
    # dplyr::arrange(p.adjust) %>% 
    # dplyr::group_by(method) %>% 
    # dplyr::slice_head(n = 10) %>% 
    dplyr::select(Description, method, FoldEnrichment, p.adjust, regulation)
  
  # ORA_devil_s = list(
  #   up_res_df = clusterProfiler::simplify(ORA_devil$up_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Up-regulated"), 
  #   down_res_df = clusterProfiler::simplify(ORA_devil$down_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Down-regulated")
  # )
  # ORA_glm_s = list(
  #   up_res_df = clusterProfiler::simplify(ORA_glm$up_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Up-regulated"), 
  #   down_res_df = clusterProfiler::simplify(ORA_glm$down_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Down-regulated")
  # )
  # ORA_nebula_s = list(
  #   up_res_df = clusterProfiler::simplify(ORA_nebula$up_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Up-regulated"), 
  #   down_res_df = clusterProfiler::simplify(ORA_nebula$down_res, c=c)@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Down-regulated")
  # )
  # 
  # ORA_devil_df_s = dplyr::bind_rows(ORA_devil_s$up_res_df, ORA_devil_s$down_res_df) %>% dplyr::mutate(method = "devil")
  # ORA_glm_df_s = dplyr::bind_rows(ORA_glm_s$up_res_df, ORA_glm_s$down_res_df) %>% dplyr::mutate(method = "glmGamPoi")
  # ORA_nebula_df_s = dplyr::bind_rows(ORA_nebula_s$up_res_df, ORA_nebula_s$down_res_df) %>% dplyr::mutate(method = "NEBULA")
  # ORA_all_s = dplyr::bind_rows(ORA_devil_df_s, ORA_glm_df_s, ORA_nebula_df_s)
  # 
  # ORA_list_df_simp = ORA_all_s %>% 
  #   tidyr::separate(GeneRatio, c("n", "N"), sep = "/") %>% 
  #   dplyr::mutate(geneRatio = as.numeric(n) / as.numeric(N)) %>% 
  #   # dplyr::arrange(p.adjust) %>% 
  #   # dplyr::group_by(method) %>% 
  #   # dplyr::slice_head(n = 10) %>% 
  #   dplyr::select(Description, method, FoldEnrichment, p.adjust, regulation)
  
  # ORA_plot = ORA_all %>% 
  #   tidyr::separate(GeneRatio, c("n", "N"), sep = "/") %>% 
  #   dplyr::mutate(geneRatio = as.numeric(n) / as.numeric(N)) %>% 
  #   dplyr::arrange(p.adjust) %>% 
  #   dplyr::group_by(method) %>% 
  #   dplyr::slice_head(n = 10) %>% 
  #   ggplot(mapping = aes(x = method, y = Description, col = p.adjust, size = as.numeric(n))) +
  #   geom_point() +
  #   facet_wrap(~regulation) +
  #   theme_bw() +
  #   labs(x = "Method", y = "GO term", col = "Adjusted p-value", size = "Gene count") +
  #   scale_color_gradient(high = "darkorange", low = "mediumpurple")
  # ORA_plot
  
  saveRDS(ORA_devil, file.path(res.dir, "ORA_devil.RDS"))
  saveRDS(ORA_glm, file.path(res.dir, "ORA_glm.RDS"))
  saveRDS(ORA_nebula, file.path(res.dir, "ORA_nebula.RDS"))
  # saveRDS(ORA_plot, file.path(res.dir, "ORA_plot.RDS"))
  saveRDS(ORA_list_df, file.path(res.dir, "ORA_list_df.RDS"))
  
  # ## Simplified ORA analysis ####
  # df_simp_ora = get_simplified_GOterms(ORA_devil, ORA_glm, ORA_nebula)
  # 
  # auc_scores = lapply(unique(df_simp_ora$model), function(m) {
  #   d = df_simp_ora %>% dplyr::filter(model == m)
  #   dplyr::tibble(model = m, AUC = auc(d$c, d$f))
  # }) %>% do.call("bind_rows", .)
  # 
  # simp_plot <- df_simp_ora %>%
  #   dplyr::select(model, n_simplified, f, c) %>%
  #   tidyr::pivot_longer(c(f, n_simplified)) %>%
  #   dplyr::mutate(name = ifelse(name=="f", "Fraction simplified", "N simplified")) %>%
  #   dplyr::filter(name == "Fraction simplified") %>%
  #   ggplot(mapping = aes(x=c, y=value, col=model)) +
  #   geom_point() +
  #   geom_line() +
  #   theme_bw() +
  #   facet_wrap(~name, scales = "free", ncol = 1,strip.position = "top") +
  #   scale_color_manual(values = method_colors) +
  #   labs(y = "Value", x="Clustering cutoff", col="")
  # simp_plot
  # 
  # saveRDS(df_simp_ora, file.path(res.dir, "ORA_df_simp.RDS"))
  # saveRDS(auc_scores, file.path(res.dir, "ORA_auc_scores.RDS"))
  # saveRDS(simp_plot, file.path(res.dir, "ORA_simp_plot.RDS"))
  # ## ORA GO plot ####
  # s_cutoff = 0.5
  # ORA_devil_s = clusterProfiler::simplify(ORA_devil, cutoff=s_cutoff)
  # ORA_glm_s = clusterProfiler::simplify(ORA_glm, cutoff=s_cutoff)
  # ORA_nebula_s = clusterProfiler::simplify(ORA_nebula, cutoff=s_cutoff)
  # 
  # GO_plot = plot_dotplot_GO(ORA_devil@result, ORA_glm@result, ORA_nebula@result)
  # GO_plot
  # saveRDS(GO_plot, file.path(res.dir, "gsea_GO_plot.RDS"))
  
  # ## TissueEnrich analysis ####
  library(TissueEnrich)
  tissue = target_tissue <- "Skeletal Muscle"
  
  ## Filter ORA results by adjusted p-value
  ORA_devil$up_res@result   <- ORA_devil$up_res@result[ORA_devil$up_res@result$p.adjust   <= pval_cut, ]
  ORA_devil$down_res@result <- ORA_devil$down_res@result[ORA_devil$down_res@result$p.adjust <= pval_cut, ]
  
  ORA_glm$up_res@result     <- ORA_glm$up_res@result[ORA_glm$up_res@result$p.adjust       <= pval_cut, ]
  ORA_glm$down_res@result   <- ORA_glm$down_res@result[ORA_glm$down_res@result$p.adjust   <= pval_cut, ]
  
  ORA_nebula$up_res@result  <- ORA_nebula$up_res@result[ORA_nebula$up_res@result$p.adjust <= pval_cut, ]
  ORA_nebula$down_res@result<- ORA_nebula$down_res@result[ORA_nebula$down_res@result$p.adjust <= pval_cut, ]
  
  # devil
  tissue_ora_devil_up   <- get_tissue_specific_res(ORA_devil$up_res,   "devil")
  tissue_ora_devil_down <- get_tissue_specific_res(ORA_devil$down_res, "devil")
  tissue_ora_devil      <- dplyr::bind_rows(tissue_ora_devil_up, tissue_ora_devil_down)
  
  ## glmGamPoi
  tissue_ora_glm_up     <- get_tissue_specific_res(ORA_glm$up_res,   "glmGamPoi")
  tissue_ora_glm_down   <- get_tissue_specific_res(ORA_glm$down_res, "glmGamPoi")
  tissue_ora_glm        <- dplyr::bind_rows(tissue_ora_glm_up, tissue_ora_glm_down)
  
  ## NEBULA
  tissue_ora_nebula_up   <- get_tissue_specific_res(ORA_nebula$up_res,   "NEBULA")
  tissue_ora_nebula_down <- get_tissue_specific_res(ORA_nebula$down_res, "NEBULA")
  tissue_ora_nebula      <- dplyr::bind_rows(tissue_ora_nebula_up, tissue_ora_nebula_down)

  best_df = lapply(list(tissue_ora_devil, tissue_ora_glm, tissue_ora_nebula), function(tissue_gse) {
    print(unique(tissue_gse$method))
    path = "response to other organism"
    d = lapply(unique(tissue_gse$path_name), function(path) {
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
    # dplyr::mutate(Tissue = factor(Tissue, levels = rev(c("Skeletal Muscle", "Generic", "Cerebral Cortex", "Liver", "Prostate", "Adipose Tissue")))) %>%
    ggplot(mapping = aes(x = Tissue, y=f, fill=method)) +
    geom_col(position = "dodge") +
    theme(legend.position = "bottom") +
    coord_flip() +
    theme_bw() +
    labs(y = "Fraction", fill="") +
    scale_fill_manual(values = method_colors)

  saveRDS(tissue_ora_devil, file.path(res.dir, "tissue_ora_devil.RDS"))
  saveRDS(tissue_ora_glm, file.path(res.dir, "tissue_ora_glm.RDS"))
  saveRDS(tissue_ora_nebula, file.path(res.dir, "tissue_ora_nebula.RDS"))
  saveRDS(tissue_specific_dist_plot, file.path(res.dir, "ora_tissue_plot.RDS"))
  
  # UMAP analysis ####
  source("utils/utils.R")
  dataset_name <- "MuscleRNA"
  data_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/muscle/rna/seurat_counts_rna.RDS"
  input_data <- read_data(dataset_name, data_path)
  input_data <- prepare_rna_input(input_data)
  
  # Gene selection based on LFC & pvalue cutoff #
  lfc_strict <- 1.1
  
  nebula_genes <- rna_deg_nebula$name
  devil_genes  <- rna_deg_devil$name
  glm_genes    <- rna_deg_glm$name
  
  nebula_private_raw <- setdiff(nebula_genes, union(devil_genes, glm_genes))
  devil_private_raw  <- setdiff(devil_genes,  union(nebula_genes, glm_genes))
  glm_private_raw    <- setdiff(glm_genes,    union(nebula_genes, devil_genes))
  glm_devil_shared_raw <- intersect(glm_genes, devil_genes)
  all_shared_raw <- Reduce(intersect, list(nebula_genes, devil_genes, glm_genes))
  
  devil_private <- rna_devil %>%
    dplyr::filter(name %in% devil_private_raw,
                  abs(lfc) > lfc_strict) %>%
    dplyr::pull(name)
  
  glm_private <- rna_glm %>%
    dplyr::filter(name %in% glm_private_raw,
                  abs(lfc) > lfc_strict) %>%
    dplyr::pull(name)
  
  nebula_private <- rna_nebula %>%
    dplyr::filter(name %in% nebula_private_raw,
                  abs(lfc) > lfc_strict) %>%
    dplyr::pull(name)
  
  glm_devil_shared <- glm_devil_shared_raw
  all_shared       <- all_shared_raw
  
  gene_sets <- list(
    nebula_private = nebula_private,
    devil_private  = devil_private,
    glm_private    = glm_private,
    all_shared     = all_shared
  )
  
  gene_list = list(
    "glmGamPoi private" = glm_private, 
    "devil private" = devil_private,
    "nebula private" = nebula_private,
    "glmGamPoi and devil" = glm_devil_shared,
    "shared" = all_shared
  )
  
  N_subsample <- 10000
  N_gene_sub = 100
  sample_idx = sample(1:ncol(input_data$counts), N_subsample, replace = FALSE)
  umap_plots = list()
  df_acc = dplyr::tibble()
  for (i in 1:length(gene_list)) {
    list_name = names(gene_list)[i]
    intersting_genes = gene_list[[i]]
    
    if(length(intersting_genes) < 50) next
    intersting_genes <- sample(intersting_genes, size = min(N_gene_sub, length(intersting_genes)))
    #intersting_genes = intersting_genes[1:min(N_gene_sub, length(intersting_genes))]
    
    # mat <- input_data$counts[intersting_genes,sample_idx] %>% as.matrix()
    mat <- input_data$counts[,sample_idx] %>% as.matrix()
    meta <- input_data$metadata[sample_idx,] 
    
    meta <- meta %>% 
      dplyr::mutate(age_pop = ifelse(age_pop == "old_pop", "Old", "Young"))
    
    # reorder by samples
    # ordered_indices <- order(meta$sample)
    # mat <- mat[, ordered_indices]
    # meta <- meta[ordered_indices, ]
    
    design_test <- c
    if (design_test == "age_only") {
      
      # Just old vs young (all cell types together)
      meta <- meta %>%
        dplyr::mutate(group = age_pop)   # "Young", "Old"
      
      design_name = "Age across cell types"
      
    } else if (design_test == "age_type1") {
      
      # Age effect within cell_type1; everything else can be "Other"
      meta <- meta %>%
        dplyr::mutate(
          group = dplyr::case_when(
            cell_type == "Type I" & age_pop == "Young" ~ "Young - Type I",
            cell_type == "Type I" & age_pop == "Old"   ~ "Old - Type I",
            TRUE                                      ~ "Other"
          )
        )
      
      design_name = "Age in Type I cells"
      
    } else if (design_test == "age_type2") {
      
      # Age effect within cell_type2; everything else "Other"
      meta <- meta %>%
        dplyr::mutate(
          group = dplyr::case_when(
            cell_type == "Type II" & age_pop == "Young" ~ "Young - Type II",
            cell_type == "Type II" & age_pop == "Old"   ~ "Old - Type II",
            TRUE                                      ~ "Other"
          )
        )
      
      design_name = "Age in Type II cells"
      
    } else if (design_test == "interaction") {
      # Full 2x2: age × cell_type
      meta <- meta %>%
        dplyr::mutate(
          group = dplyr::case_when(
            cell_type == "Type I" & age_pop == "Young" ~ "Young - Type I",
            cell_type == "Type I" & age_pop == "Old"   ~ "Old - Type I",
            cell_type == "Type II" & age_pop == "Young" ~ "Young Type II",
            cell_type == "Type II" & age_pop == "Old"   ~ "Old Type II",
            TRUE                                      ~ "Other"
          )
        )
      
      design_name = "Iteraction term"
    }
    
    # subset out "Other" *before* constructing Seurat object
    keep <- meta$group != "Other"
    mat  <- mat[, keep]
    meta <- meta[keep, ]
    
    seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta)
    seurat_obj = subset(seurat_obj, features = intersting_genes)
    
    seurat_obj <- NormalizeData(seurat_obj)
    #seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = intersting_genes)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
    
    df_umap = seurat_obj@reductions$umap@cell.embeddings %>% as.tibble()
    df_umap$group = seurat_obj$group
    df_umap$cell_type = seurat_obj$cell_type
    df_umap$age_cluster = seurat_obj$age_cluster
    df_umap$sample = seurat_obj$sample
    
    # df_umap %>% 
    #   ggplot(mapping = aes(x = umap_1, y = umap_2, colour = sample)) +
    #   geom_point(size = .5) +
    #   theme_bw() + 
    #   ggtitle(list_name, design_name) +
    #   labs(x = "UMAP 1", y = "UMAP 2", col = "Group")
    
    p = df_umap %>% 
      ggplot(mapping = aes(x = umap_1, y = umap_2, colour = group)) +
      geom_point(size = .5) +
      theme_bw() + 
      ggtitle(list_name, design_name) +
      labs(x = "UMAP 1", y = "UMAP 2", col = "Group")
    
    print(p)
    umap_plots[[i]] = p
    
    # # Accuracy
    library(caret)
    library(glmnet)
    # Prepare data

    scaled_mat <- seurat_obj@assays$RNA$scale.data
    X = seurat_obj@reductions$umap@cell.embeddings
    X = seurat_obj@reductions$pca@cell.embeddings[,1:2]
    y <- factor(seurat_obj$group)

    set.seed(123)
    cv_model <- cv.glmnet(
      X, y,
      family       = "multinomial",
      alpha        = 0,
      nfolds       = 5,
      type.measure = "class"
    )

    pred <- predict(cv_model, X, s = "lambda.min", type = "class")
    pred <- as.factor(pred[, 1])  # ensure it's a 1D factor
    
    cm  <- confusionMatrix(pred, y)
    acc <- cm$overall["Accuracy"]

    df_acc = dplyr::bind_rows(df_acc, dplyr::tibble(ACC = acc, list = list_name))
  }
  
  names(umap_plots) = names(gene_list)
  saveRDS(umap_plots, file.path(res.dir, "umaps.RDS"))
  saveRDS(df_acc, file.path(res.dir, "df_acc.RDS"))
  
  # Subsampling analysis ####
  
  sub_metrics = dplyr::tibble()
  sub_df = dplyr::tibble()
  df_jacc = dplyr::tibble()
  
  m = "devil"
  for (m in c("devil", "glmGamPoi", "nebula")) {
    sanitize_p <- function(p) {
      # Replace NA with 1
      p[is.na(p)] <- 1

      # Replace zeros with tiny value
      p[p == 0] <- 1e-300

      # Clip to valid range
      p[p < 1e-300] <- 1e-300
      p[p > 1] <- 1
      
      return(p)
    }
    
    compute_rankmetric <- function(df) {
      df %>%
        mutate(
          pval = sanitize_p(pval),
          #RankMetric = -log(pval) * lfc
          RankMetric = -log(pval) * sign(lfc)
          #RankMetric = -log(pval) * abs(lfc)
          # RankMetric = -log(pval)
        ) %>%
        arrange(desc(RankMetric))
    }
    
    rna_full <- readRDS(paste0("results/MuscleRNA/full/", m, "_", c, ".RDS")) %>%
      compute_rankmetric()
    
    rna_sub <- readRDS(paste0("results/MuscleRNA/subsampled/", m, "_", c, ".RDS")) %>%
      compute_rankmetric()
    
    if ("geneID" %in% colnames(rna_full)) rna_full$name = rna_full$geneID
    if ("geneID" %in% colnames(rna_sub)) rna_sub$name = rna_sub$geneID
    
    library(dplyr)
    # library(RankAggreg)
    # library(rbo)
    
    # Merge full + sub on gene
    df = rna_full %>%
      dplyr::select(name, RankMetric, lfc) %>%
      dplyr::mutate(rank_full = rank(-RankMetric, ties.method = "average")) %>%
      inner_join(
        rna_sub %>% 
          dplyr::select(name, RankMetric, lfc) %>% 
          dplyr::mutate(rank_sub = rank(-RankMetric, ties.method = "average")),
        by = "name"
      )
    
    # Spearman + Kendall
    spearman = cor(df$rank_full, df$rank_sub, method = "spearman")
    kendall  = cor(df$rank_full, df$rank_sub, method = "kendall")
    
    # Rank displacement
    df = df %>% dplyr::mutate(rank_shift = rank_sub - rank_full)
    mean_abs_shift = mean(abs(df$rank_shift))
    median_abs_shift = median(abs(df$rank_shift))
    
    # Top-k overlap
    j_vec = lapply(1:2000, function(k) {
      top_full = df$name[order(df$rank_full)][1:k]
      top_sub  = df$name[order(df$rank_sub)][1:k]
      jaccard_topk = length(intersect(top_full, top_sub)) /
        length(union(top_full, top_sub))  
      jaccard_topk
    }) %>% unlist()
    # plot(1:length(j_vec), j_vec)
    df_jacc = dplyr::bind_rows(
      df_jacc,
      dplyr::tibble(N = 1:length(j_vec), J = j_vec, method = m)
    )
    
    k = 50
    top_full = df$name[order(df$rank_full)][1:k]
    top_sub  = df$name[order(df$rank_sub)][1:k]
    jaccard_topk = length(intersect(top_full, top_sub)) /
      length(union(top_full, top_sub))
    
    sub_df = dplyr::bind_rows(
      sub_df,
      df %>% 
        dplyr::mutate(method = m) %>% 
        dplyr::mutate(type = c))
    
    sub_metrics = dplyr::bind_rows(
      sub_metrics,
      dplyr::tibble(method = m, type = c, spearman = spearman, kendall = kendall, jaccard_topk = jaccard_topk, mean_abs_shift, median_abs_shift)
    )
  }
  
  rank_plot_corr = sub_df %>% 
    ggplot(aes(rank_full, rank_sub, col = method)) +
    geom_point(alpha = .3) +
    geom_abline(color = "black") +
    theme_bw() +
    labs(title = "Rank stability: Full vs Subsample") +
    facet_wrap(~method) +
    scale_color_manual(values = method_colors)
  
  rank_shift_plot = sub_df %>% 
    ggplot(aes(x = rank_shift, col = method)) +
    geom_density() +
    theme_bw() +
    labs(x = "Rank shift", y = "Density") +
    scale_color_manual(values = method_colors)
  
  plot_jacc = df_jacc %>% 
    ggplot(mapping = aes(x = N, y = J, col = method)) +
    geom_line() +
    theme_bw() +
    labs(x = "Top K", y = "Jaccard Index", col = "Method") +
    scale_color_manual(values = method_colors)
  
  saveRDS(sub_df, file.path(res.dir, "sub_df.RDS"))
  saveRDS(sub_metrics, file.path(res.dir, "sub_metrics_df.RDS"))
  saveRDS(rank_plot_corr, file.path(res.dir, "sub_rank_plot_corr.RDS"))
  saveRDS(rank_shift_plot, file.path(res.dir, "sub_rank_shift_plot.RDS"))
  saveRDS(plot_jacc, file.path(res.dir, "plot_jacc.RDS"))
}

print("All done!")
stop()

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
  deg1 = rna_deg_devil$name,
  deg2 = rna_deg_glm$name, 
  gseGO_res1 = gseGO_devil, 
  gseGO_res2 = gseGO_glm,
  name1 = "devil private", 
  name2 = "glmGamPoi private"
)







# Classification accuracy ####
i = 5
list_name = names(gene_list)[i]
intersting_genes = gene_list[[i]]

# 1. Subset expression and metadata
mat <- input_data$counts[intersting_genes, sample_idx] %>% as.matrix()
meta <- input_data$metadata[sample_idx,] 

meta <- meta %>% 
  dplyr::mutate(age_pop = ifelse(age_pop == "old_pop", "Old", "Young"))

design_test <- c
meta$cell_type

if (design_test == "age_only") {
  
  # Just old vs young (all cell types together)
  meta <- meta %>%
    dplyr::mutate(group = age_pop)   # "Young", "Old"
  
  design_name = "Age across cell types"
  
} else if (design_test == "age_type1") {
  
  # Age effect within cell_type1; everything else can be "Other"
  meta <- meta %>%
    dplyr::mutate(
      group = dplyr::case_when(
        cell_type == "Type I" & age_pop == "Young" ~ "Young - Type I",
        cell_type == "Type I" & age_pop == "Old"   ~ "Old - Type I",
        TRUE                                      ~ "Other"
      )
    )
  
  design_name = "Age in Type I cells"
  
} else if (design_test == "age_type2") {
  
  # Age effect within cell_type2; everything else "Other"
  meta <- meta %>%
    dplyr::mutate(
      group = dplyr::case_when(
        cell_type == "Type II" & age_pop == "Young" ~ "Young - Type II",
        cell_type == "Type II" & age_pop == "Old"   ~ "Old - Type II",
        TRUE                                      ~ "Other"
      )
    )
  
  design_name = "Age in Type II cells"
  
} else if (design_test == "interaction") {
  # Full 2x2: age × cell_type
  meta <- meta %>%
    dplyr::mutate(
      group = dplyr::case_when(
        cell_type == "Type I" & age_pop == "Young" ~ "Young - Type I",
        cell_type == "Type I" & age_pop == "Old"   ~ "Old - Type I",
        cell_type == "Type II" & age_pop == "Young" ~ "Young Type II",
        cell_type == "Type II" & age_pop == "Old"   ~ "Old Type II",
        TRUE                                      ~ "Other"
      )
    )
  
  design_name = "Iteraction term"
}

# -------------------------------
# Compute classification accuracy
# -------------------------------

library(caret)
library(glmnet)

keep <- meta$group != "Other"
mat  <- mat[, keep]
meta <- meta[keep, ]

# Prepare data
X <- t(mat[intersting_genes , keep])
y <- factor(meta$group)

# Cross-validated glmnet (very fast, handles high dimensions)
set.seed(123)
cv_model <- cv.glmnet(X, y, 
                      family = "multinomial",
                      alpha = 1,  # 1 = lasso, 0 = ridge, 0.5 = elastic net
                      nfolds = 5,
                      type.measure = "class")

# Predict and evaluate
predictions <- predict(cv_model, X, s = "lambda.min", type = "class")
confusionMatrix(factor(predictions), y)

# AUC
auc_score <- max(cv_model$cvm)
print(paste("AUC:", round(auc_score, 3)))
