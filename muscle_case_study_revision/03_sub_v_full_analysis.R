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

methods <- c("devil", "glmGamPoi", "NEBULA")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

c = conditions[2]

for (c in conditions) {
  res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/", c)
  dir.create(res.dir, recursive = TRUE)
  
  sub_metrics = dplyr::tibble()
  sub_df = dplyr::tibble()
  df_jacc = dplyr::tibble()
  
  m = "devil"
  n_k = sum(readRDS(paste0("results/MuscleRNA/full/", m, "_", c, ".RDS"))$pval == 0)
  
  for (m in c("devil", "glmGamPoi", "nebula")) {
    sanitize_p <- function(p) {
      # Replace NA with 1
      p[is.na(p)] <- 1
      
      # Replace zeros with tiny value
      p[p == 0] <- min(p[p != 0])
      
      # Clip to valid range
      #p[p < 1e-300] <- 1e-300
      # p[p == 0] <- 1e-300 * runif(sum(p == 0), 0.9, 1.1)
      p[p > 1] <- 1
      
      return(p)
    }
    
    compute_rankmetric <- function(df) {
      df %>%
        mutate(
          pval = sanitize_p(pval),
          #RankMetric = -log(pval) * lfc
          #RankMetric = -log(pval) * sign(lfc)
          RankMetric = -log(pval) * lfc
          #RankMetric = -log(pval)
          #RankMetric = lfc
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
    # df = rna_full %>%
    #   dplyr::select(name, RankMetric, lfc) %>%
    #   dplyr::mutate(rank_full = rank(-RankMetric, ties.method = "average")) %>%
    #   inner_join(
    #     rna_sub %>% 
    #       dplyr::select(name, RankMetric, lfc) %>% 
    #       dplyr::mutate(rank_sub = rank(-RankMetric, ties.method = "average")),
    #     by = "name"
    #   )
    df <- rna_full %>%
      dplyr::select(name, RankMetric, lfc) %>%
      dplyr::arrange(desc(RankMetric), desc(abs(lfc))) %>%
      dplyr::mutate(rank_full = row_number()) %>%
      dplyr::inner_join(
        rna_sub %>%
          dplyr::select(name, RankMetric, lfc) %>%
          dplyr::arrange(desc(RankMetric), desc(abs(lfc))) %>%
          dplyr::mutate(rank_sub = row_number()),
        by = "name"
      )
    
    # Spearman + Kendall
    spearman = cor(df$rank_full, df$rank_sub, method = "spearman")
    kendall  = cor(df$rank_full, df$rank_sub, method = "kendall")
    
    rna_full %>% 
      dplyr::arrange(-RankMetric)
    
    # Rank displacement
    df = df %>% dplyr::mutate(rank_shift = rank_sub - rank_full)
    mean_abs_shift = mean(abs(df$rank_shift))
    median_abs_shift = median(abs(df$rank_shift))
    
    # Top-k overlap
    ks = c(1:999, seq(1000, 10000, by = 50))
    j_vec = lapply(ks, function(k) {
      top_full = df$name[order(df$rank_full)][1:k]
      top_sub  = df$name[order(df$rank_sub)][1:k]
      jaccard_topk = length(intersect(top_full, top_sub)) /
        length(union(top_full, top_sub))  
      jaccard_topk
    }) %>% unlist()
    # plot(1:length(j_vec), j_vec)
    df_jacc = dplyr::bind_rows(
      df_jacc,
      dplyr::tibble(N = ks, J = j_vec, method = m)
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
  
  # 
  gsea_go_list_full = readRDS(paste0("results/MuscleRNA/per_contrast_vector_analysis/full/", c, "/gsea_GO_list_df.RDS"))
  gsea_go_list_sub = readRDS(paste0("results/MuscleRNA/per_contrast_vector_analysis/subsampled/", c, "/gsea_GO_list_df.RDS"))
  
  # Create Venn plots for each method comparing full vs sub
  library(patchwork)
  
  venn_plots <- list()
  df_venn_bar = dplyr::tibble()
  m = "NEBULA"
  for (m in methods) {
    x <- list(
      full = gsea_go_list_full %>%
        dplyr::filter(method == m) %>%
        dplyr::pull(Description),
      sub = gsea_go_list_sub %>%
        dplyr::filter(method == m) %>%
        dplyr::pull(Description)
    )
    
    df_venn_bar = dplyr::bind_rows(
      df_venn_bar,
      dplyr::tibble(method = m, all = length(intersect(x$full, x$sub)), sub = sum(!x$sub %in% x$full), full = sum(!x$full %in% x$sub))
    )
    
    venn_plots[[m]] <- ggVennDiagram(x, color = 1, lwd = 0.8) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggtitle(m) +
      theme(legend.position = "none")
  }
  
  pal <- c(
    "Subsample only"  = "#BFD3E6",  # muted light blue
    "Shared"  = "#D81B60",  # vivid magenta highlight
    "Full only" = "#8b96c5"   # deep indigo-purple
  )
  
  venn_bar_plot = df_venn_bar %>% 
    tidyr::pivot_longer(!method) %>% 
    group_by(method) %>%
    mutate(prop = value / sum(value)) %>%
    mutate(
      name = factor(name, levels = c("sub","all","full"),
                    labels = c("Subsample only",
                               "Shared",
                               "Full only"))
    ) %>% 
    ungroup() %>% 
    ggplot(aes(x = method, y = prop, fill = name)) +
    geom_bar(stat = "identity") +
    #geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")  + 
    labs(x = "DE method", y = "Share of GSEA GO terms", fill = NULL) +
    coord_flip() +
    scale_fill_manual(values = pal)
  venn_bar_plot
  
  # Combine all Venn plots into a single figure
  combined_venn_gsea <- wrap_plots(venn_plots, ncol = 3)
  
  ora_go_list_full = readRDS(paste0("results/MuscleRNA/per_contrast_vector_analysis/full/", c, "/ORA_list_df.RDS"))
  ora_go_list_sub = readRDS(paste0("results/MuscleRNA/per_contrast_vector_analysis/subsampled/", c, "/ORA_list_df.RDS"))
  
  # Create Venn plots for each method comparing full vs sub
  library(patchwork)
  
  venn_plots <- list()
  
  for (m in methods) {
    x <- list(
      full = ora_go_list_full %>%
        dplyr::filter(method == m) %>%
        dplyr::pull(Description),
      sub = ora_go_list_sub %>%
        dplyr::filter(method == m) %>%
        dplyr::pull(Description)
    )
    
    venn_plots[[m]] <- ggVennDiagram(x, color = 1, lwd = 0.8) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggtitle(m) +
      theme(legend.position = "none")
  }
  
  # Combine all Venn plots into a single figure
  combined_venn_ora <- wrap_plots(venn_plots, ncol = 3)
  print(combined_venn_ora)
  
  # Create a list to store all GO terms by method and dataset type for UpSet plot
  go_sets <- list()
  
  for (m in methods) {
    # Extract GO terms for full dataset
    go_sets[[paste0(m, " (full)")]] <- gsea_go_list_full %>%
      dplyr::filter(method == m) %>%
      dplyr::pull(Description)
    
    # Extract GO terms for subsampled dataset
    go_sets[[paste0(m, " (sub)")]] <- gsea_go_list_sub %>%
      dplyr::filter(method == m) %>%
      dplyr::pull(Description)
  }
  
  # Create UpSet plot
  library(UpSetR)
  upset_plot = upset(fromList(go_sets), 
        order.by = "freq",
        nsets = length(go_sets),
        nintersects = 40,
        mb.ratio = c(0.55, 0.45),
        text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2))
  
  saveRDS(sub_df, file.path(res.dir, "sub_df.RDS"))
  saveRDS(sub_metrics, file.path(res.dir, "sub_metrics_df.RDS"))
  saveRDS(rank_plot_corr, file.path(res.dir, "sub_rank_plot_corr.RDS"))
  saveRDS(rank_shift_plot, file.path(res.dir, "sub_rank_shift_plot.RDS"))
  saveRDS(plot_jacc, file.path(res.dir, "plot_jacc.RDS"))
  saveRDS(venn_bar_plot, file.path(res.dir, "venn_bar_plot_gsea.RDS"))
  saveRDS(combined_venn_gsea, file.path(res.dir, "combined_venn_gsea.RDS"))
  saveRDS(combined_venn_ora, file.path(res.dir, "combined_venn_ora.RDS"))
  saveRDS(upset_plot, file.path(res.dir, "upset_plot.RDS"))
}
