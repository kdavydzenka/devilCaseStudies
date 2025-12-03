### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
#set.seed(1234)
#source("utils/utils_analysis.R")
source("utils/go_and_ora_classes.R")

TOP_K = 10
method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "devil" = "#099668"
)

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

c = conditions[3]

for (c in conditions) {
  res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/", c)
  
  files <- list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
  
  for (f in files) {
    # Extract filename without extension
    obj_name <- tools::file_path_sans_ext(basename(f))
    
    # Read the RDS
    obj <- readRDS(f)
    
    # Assign to the global environment (or specify an environment)
    assign(obj_name, obj)
  }
  
  venn_plot
  volcanos
  
  gsea_simp_plot
  gsea_auc_scores
  
  gsea_GO_list_df_top = gsea_GO_list_df %>% 
    dplyr::group_by(method) %>%
    #dplyr::arrange(p.adjust) %>% 
    dplyr::arrange(-log(p.adjust) * abs(enrichmentScore)) %>%
    dplyr::slice_head(n = TOP_K)
  gsea_GO_list_df_top$class = go_classification[gsea_GO_list_df_top$Description]
  
  gsea_GO_list_df_top %>% 
    ggplot(mapping = aes(x = method, y = Description, col = p.adjust, size = setSize)) +
    geom_point() +
    facet_grid(class~regulation, scales = "free", space = "free") +
    scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p-value") +
    labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Count") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  
  ORA_list_df_top = ORA_list_df %>% 
    dplyr::group_by(method) %>%
    #dplyr::arrange(p.adjust) %>% 
    dplyr::arrange(-log(p.adjust) * abs(FoldEnrichment)) %>%
    dplyr::slice_head(n = TOP_K)
  ORA_list_df_top$class = go_classification[ORA_list_df_top$Description]
  
  ORA_list_df_top %>% 
    ggplot(mapping = aes(x = method, y = Description, col = p.adjust, size = FoldEnrichment)) +
    geom_point() +
    facet_grid(class~regulation, , scales = "free", space = "free") +
    labs(x = "Method", y = "GO term", col = "Adjusted p-value", size = "FoldEnrichment") +
    scale_color_gradient(high = "darkorange", low = "mediumpurple") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  
  sub_metrics_df
  sub_rank_plot_corr
  sub_rank_shift_plot
  
  df_acc
  
}
