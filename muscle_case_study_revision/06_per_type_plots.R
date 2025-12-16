
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
#conditions <- c("age_only", "age_type1", "age_type2", "interaction")
conditions <- c("age_type1", "age_type2", "interaction")

c = conditions[2]

for (c in conditions) {
  res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/full/", c)
  files <- list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
  for (f in files) {
    obj_name <- tools::file_path_sans_ext(basename(f))
    obj <- readRDS(f)
    assign(obj_name, obj)
  }
  
  res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/", c)
  files <- list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
  for (f in files) {
    obj_name <- tools::file_path_sans_ext(basename(f))
    obj <- readRDS(f)
    assign(obj_name, obj)
  }
  
  venn_plot
  venn_bar_plot
  df_acc
  
  volcanos
  
  gsea_simp_plot
  gsea_auc_scores
  
  ORA_simp_plot
  ORA_auc_scores
  
  gsea_GO_list_df_top = gsea_GO_list_df_simp %>% 
    dplyr::group_by(method) %>%
    dplyr::arrange(p.adjust) %>% 
    #dplyr::arrange(-log(p.adjust) * abs) %>% 
    dplyr::slice_head(n = TOP_K)
  gsea_GO_list_df_top$class = go_classification[gsea_GO_list_df_top$Description]
  
  gsea_plot = gsea_GO_list_df_top %>% 
    ggplot(mapping = aes(x = method, y = Description, col = p.adjust, size = setSize)) +
    geom_point() +
    facet_grid(class~regulation, scales = "free", space = "free") +
    scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p-value") +
    labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Count") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  gsea_plot
  
  ORA_list_df_top = ORA_list_df %>% 
    dplyr::group_by(method) %>%
    dplyr::arrange(p.adjust) %>% 
    #dplyr::arrange(-log(p.adjust) * abs(FoldEnrichment)) %>%
    dplyr::slice_head(n = TOP_K)
  ORA_list_df_top$class = go_classification[ORA_list_df_top$Description]

  ORA_plot = ORA_list_df_top %>% 
    ggplot(mapping = aes(x = method, y = Description, col = p.adjust, size = FoldEnrichment)) +
    geom_point() +
    facet_grid(class~regulation, , scales = "free", space = "free") +
    labs(x = "Method", y = "GO term", col = "Adjusted p-value", size = "FoldEnrichment") +
    scale_color_gradient(high = "darkorange", low = "mediumpurple") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  ORA_plot
  
  sub_metrics_df
  sub_rank_shift_plot
  
  plot_jacc
  umaps$`glmGamPoi and devil`
  umaps$`glmGamPoi private`
  umaps$`nebula private`
  
  combined_venn_gsea
    
  gsea_tissue_plot
  ora_tissue_plot
  
  des = "
  AB
  CD
  EF
  EF
  EF
  GH
  IL
  MN"
  
  p = volcanos + venn_plot + 
    gsea_simp_plot + ORA_simp_plot + 
    gsea_plot + ORA_plot +
    gsea_tissue_plot + ora_tissue_plot +
    combined_venn_gsea + combined_venn_ora + 
    sub_rank_shift_plot + plot_jacc +
    plot_layout(design = des)
  
  dir.create(paste0("figures/",c))
  ggsave(paste0("figures/",c, "/report.pdf"), plot = p, width = 30, height = 30, units = "in")
}
