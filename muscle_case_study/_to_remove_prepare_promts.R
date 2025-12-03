### Results downstream analysis ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
          "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
#set.seed(1234)
source("utils/utils_analysis.R")
source("utils/go_classification.R")

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "devil" = "#099668"
)

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

c = conditions[4]

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
  
  TOP_K = 10
  
  gsea_GO_list_df_top = gsea_GO_list_df %>% 
    dplyr::group_by(method) %>%
    #dplyr::arrange(p.adjust) %>% 
    dplyr::arrange(-log(p.adjust) * abs(enrichmentScore)) %>%
    dplyr::slice_head(n = TOP_K)
  gsea_GO_list_df_top %>% 
    dplyr::select(method, Description) %>% 
    print(n = 3 * TOP_K)
  
  ORA_list_df_top = ORA_list_df %>% 
    dplyr::group_by(method) %>%
    #dplyr::arrange(p.adjust) %>% 
    dplyr::arrange(-log(p.adjust) * abs(FoldEnrichment)) %>%
    dplyr::slice_head(n = TOP_K)
  ORA_list_df_top %>% 
    dplyr::select(method, Description) %>% 
    print(n = 3 * TOP_K)
}  