setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "glue", "viridis", "smplot2", "gridExtra",
          "ggpubr", "ggrepel", "ggpointdensity","clusterProfiler", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils_analysis.R")

# ORA enrichment - per gene category

classified_all <- readRDS("plot/revision/data_to_plot/classified_genes.RDS")

de_genes_classified <- classified_all %>%
  dplyr::filter(category != "Not significant")

gene_lists <- de_genes_classified %>%
  group_by(method, category) %>%
  summarise(genes = list(unique(name)), .groups = "drop") %>%
  unite(label, method, category, sep = "_") %>%
  tibble::deframe()

saveRDS(gene_lists, file = "results/ORA/gene_lists_method_category.RDS")

#background <- unique(de_genes_classified$name)

run_ORA <- function(genes, label) {
  enrichGO(
    gene          = genes,
    #universe = background,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  ) %>% 
    as.data.frame() %>% 
    mutate(GeneSet = label)
}

ora_results <- lapply(names(gene_lists), function(label) {
  run_ORA(gene_lists[[label]], label)
})

names(ora_results) <- names(gene_lists)

saveRDS(ora_results, file = "results/ORA/ora_results_list.RDS")

ora_results <- readRDS("results/ORA/ora_results_list.RDS")

ora_results <- ora_results %>% bind_rows()
ora_sig <- ora_results %>% 
  dplyr::filter(p.adjust < 0.05)

sel_path <- c("muscle_system_process",
              "muscle_contraction",
              "muscle_tissue_development",
              "homophilic cell adhesion via plasma membrane adhesion molecules",
              "striated muscle contraction",
              "regulation of actin filament-based process",
              "muscle cell differentiation",
              "myofibril assembly",
              "regulation of protein transport",
              "response to oxidative stress",
              "actomyosin structure organization",
              "regulation of protein kinase activity",
              "ferroptosis",
              "transforming growth factor beta receptor superfamily signaling pathway",
              "cell surface receptor protein serine/threonine kinase signaling pathway",
              "regulation of p38MAPK cascade",
              "response to metal ion",
              "cardiac muscle tissue development",
              "regulation of cell-substrate adhesion",
              "regulation of presynaptic membrane potential",
              "actin filament-based movement",
              "gliogenesis",
              "positive regulation of angiogenesis",
              "response to reactive oxygen species",
              "myofibril assembly",
              "cell-substrate adhesion",
              "muscle system process")

res_sel <- ora_sig %>%
  dplyr::filter(Description %in% sel_path)

path_cat_ora <- tibble::tribble(
  ~Description, ~Category,
  
  # # Muscle / Structural Development
  "muscle system process", "Muscle / Structural Development",
  "muscle contraction", "Muscle / Structural Development",
  "muscle tissue development", "Muscle / Structural Development",
  "striated muscle contraction", "Muscle / Structural Development",
  "muscle cell differentiation", "Muscle / Structural Development",
  "myofibril assembly", "Muscle / Structural Development",
  "actomyosin structure organization", "Muscle / Structural Development",
  "actin filament-based movement", "Muscle / Structural Development",
  "cardiac muscle tissue development", "Muscle / Structural Development",
  
  # Cell Adhesion & Cytoskeletal Regulation
  "homophilic cell adhesion via plasma membrane adhesion molecules", "Cell Adhesion / Cytoskeleton",
  "regulation of actin filament-based process", "Cell Adhesion / Cytoskeleton",
  "regulation of cell-substrate adhesion", "Cell Adhesion / Cytoskeleton",
  "cell-substrate adhesion", "Cell Adhesion / Cytoskeleton",
  
  # Signaling & Kinase Cascades
  "regulation of protein kinase activity", "Signaling / Kinase",
  "transforming growth factor beta receptor superfamily signaling pathway", "Signaling / Kinase",
  "cell surface receptor protein serine/threonine kinase signaling pathway", "Signaling / Kinase",
  "regulation of p38MAPK cascade", "Signaling / Kinase",
  "regulation of protein transport", "Signaling / Kinase",
  
  # Oxidative Stress / Metal / Ferroptosis
  "response to oxidative stress", "Stress Response",
  "response to reactive oxygen species", "Stress Response",
  "response to metal ion", "Stress Response",
  "ferroptosis", "Stress Response",
  
  # Neuro / Angiogenesis
  "gliogenesis", "Neural Development",
  "regulation of presynaptic membrane potential", "Neural Development",
  "positive regulation of angiogenesis", "Angiogenesis"
)

res_sel <- res_sel %>%
  left_join(path_cat_ora, by = "Description") 

heatmap_data <- res_sel %>%
  dplyr::filter(Description %in% sel_path) %>%
  dplyr::select(Category, GeneSet, Description, FoldEnrichment)

heatmap_ora <- ggplot(
  heatmap_data,
  aes(x = GeneSet, y = Description, fill = FoldEnrichment)
) +
  geom_tile(color = "grey60", linewidth = 0.25) +
  
  facet_grid(
    rows = vars(Category),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0,
    name = "Enrichment\nScore",
    oob = scales::squish
  ) +
  
  labs(
    x = "Method / Test Type",
    y = "GO Biological Process"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    # Axis text
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      size = 11, color = "black"
    ),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    
    # Facet label formatting
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.text.y.left = element_text(
      angle = 0, face = "bold", size = 12, color = "black"
    ),
    
    # Panel/grid
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5
    ),
    panel.spacing.y = unit(1.0, "lines"),
    
    # Margins
    plot.margin = margin(10, 15, 10, 10),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
heatmap_ora

ggsave("plot/revision/heatmap_ora.png", dpi = 400, width = 11.0, height = 7.0, plot = heatmap_gsea)
