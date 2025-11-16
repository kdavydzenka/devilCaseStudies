### Revised downstream analysis using different design combinations ###

rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "glue", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ggVennDiagram",
          "clusterProfiler", "org.Hs.eg.db", "ComplexUpset", "TissueEnrich")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils_analysis.R")


###------------Load results---------------###

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age_only", "age_type1", "age_type2", "interaction")

read_rna <- function(method, condition, rename = TRUE) {
  path <- glue::glue("results/MuscleRNA/full/{method}_{condition}.RDS")
  dat <- readRDS(path)
  #if (rename && "name" %in% names(dat)) dat <- dat %>% rename(geneID = name)
  dat
}

res_data <- purrr::cross_df(list(method = methods, condition = conditions)) %>%
  mutate(
    data = map2(method, condition, ~read_rna(.x, .y)),
    name = paste(.data$method, .data$condition, sep = "_")
  ) %>%
  dplyr::select(name, data) %>%
  deframe()


# nebula 

nebula_names <- c("nebula_age_only", "nebula_age_type1", "nebula_age_type2", "nebula_interaction")

for (nm in nebula_names) {
  res_data[[nm]] <- res_data[[nm]] %>%
    mutate(lfc = lfc / log(2)) %>%
    dplyr::select(name, pval, adj_pval, lfc)
}


# Setup

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue',
  "devil" = "#099668"
)

de_colors <- c("Down-reg" = "steelblue", 
               "Up-reg" = "indianred", 
               "n.s." = "grey")


###----------Volcano plot - explorative-----------###

lfc_cut <- 1.0
pval_cut <- .05

outliers_to_remove <- c("")

prep_rna_data <- function(data, method, condition, pval_col = "adj_pval") {
  
  pvals <- data[[pval_col]]
  min_nonzero <- suppressWarnings(min(pvals[pvals > 0], na.rm = TRUE))
  if (!is.finite(min_nonzero)) min_nonzero <- 1e-300
  
  data %>% 
    dplyr::mutate(
      method = method,
      condition  = condition,
      "{pval_col}" := if_else(.data[[pval_col]] == 0, min_nonzero, .data[[pval_col]])
    ) |>
    dplyr::filter(!name %in% outliers_to_remove)
}

rna_join <- purrr::cross_df(list(method = methods, condition = conditions)) %>% 
  purrr::pmap_dfr(function(method, condition) {
    df <- res_data[[glue::glue("{method}_{condition}")]]
    prep_rna_data(df, method, condition)
  })


rna_join <- rna_join %>% 
  mutate(
    isDE = abs(lfc) >= lfc_cut & adj_pval <= pval_cut,
    DEtype = case_when(
      !isDE ~ "n.s.",
      lfc > 0 ~ "Up-reg",
      TRUE ~ "Down-reg"
    )
  )

p_all <- rna_join %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(color = DEtype), size = 1, alpha = 0.4) +
  scale_color_manual(values = de_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)), 
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  ggh4x::facet_nested(factor(method, levels = c("devil", "glmGamPoi", "nebula"))
                      ~factor(condition, levels = c("age_only", "age_type1", "age_type2", "interaction")), 
                      scales ="free", independent = "y")+
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "DEtype") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'right',
                 legend.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12, color = "black"),  
                 strip.text = element_text(size = 12, face = "plain", color = "black"),
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, color = "black"))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
p_all


ggsave("plot/revision/volcanos_full_joint.png", dpi = 400, width = 12.0, height = 15.0, plot = p_all)


### -------------Classify genes----------------- ###

classify_genes <- function(method_name, res_list) {
  
  # Subset elements for specific method only
  sub <- res_list[stringr::str_detect(names(res_list), paste0("^", method_name, "_"))]
  
  # Read in the four expected test types
  res_age   <- sub[[paste0(method_name, "_age_only")]]   %>% dplyr::select(name, adj_pval, lfc)
  res_type1 <- sub[[paste0(method_name, "_age_type1")]]  %>% dplyr::select(name, adj_pval, lfc)
  res_type2 <- sub[[paste0(method_name, "_age_type2")]]  %>% dplyr::select(name, adj_pval, lfc)
  res_inter <- sub[[paste0(method_name, "_interaction")]]%>% dplyr::select(name, adj_pval, lfc)
  
  # Rename columns
  colnames(res_age)   <- c("name", "padj_Age",   "lfc_Age")
  colnames(res_type1) <- c("name", "padj_Type1", "lfc_Type1")
  colnames(res_type2) <- c("name", "padj_Type2", "lfc_Type2")
  colnames(res_inter) <- c("name", "padj_Inter", "lfc_Inter")
  
  
  df <- purrr::reduce(list(res_age, res_type1, res_type2, res_inter), full_join, by = "name")
  
  # Classification logic
  df %>%
    dplyr::mutate(
      sig_Age   = padj_Age   < pval_cut & abs(lfc_Age)   >= lfc_cut,
      sig_Type1 = padj_Type1 < pval_cut & abs(lfc_Type1) >= lfc_cut,
      sig_Type2 = padj_Type2 < pval_cut & abs(lfc_Type2) >= lfc_cut,
      sig_Inter = padj_Inter < pval_cut & abs(lfc_Inter) >= lfc_cut,
      same_dir  = sign(lfc_Type1) == sign(lfc_Type2),
      
      category = dplyr::case_when(
        sig_Type1 & sig_Type2 & same_dir & !sig_Inter ~ "Shared aging",
        sig_Type1 & !sig_Type2 ~ "Type I specific",
        !sig_Type1 & sig_Type2 ~ "Type II specific",
        sig_Type1 & sig_Type2 & !same_dir ~ "Divergent regulation",
        !sig_Type1 & !sig_Type2 & sig_Inter ~ "Interaction only",
        !sig_Type1 & !sig_Type2 & sig_Age ~ "Age only",
        TRUE ~ "Not significant"
      ),
      method = method_name
    )
}

# Apply to all methods

methods <- unique(stringr::str_extract(names(res_data), "^[^_]+"))
classified_list <- lapply(methods, function(m) classify_genes(m, res_data))
names(classified_list) <- methods

classified_all  <- dplyr::bind_rows(classified_list)
de_genes <- classified_all %>%
  dplyr::filter(category != "Not significant")

de_summary <- de_genes %>%
  dplyr::count(method, category) %>%
  group_by(method) %>%
  dplyr::mutate(
    total = sum(n),
    prop = n / total * 100
  ) %>%
  ungroup()

saveRDS(classified_all, file = "plot/revision/data_to_plot/classified_genes.RDS")


# Stacked barplot

category_colors <- c(
  "Shared aging" = "#AD002AB2",
  "Type I specific" = "#E18727B2",
  "Type II specific" = "#20854Eb2",
  "Interaction only" = "#00468BB2",
  "Age only" = "#FED43999",
  "Divergent regulation" = "#FD88CC99",
  "Not significant" = "darkgray"
)

de_summary$category <- factor(de_summary$category,
                              levels = c("Shared aging", "Type I specific", "Type II specific",
                                         "Interaction only", "Age only", "Divergent regulation")
)

bar <- ggplot(de_summary, aes(x = method, y = prop, fill = category)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = category_colors) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = NULL,
    y = "Proportion of DE genes (%)",
    fill = "Category"
  )

bar

ggsave("plot/revision/barplot_gene_categ.png", dpi = 400, width = 6.0, height = 4.0, plot = bar)


# Upset plot

de_genes <- de_genes %>%
  dplyr::filter(category != "Divergent regulation")

upset_data <- de_genes %>%
  dplyr::select(name, method, category) %>%
  dplyr::mutate(flag = 1) %>%
  pivot_wider(
    names_from = method,
    values_from = flag,
    values_fill = 0
  )

cu_themes <- list(
  intersections_matrix = theme(
    axis.title.x = element_text(),
    axis.title.y = element_text()
  ),
  intersections   = theme(),
  overall_sizes   = theme(),
  sets            = theme()
)

categories <- unique(upset_data$category)
print(categories)

upset_data$category <- factor(
  upset_data$category,
  levels = c("Shared aging", "Type I specific", "Type II specific",
             "Interaction only", "Age only")
)

plots <- lapply(categories, function(cat) {
  cat_data <- upset_data %>% filter(category == cat)
  message("Plotting category: ", cat)
  
  ComplexUpset::upset(
    cat_data,
    intersect = c("devil", "glmGamPoi", "nebula"),
    base_annotations = "auto",
    themes = cu_themes,
    name = cat,
    min_size = 1,
    width_ratio = 0.2
  )
})

joint_upset <- wrap_plots(plots, ncol = 2)
joint_upset

ggsave("plot/revision/upset_gene_categ.png", dpi = 400, width = 12.0, height = 12.0, plot = joint_upset)


# Volcano plot per test, color per category

volcano_long <- classified_all %>%
  dplyr::select(name, method, category,
         lfc_Age, padj_Age,
         lfc_Type1, padj_Type1,
         lfc_Type2, padj_Type2,
         lfc_Inter, padj_Inter) %>%
  pivot_longer(
    cols = starts_with("lfc_"),
    names_to = "test_type",
    names_prefix = "lfc_",
    values_to = "lfc"
  ) %>%
  dplyr::mutate(
    padj = case_when(
      test_type == "Age" ~ padj_Age,
      test_type == "Type1" ~ padj_Type1,
      test_type == "Type2" ~ padj_Type2,
      test_type == "Inter" ~ padj_Inter
    )
  ) %>%
  dplyr::filter(!is.na(lfc), !is.na(padj))


volcano_long$test_type <- factor(volcano_long$test_type,
                                 levels = c("Age", "Type1", "Type2", "Inter"),
                                 labels = c("Age-only", "Type I", "Type II", "Interaction")
)

p_volcanos <- volcano_long %>%
  ggplot(mapping = aes(x = lfc, y = -log10(padj))) +
  geom_point(data = subset(volcano_long, category %in% c("Shared aging")),
             aes(col = category), size = 1.0, alpha = 0.3) +
  geom_point(data = subset(volcano_long, category %in% c("TypeI specific", "TypeII specific", "Interaction only",
                                                              "Age global only", "Divergent regulation",
                                                              "Not significant")),
             aes(col = category), size = 1.0, alpha = 0.5) +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(volcano_long$lfc)), 
                                  ceiling(max(volcano_long$lfc)), by = 2)) +
  ggh4x::facet_nested(factor(method, levels = c("devil", "glmGamPoi", "nebula"))
                      ~factor(test_type, levels = c("Age-only", "Type I", "Type II", "Interaction")), 
                      scales ="free", independent = "y")+
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene category") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'right',
                 legend.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12, color = "black"),  
                 strip.text = element_text(size = 12, face = "plain", color = "black"),
                 axis.text = element_text(size = 12, color = "black"),
                 axis.title = element_text(size = 12, color = "black"))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

p_volcanos

ggsave("plot/revision/volcanos_full_geneCat.png", dpi = 400, width = 14.0, height = 10.0, plot = p_volcanos)


###------------Gene Set Enrichment analysis-------------### 

# Use all genes per test ranked by score

methods <- names(classified_list)
tests <- c("Type1", "Type2", "Inter")

gsea_all <- list()

for (m in methods) {
  df <- classified_list[[m]]
  
  for (t in tests) {
    padj_col <- paste0("padj_", t)
    lfc_col  <- paste0("lfc_",  t)
    combo_name <- paste(m, t, sep = "_")
    
    message("Running GSEA for: ", combo_name)
    gsea_all[[combo_name]] <- enrichmentGO(df, padj_col, lfc_col)
  }
}

devil_Age <- readRDS("results/gsea_GO/old/gseGO_devil.RDS")
glmGamPoi_Age <- readRDS("results/gsea_GO/old/gseGO_glm.RDS")
nebula_Age <- readRDS("results/gsea_GO/old/gseGO_nebula.RDS")

gsea_all <- c(
  gsea_all,
  list(
    devil_Age = devil_Age,
    glmGamPoi_Age = gseGO_glm,
    nebula_Age = nebula_Age
  )
)

saveRDS(gsea_all, file = "results/gsea_GO/new/gsea_all_tests.RDS")

# Simplified

s_cutoff = 0.6

gsea_all_simplified <- lapply(gsea_all, function(x) {
  clusterProfiler::simplify(x, cutoff = s_cutoff)
})

saveRDS(gsea_all_simplified, file = "results/gsea_GO/new/gsea_all_tests_s.RDS")


# Plot Enrichment

gsea_all_simplified <- readRDS("results/gsea_GO/new/gsea_all_tests_s.RDS")
names(gsea_all_simplified) <- c("devil_Type I", "devil_Type II", "devil_Interaction",
                                "glmGamPoi_Type I", "glmGamPoi_Type II", "glmGamPoi_Interaction",
                                "nebula_Type I", "nebula_Type II", "nebula_Interaction",
                                "devil_Age", "glmGamPoi_Age", "nebula_Age")

results_list <- lapply(gsea_all_simplified, function(x) x@result)

results_df <- dplyr::bind_rows(results_list, .id = "source")
results_df <- results_df %>%
  mutate(Regulation = ifelse(NES > 0, "Up", "Down"))

selected_pathways <- c("myofibril assembly",
                       "actin-mediated cell contraction",
                       "muscle cell development",
                       "homophilic cell adhesion via plasma membrane adhesion molecules",
                       "positive regulation of cytokine production",
                       "ncRNA processing",
                       "cellular nitrogen compound catabolic process",
                       "immune system process",
                       "myelin maintenance",
                       "myelination in peripheral nervous system",
                       "peripheral nervous system axon ensheathment",
                       "striated muscle cell development",
                       "regulation of ubiquitin protein ligase activity",
                       "blood vessel endothelial cell migration",
                       "regulation of leukocyte migration",
                       "regulation of cell division",
                       "immune response-activating cell surface receptor signaling pathway",
                       "leukocyte apoptotic process",
                       "regulation of G2/M transition of mitotic cell cycle",
                       "tRNA metabolic process",
                       "regulation of apoptotic signaling pathway",
                       "ribosome biogenesis",
                       "response to endoplasmic reticulum stress",
                       "cellular component assembly involved in morphogenesis",
                       "regulation of substrate adhesion-dependent cell spreading",
                       "T cell mediated immunity",
                       "skeletal muscle cell differentiation",
                       "cellular response to hypoxia",
                       "epithelial cell proliferation",
                       "regulation of voltage-gated calcium channel activity",
                       "regulation of myeloid leukocyte mediated immunity",
                       "nucleosome organization",
                       "aerobic electron transport chain",
                       "positive regulation of muscle contraction",
                       "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress",
                       "cellular oxidant detoxification",
                       "transmission of nerve impulse",
                       "response to interleukin-4",
                       "mast cell activation involved in immune response",
                       "lysosome organization",
                       "leukocyte migration",
                       "leukocyte cell-cell adhesion",
                       "negative regulation of cyclin-dependent protein serine/threonine kinase activity",
                       "mitochondrial ATP synthesis coupled electron transport",
                       "mitochondrial respiratory chain complex assembly",
                       "negative regulation of interleukin-1 production",
                       "regulation of endothelial cell apoptotic process",
                       "positive regulation of cytokine-mediated signaling pathway",
                       "signal transduction in response to DNA damage",
                       "leukocyte apoptotic process")

results_selected <- results_df %>%
  dplyr::filter(Description %in% selected_pathways)


path_cat_gsea <- tibble::tribble(
  ~Description, ~Category,
  # Muscle / Structural Development
  "myofibril assembly", "Muscle / Structural Development",
  "actin-mediated cell contraction", "Muscle / Structural Development",
  "muscle cell development", "Muscle / Structural Development",
  "striated muscle cell development", "Muscle / Structural Development",
  "skeletal muscle cell differentiation", "Muscle / Structural Development",
  "positive regulation of muscle contraction", "Muscle / Structural Development",
  
  # Cell Adhesion & Morphogenesis
  "homophilic cell adhesion via plasma membrane adhesion molecules", "Cell Adhesion / Morphogenesis",
  "cellular component assembly involved in morphogenesis", "Cell Adhesion / Morphogenesis",
  "regulation of substrate adhesion-dependent cell spreading", "Cell Adhesion / Morphogenesis",
  
  # Immune / Inflammatory Response
  "immune system process", "Immune / Inflammatory Response",
  "positive regulation of cytokine production", "Immune / Inflammatory Response",
  "regulation of leukocyte migration", "Immune / Inflammatory Response",
  "immune response-activating cell surface receptor signaling pathway", "Immune / Inflammatory Response",
  "T cell mediated immunity", "Immune / Inflammatory Response",
  "mast cell activation involved in immune response", "Immune / Inflammatory Response",
  "response to interleukin-4", "Immune / Inflammatory Response",
  "regulation of myeloid leukocyte mediated immunity", "Immune / Inflammatory Response",
  "leukocyte migration", "Immune / Inflammatory Response",
  "leukocyte cell-cell adhesion", "Immune / Inflammatory Response",
  "positive regulation of cytokine-mediated signaling pathway", "Immune / Inflammatory Response",
  "negative regulation of interleukin-1 production", "Immune / Inflammatory Response",
  
  # Apoptosis 
  "leukocyte apoptotic process", "Apoptosis",
  "regulation of apoptotic signaling pathway", "Apoptosis",
  "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress", "Apoptosis",
  "regulation of endothelial cell apoptotic process", "Apoptosis",
  
  # Cell Cycle / Division
  "regulation of cell division", "Cell Cycle / Division",
  "regulation of G2/M transition of mitotic cell cycle", "Cell Cycle / Division",
  "negative regulation of cyclin-dependent protein serine/threonine kinase activity", "Cell Cycle / Division",
  "signal transduction in response to DNA damage", "Cell Cycle / Division",
  
  # RNA / Protein Metabolism
  "ncRNA processing", "RNA / Protein Metabolism",
  "tRNA metabolic process", "RNA / Protein Metabolism",
  "ribosome biogenesis", "RNA / Protein Metabolism",
  "regulation of ubiquitin protein ligase activity", "RNA / Protein Metabolism",
  "response to endoplasmic reticulum stress", "RNA / Protein Metabolism",
  
  # Mitochondrial / Oxidative Metabolism
  "aerobic electron transport chain", "Mitochondrial / Oxidative Metabolism",
  "mitochondrial ATP synthesis coupled electron transport", "Mitochondrial / Oxidative Metabolism",
  "mitochondrial respiratory chain complex assembly", "Mitochondrial / Oxidative Metabolism",
  "cellular oxidant detoxification", "Mitochondrial / Oxidative Metabolism",
  
  # Nervous System / Myelination
  "myelin maintenance", "Nervous System / Myelination",
  "myelination in peripheral nervous system", "Nervous System / Myelination",
  "peripheral nervous system axon ensheathment", "Nervous System / Myelination",
  "transmission of nerve impulse", "Nervous System / Myelination",
  
  # Others
  "blood vessel endothelial cell migration", "Vascular Biology",
  "cellular nitrogen compound catabolic process", "Metabolism",
  "nucleosome organization", "Chromatin / Epigenetic Regulation",
  "cellular response to hypoxia", "Stress Response",
  "lysosome organization", "Organelle Organization",
  "epithelial cell proliferation"
)


results_selected <- results_selected %>%
  left_join(path_cat_gsea, by = "Description") 

results_selected$Category[is.na(results_selected$Category)] <- "Other"

heatmap_data <- results_selected %>%
  dplyr::filter(Description %in% selected_pathways) %>%
  dplyr::select(Category, source, Description, enrichmentScore, Regulation)

heatmap_gsea <- ggplot(
  heatmap_data,
  aes(x = source, y = Description, fill = enrichmentScore)
) +
  geom_tile(color = "grey40", linewidth = 0.3) +  
  facet_grid(
    Category ~ Regulation,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name = "Enrichment\nScore",
    limits = c(-1, 1)
  ) +
  labs(
    x = "Method / Test Type",
    y = "GO Biological Process",
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # Axes
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    
    # Facets
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12, color = "black"),
    strip.text.x = element_text(face = "bold", size = 12, color = "black"),
    
    # Panel borders and spacing
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),  # << facet border
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.4, "lines"),
    plot.margin = margin(10, 10, 10, 10),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("plot/revision/heatmap_gsea.png", dpi = 400, width = 15.0, height = 10.0, plot = heatmap_gsea)


###--------------Tissue Enrichment---------------###

get_tissue_specific_res <- function(gse_result, method) {
  # Extract significant GO terms
  significant_terms <- subset(gse_result@result, qvalue < 0.05)
  if (nrow(significant_terms) == 0) return(NULL)
  
  # Collect gene sets
  all_gene_sets <- lapply(significant_terms$ID, function(go_id) {
    genes <- DOSE::geneInCategory(gse_result)[[go_id]]
    genes
  })
  names(all_gene_sets) <- significant_terms$Description
  
  # Run TissueEnrich on each gene set safely
  results <- lapply(seq_along(all_gene_sets), function(i) {
    gene_set <- all_gene_sets[[i]]
    if (length(gene_set) < 5) return(NULL)
    
    gs <- GeneSet(geneIds = gene_set,
                  organism = "Homo Sapiens",
                  geneIdType = SymbolIdentifier())
    tryCatch(
      suppressMessages(TissueEnrich::teEnrichment(gs)),
      error = function(e) NULL
    )
  })
  names(results) <- names(all_gene_sets)
  
  # Extract and combine results
  res_df_list <- lapply(seq_along(results), function(i) {
    te_result <- results[[i]]
    if (is.null(te_result) || length(te_result) == 0) return(NULL)
    
    seOut <- te_result[[1]]
    enrichmentOutput <- setNames(
      data.frame(assay(seOut), row.names = rowData(seOut)[, 1]),
      colData(seOut)[, 1]
    )
    enrichmentOutput$Tissue <- row.names(enrichmentOutput)
    enrichmentOutput$path_name <- names(results)[i]
    enrichmentOutput
  })
  
  res_df_list <- Filter(Negate(is.null), res_df_list)
  if (length(res_df_list) == 0) return(NULL)
  
  res_df <- dplyr::bind_rows(res_df_list)
  res_df$method <- method
  return(res_df)
}

tissue_specific_results <- list()

for (method in names(gsea_all_simplified)) {
  gse_result <- gsea_all_simplified[[method]]
  tissue_specific_results[[method]] <- get_tissue_specific_res(gse_result, method)
}

tissue = target_tissue <- "Skeletal Muscle"

pval_cut = .05
best_df <- lapply(tissue_specific_results, function(tissue_gse) {
  if (is.null(tissue_gse)) return(NULL)
  
  print(unique(tissue_gse$method))
  
  # Main loop across each pathway
  best_per_path <- lapply(unique(tissue_gse$path_name), function(path) {
    r <- tissue_gse %>%
      dplyr::filter(path_name == path) %>%
      dplyr::filter(Log10PValue >= -log10(pval_cut))
    
    if (nrow(r) > 0) {
      r <- na.omit(r)
      r <- r %>% dplyr::filter(fold.change == max(fold.change))
      if (nrow(r) > 1) {
        r <- r %>% dplyr::filter(Tissue.Specific.Genes == max(Tissue.Specific.Genes))
        r <- r[1,]
      }
    } else {
      r <- tissue_gse %>%
        dplyr::filter(path_name == path) %>%
        dplyr::sample_n(1) %>%
        dplyr::mutate(Tissue = "Generic")
    }
    r
  }) %>% dplyr::bind_rows()
  
  best_per_path
}) %>% dplyr::bind_rows()

saveRDS(best_df, file = "plot/revision/data_to_plot/tissue_enrichment.RDS")

# Plot

best_df <- readRDS("plot/revision/data_to_plot/tissue_enrichment.RDS")

plot_tissue_distribution <- function(data, method_colors = NULL) {
  
  df_summary <- data %>% 
    dplyr::group_by(method, Tissue) %>% 
    dplyr::summarise(n = n(), .groups = "drop_last") %>% 
    dplyr::mutate(f = n / sum(n)) %>% 
    dplyr::mutate(
      Tissue = factor(
        Tissue,
        levels = c(
          "Skeletal Muscle", "Generic", "Cerebral Cortex",
          "Adipose Tissue", "Cervix, uterine", "Liver", "Heart Muscle"
        )
      )
    )
  
  ggplot(df_summary, aes(x = Tissue, y = f, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "Tissue",
      y = "Fraction",
      title = ""
    ) +
    scale_fill_manual(values = method_colors)
}

method_colors <- c(
  "devil_Age" = "#099668",
  "devil_Type I" = "#2AA198",
  "devil_Type II" = "#859900",
  "glmGamPoi_Age" = "#EAB578",
  "glmGamPoi_Type I" = "#FFC300",
  "glmGamPoi_Type II" = "#FDE72F",
  "glmGamPoi_Interaction" = "#FF9473",
  "nebula_Age"= "steelblue",
  "nebula_Type I" = "#56B4E9",
  "nebula_Type II" = "#88CCEE"
)
tissue_specific_dist_plot <- plot_tissue_distribution(best_df, method_colors)
tissue_specific_dist_plot

ggsave("plot/revision/tissue_dist_plot.png", dpi = 400, width = 9.0, height = 8.0, plot = tissue_specific_dist_plot)


