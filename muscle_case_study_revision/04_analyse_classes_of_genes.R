### Revised downstream analysis using different design combinations ###

rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "glue", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ggVennDiagram",
          "clusterProfiler", "org.Hs.eg.db", "ComplexUpset", "TissueEnrich")
sapply(pkgs, require, character.only = TRUE)
source("utils/utils_analysis.R")
source("utils/utils.R")


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

# ggsave("plot/revision/volcanos_full_joint.png", dpi = 400, width = 14.0, height = 10.0, plot = p_all)


### -------------Classify genes----------------- ###
classify_genes <- function(method_name, res_list,
                           pval_cut = 0.05,
                           lfc_cut  = 1.00) {
  
  # Subset elements for specific method only
  sub <- res_list[stringr::str_detect(names(res_list), paste0("^", method_name, "_"))]
  
  # Check that all expected tests are present
  expected_suffix <- c("age_only", "age_type1", "age_type2", "interaction")
  expected_names  <- paste0(method_name, "_", expected_suffix)
  missing <- setdiff(expected_names, names(sub))
  if (length(missing) > 0) {
    stop("Missing results for: ", paste(missing, collapse = ", "))
  }
  
  # Read in the four expected test types
  res_type1 <- sub[[paste0(method_name, "_age_type1")]]  %>%
    dplyr::select(name, adj_pval, lfc)
  res_type2 <- sub[[paste0(method_name, "_age_type2")]]  %>%
    dplyr::select(name, adj_pval, lfc)
  res_inter <- sub[[paste0(method_name, "_interaction")]] %>%
    dplyr::select(name, adj_pval, lfc)
  
  # Rename columns
  colnames(res_type1) <- c("name", "padj_Type1", "lfc_Type1")
  colnames(res_type2) <- c("name", "padj_Type2", "lfc_Type2")
  colnames(res_inter) <- c("name", "padj_Inter", "lfc_Inter")
  
  df <- purrr::reduce(list(res_type1, res_type2, res_inter),
                      dplyr::full_join, by = "name")
  
  df = df %>%
    dplyr::mutate(
      sig_Type1 = !is.na(padj_Type1) & padj_Type1 < pval_cut & !is.na(lfc_Type1) & abs(lfc_Type1) >= lfc_cut,
      sig_Type2 = !is.na(padj_Type2) & padj_Type2 < pval_cut & !is.na(lfc_Type2) & abs(lfc_Type2) >= lfc_cut,
      sig_Inter = !is.na(padj_Inter) & padj_Inter < pval_cut & !is.na(lfc_Inter) & abs(lfc_Inter) >= lfc_cut,
      
      same_dir  = dplyr::if_else(
        is.na(lfc_Type1) | is.na(lfc_Type2),
        FALSE,
        sign(lfc_Type1) == sign(lfc_Type2)
      ),
      
      category = dplyr::case_when(
        sig_Type1 & sig_Type2 & same_dir & !sig_Inter ~ "Shared aging",
        sig_Type1 & sig_Type2 & same_dir & sig_Inter ~ "Shared aging + interaction",
        sig_Type1 & !sig_Type2              ~ "Type I specific",
        !sig_Type1 & sig_Type2              ~ "Type II specific",
        sig_Type1 & sig_Type2 & !same_dir   ~ "Divergent regulation",
        !sig_Type1 & !sig_Type2 & sig_Inter ~ "Interaction only",
        TRUE                                ~ "Not significant"
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

#saveRDS(classified_all, file = "plot/revision/data_to_plot/classified_genes.RDS")


# Stacked barplot

category_colors <- c(
  "Shared aging" = "#AD002AB2",
  "Shared aging + interaction" = "#FED43999",
  "Type I specific" = "#E18727B2",
  "Type II specific" = "#20854Eb2",
  "Interaction only" = "#00468BB2",
  "Divergent regulation" = "#FD88CC99",
  "Not significant" = "darkgray"
)

de_summary$category <- factor(de_summary$category)

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

# Volcano plot per test, color per category

volcano_long <- classified_all %>%
  dplyr::select(name, method, category,
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

specific_genes = c("SLC2A4", "MYH2", "MYH1", "PER2", "CTSD", "SOCS3", "TGFB2", "PCDHGA1", "SAA2", "ID1", "TNNT2", "PER2", "AMPD3", "JUN")

volcano_long <- volcano_long %>% dplyr::mutate(is_aging_marker = name %in% specific_genes)
volcano_long_devil = volcano_long %>% dplyr::filter(method == "devil")

p_volcano_devil <- volcano_long_devil %>%
  dplyr::filter(method == "devil") %>% 
  ggplot(mapping = aes(x = lfc, y = -log10(padj))) +
  geom_point(data = volcano_long_devil,
             aes(col = category), size = 1.0, alpha = 0.5) +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(volcano_long_devil$lfc)), 
                                  ceiling(max(volcano_long_devil$lfc)), by = 2)) +
  ggh4x::facet_nested(factor(test_type, levels = c("Age-only", "Type I", "Type II", "Interaction"))~
                        factor(method, levels = c("devil", "glmGamPoi", "nebula")), 
                      scales ="free", independent = "y")+
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene category") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  geom_label_repel(
    data = subset(volcano_long_devil, is_aging_marker & padj < pval_cut & abs(lfc) > lfc_cut),
    aes(label = name),
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
  )
p_volcano_devil 


p_volcano_all <- volcano_long %>%
  ggplot(mapping = aes(x = lfc, y = -log10(padj))) +
  geom_point(data = volcano_long,
             aes(col = category), size = 1.0, alpha = 0.5) +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(volcano_long$lfc)), 
                                  ceiling(max(volcano_long$lfc)), by = 2)) +
  ggh4x::facet_nested(factor(test_type, levels = c("Age-only", "Type I", "Type II", "Interaction"))~
                        factor(method, levels = c("devil", "glmGamPoi", "nebula")), 
                      scales ="free", independent = "y")+
  labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "Gene category") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  geom_label_repel(
    data = subset(volcano_long, is_aging_marker & padj < pval_cut & abs(lfc) > lfc_cut),
    aes(label = name),
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
  )
p_volcano_all

saveRDS(p_volcano_devil, "figures/volcano_devil.RDS")
saveRDS(p_volcano_all, "figures/volcano_all.RDS")
