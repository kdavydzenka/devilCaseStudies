### Revised downstream analysis using different design combinations ###

rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "glue", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram')
sapply(pkgs, require, character.only = TRUE)


###----------------- Load results-----------------------###

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
    select(name, pval, adj_pval, lfc)
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


###------------Volcano plot--------------###

lfc_cut <- 1.0
pval_cut <- .05

outliers_to_remove <- c("")

prep_rna_data <- function(data, method, lfc_col = "lfc", pval_col = "adj_pval") {
  # get column symbols
  pval_sym <- sym(pval_col)
  
  # find minimal non-zero adjusted p-value
  min_nonzero <- data %>%
    filter(!!pval_sym > 0) %>%
    summarise(min_val = min(!!pval_sym, na.rm = TRUE)) %>%
    pull(min_val)
  
  # clean and annotate
  data %>%
    mutate(
      method = method,
      !!pval_col := if_else(
        !!pval_sym == 0,
        min_nonzero,
        !!pval_sym
      )
    ) %>%
    dplyr::filter(!name %in% outliers_to_remove)
}


methods <- c("devil", "glmGamPoi", "nebula")

rna_join <- map_df(methods, function(m) {
  prep_rna_data(res_data[[glue("{m}_interaction")]], m)
})

rna_join <- rna_join %>%
  mutate(
    isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
    DEtype = case_when(
      !isDE ~ "n.s.",
      lfc > 0 ~ "Up-reg",
      TRUE ~ "Down-reg"
    )
  )


p4 <- ggplot(rna_join, aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(color = DEtype), size = 1.5, alpha = 0.3) +
  scale_color_manual(values = de_colors) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 12) +
  facet_wrap(~factor(method, levels = methods), nrow = 1, scales = "free") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    name = expression(Log[2]~FC)
  ) +
  labs(
    y = expression(-log[10]~adjusted~P[value]),
    color = "DE type",
    title = "Interaction"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

p4


join_p <- (p1 / p2 / p3 / p4)

ggsave("plot/revision/volcanos_full_joint.png", dpi = 400, width = 12.0, height = 15.0, plot = join_p)


### ---------Classify genes-------------- ###

classify_genes <- function(method_name, res_list) {
  
  # Subset elements for specific method only
  sub <- res_list[stringr::str_detect(names(res_list), paste0("^", method_name, "_"))]
  
  # Read in the four expected test types
  res_age   <- sub[[paste0(method_name, "_age_only")]]   %>% select(name, adj_pval, lfc)
  res_type1 <- sub[[paste0(method_name, "_age_type1")]]  %>% select(name, adj_pval, lfc)
  res_type2 <- sub[[paste0(method_name, "_age_type2")]]  %>% select(name, adj_pval, lfc)
  res_inter <- sub[[paste0(method_name, "_interaction")]]%>% select(name, adj_pval, lfc)
  
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
        sig_Type1 & !sig_Type2 ~ "TypeI specific",
        !sig_Type1 & sig_Type2 ~ "TypeII specific",
        sig_Type1 & sig_Type2 & !same_dir ~ "Divergent regulation",
        !sig_Type1 & !sig_Type2 & sig_Inter ~ "Interaction only",
        !sig_Type1 & !sig_Type2 & sig_Age ~ "Age global only",
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


#----------- Stacked barplot-------------#

category_colors <- c(
  "Shared aging" = "#fba59b",
  "TypeI specific" = "#a57c90",
  "TypeII specific" = "#94afb8",
  "Interaction only" = "#F0E442",
  "Age global only" = "#B6C687",
  "Divergent regulation" = "#FCC88A",
  "Not significant" = "lightgray"
)

de_summary$category <- factor(de_summary$category,
                              levels = c("Shared aging", "TypeI specific", "TypeII specific",
                                         "Interaction only", "Age global only", "Divergent regulation")
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


# ---Volcano plot per test, color per category------#

ggplot(devil, aes(x = lfc_Inter, y = -log10(padj_Inter), color = category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Interaction volcano plot (Age × Cell type)",
    x = "log₂ Fold Change (Interaction effect)",
    y = "−log₁₀ Adjusted p-value"
  )


volcano_long <- classified_all %>%
  select(name, method, category,
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
  mutate(
    padj = case_when(
      test_type == "Age" ~ padj_Age,
      test_type == "Type1" ~ padj_Type1,
      test_type == "Type2" ~ padj_Type2,
      test_type == "Inter" ~ padj_Inter
    )
  ) %>%
  filter(!is.na(lfc), !is.na(padj))


volcano_long$test_type <- factor(volcano_long$test_type,
                                 levels = c("Age", "Type1", "Type2", "Inter"),
                                 labels = c("Age-only", "Type I", "Type II", "Interaction")
)

volcano <- ggplot(volcano_long, aes(x = lfc, y = -log10(padj), color = category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  facet_grid(method ~ test_type, scales = "free") +
  #scale_x_continuous(
    #breaks = scales::pretty_breaks(n = 6),
    #name = expression(Log[2]~FC)
  #) +
  labs(
    y = expression(-log[10]~adjusted~P[value]),
    x = "log₂ Fold Change",
    color = "Gene category",
    title = ""
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
volcano

ggsave("plot/revision/volcanos_full_geneCat.png", dpi = 400, width = 14.0, height = 12.0, plot = volcano)



### Compare in full lfc_age vs lfc_age_cellType, pval_age vs pval_age_cellType ###

#compare_age_cell <- function(method, data_list,
                                    #lfc_col = "lfc",
                                    #pval_col = "adj_pval") {
  #age <- data_list[[glue("rna_{method}_age")]]
  #age_cell <- data_list[[glue("rna_{method}_interaction")]]
  
  #merged <- inner_join(
    #age %>% select(geneID, lfc_age = all_of(lfc_col), pval_age = all_of(pval_col)),
    #age_cell %>% select(geneID, lfc_cell = all_of(lfc_col), pval_cell = all_of(pval_col)),
    #by = "geneID"
  #)
  
  # Plot 1: log2FC comparison
  #p1 <- ggplot(merged, aes(x = lfc_age, y = lfc_cell)) +
    #geom_point(alpha = 0.3, color = "black") +
    #geom_smooth(method = "lm", se = T, color = "darkred") +
    #labs(
      #title = glue("LFC: {method}"),
      #x = "log2FC (Age)",
      #y = "log2FC (Interaction)"
    #) +
    #theme_minimal()
  
  # Plot 2: adjusted p-value comparison 
  #p2 <- ggplot(merged, aes(x = -log10(pval_age), y = -log10(pval_cell))) +
    #geom_point(alpha = 0.3, color = "black") +
    #geom_smooth(method = "lm", se = T, color = "darkred") +
    #labs(
      #title = glue("adj_pval: {method}"),
      #x = "-log10(adj_pval, Age)",
      #y = "-log10(adj_pval, Interaction)"
    #) +
    #theme_minimal()
  
  #list(lfc_plot = p1, pval_plot = p2)
#}


# Run comparisons for all methods

#plots <- map(c("devil", "glmGamPoi", "nebula"),
             #~compare_age_cell(.x, rna_data))
#names(plots) <- c("devil", "glmGamPoi", "nebula")

# Access plots
#plots[["devil"]]$lfc_plot  
#plots[["devil"]]$pval_plot 

#plots[["glmGamPoi"]]$lfc_plot  
#plots[["glmGamPoi"]]$pval_plot  

#plots[["nebula"]]$lfc_plot  
#plots[["nebula"]]$pval_plot  

#final_plot <- (
  #plots[["devil"]]$lfc_plot + plots[["devil"]]$pval_plot
#) /
  #(
    #plots[["glmGamPoi"]]$lfc_plot + plots[["glmGamPoi"]]$pval_plot
  #) /
  #(
    #plots[["nebula"]]$lfc_plot + plots[["nebula"]]$pval_plot
  #)

#final_plot

#ggsave("plot/revision/scatter_full_age_vs_interaction.png", dpi = 400, width = 10.0, height = 10.0, plot = final_plot)



