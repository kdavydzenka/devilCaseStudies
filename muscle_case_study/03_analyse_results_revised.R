### Revised downstream analysis using different design combinations ###

rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "glue", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram')
sapply(pkgs, require, character.only = TRUE)


### MuscleRNA results ###

methods <- c("devil", "glmGamPoi", "nebula")
conditions <- c("age", "interaction")

read_rna <- function(method, condition, rename = TRUE) {
  path <- glue::glue("results/MuscleRNA/full/{method}_{condition}_rna.RDS")
  dat <- readRDS(path)
  if (rename && "name" %in% names(dat)) dat <- dat %>% rename(geneID = name)
  dat
}

rna_data <- purrr::cross_df(list(method = methods, condition = conditions)) %>%
  mutate(
    data = map2(method, condition, ~read_rna(.x, .y, rename = .x != "nebula")),
    name = paste("rna", .data$method, .data$condition, sep = "_")
  ) %>%
  dplyr::select(name, data) %>%
  deframe()


# nebula - age (full)
rna_data$rna_nebula_age <- rna_data$rna_nebula_age %>%
  dplyr::mutate(
    geneID = name,  
    lfc = lfc / log(2)  
  ) %>% 
  dplyr::select(geneID,pval,adj_pval,lfc)

# nebula - age+cellType (subsampled)
rna_nebula_age_cellType <- rna_data$rna_nebula_age_cellType %>%
  dplyr::mutate(
    geneID = gene,  
    pval = `p_age_cluster1`,  
    adj_pval = p.adjust(`p_age_cluster1`, method = "BH"),
    lfc = `logFC_age_cluster1` / log(2)  
  ) %>% 
  dplyr::select(geneID,pval,adj_pval,lfc)

# nebula - interaction 
rna_data$rna_nebula_interaction <- rna_data$rna_nebula_interaction %>%
  dplyr::mutate(
    geneID = gene,  
    pval = `p_age_cluster1:cell_typeType II`,  
    adj_pval = p.adjust(`p_age_cluster1:cell_typeType II`, method = "BH"),
    lfc = `logFC_age_cluster1:cell_typeType II` / log(2)  
  ) %>% 
  dplyr::select(geneID,pval,adj_pval,lfc)
  
  
# Setup

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  'steelblue',
  "devil" = "#099668"
)

de_colors <- c("Down-reg" = "steelblue", 
               "Up-reg" = "indianred", 
               "n.s." = "grey")


### Volcano plot ###

lfc_cut <- 1.0
pval_cut <- .05

outliers_to_remove <- c("KCTD1")

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
    filter(!geneID %in% outliers_to_remove)
}


#rna_devil$adj_pval[rna_devil$adj_pval == 0] <- min(rna_devil$adj_pval[rna_devil$adj_pval != 0])
#rna_devil <- rna_devil[rna_devil$geneID != "KCTD1", ] # remove outlier
#rna_nebula$adj_pval[rna_nebula$adj_pval == 0] <- min(rna_nebula$adj_pval[rna_nebula$adj_pval != 0])
#rna_glm$adj_pval[rna_glm$adj_pval == 0] <- min(rna_glm$adj_pval[rna_glm$adj_pval != 0])

methods <- c("devil", "glmGamPoi", "nebula")

rna_join <- map_df(methods, function(m) {
  prep_rna_data(rna_data[[glue("rna_{m}_interaction")]], m)
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


p_interaction <- ggplot(rna_join, aes(x = lfc, y = -log10(adj_pval))) +
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

p_interaction

join_p <- (p_age / p_interaction)

ggsave("plot/revision/full/volcanos_age_int.png", dpi = 400, width = 12.0, height = 10.0, plot = join_p)


### Compare in full lfc_age vs lfc_inter, pval_age vs pval_inter ###

compare_age_interaction <- function(method, data_list,
                                    lfc_col = "lfc",
                                    pval_col = "adj_pval") {
  age <- data_list[[glue("rna_{method}_age")]]
  inter <- data_list[[glue("rna_{method}_interaction")]]
  
  merged <- inner_join(
    age %>% select(geneID, lfc_age = all_of(lfc_col), pval_age = all_of(pval_col)),
    inter %>% select(geneID, lfc_inter = all_of(lfc_col), pval_inter = all_of(pval_col)),
    by = "geneID"
  )
  
  # Plot 1: log2FC comparison
  p1 <- ggplot(merged, aes(x = lfc_age, y = lfc_inter)) +
    geom_point(alpha = 0.3, color = "black") +
    geom_smooth(method = "lm", se = T, color = "darkred") +
    labs(
      title = glue("LFC: {method}"),
      x = "log2FC (Age)",
      y = "log2FC (Interaction)"
    ) +
    theme_minimal()
  
  # Plot 2: adjusted p-value comparison 
  p2 <- ggplot(merged, aes(x = -log10(pval_age), y = -log10(pval_inter))) +
    geom_point(alpha = 0.3, color = "black") +
    geom_smooth(method = "lm", se = T, color = "darkred") +
    labs(
      title = glue("adj_pval: {method}"),
      x = "-log10(adj_pval, Age)",
      y = "-log10(adj_pval, Interaction)"
    ) +
    theme_minimal()
  
  list(lfc_plot = p1, pval_plot = p2)
}


# Run comparisons for all methods

plots <- map(c("devil", "glmGamPoi", "nebula"),
             ~compare_age_interaction(.x, rna_data))
names(plots) <- c("devil", "glmGamPoi", "nebula")

# Access plots
#plots[["devil"]]$lfc_plot  
#plots[["devil"]]$pval_plot 

#plots[["glmGamPoi"]]$lfc_plot  
#plots[["glmGamPoi"]]$pval_plot  

#plots[["nebula"]]$lfc_plot  
#plots[["nebula"]]$pval_plot  

final_plot <- (
  plots[["devil"]]$lfc_plot + plots[["devil"]]$pval_plot
) /
  (
    plots[["glmGamPoi"]]$lfc_plot + plots[["glmGamPoi"]]$pval_plot
  ) /
  (
    plots[["nebula"]]$lfc_plot + plots[["nebula"]]$pval_plot
  )

final_plot

ggsave("plot/revision/full/scatter_full.png", dpi = 400, width = 10.0, height = 10.0, plot = final_plot)

