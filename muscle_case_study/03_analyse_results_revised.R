### Revised downstream analysis using different design combinations ###

rm(list = ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/muscle_case_study/")

pkgs <- c("ggplot2", "dplyr","tidyr", "purrr", "tibble", "viridis", "smplot2", "Seurat", "gridExtra",
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
rna_nebula_age <- rna_data$rna_nebula_age %>%
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
rna_nebula_int <- rna_data$rna_nebula_interaction %>%
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


# Gene selection based on LFC & pvalue cutoff 

lfc_cut <- 1.0
pval_cut <- .05

rna_devil <- rna_data$rna_devil_age %>%
  #dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "devil")

rna_glm <- rna_data$rna_glmGamPoi_age %>%
  #dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "glmGamPoi")

rna_nebula <- rna_nebula_age %>%
  #dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>%
  dplyr::mutate(method = "nebula")


### Volcano plot ###

rna_devil$adj_pval[rna_devil$adj_pval == 0] <- min(rna_devil$adj_pval[rna_devil$adj_pval != 0])
rna_devil <- rna_devil[rna_devil$geneID != "KCTD1", ] # remove outlier
rna_nebula$adj_pval[rna_nebula$adj_pval == 0] <- min(rna_nebula$adj_pval[rna_nebula$adj_pval != 0])
rna_glm$adj_pval[rna_glm$adj_pval == 0] <- min(rna_glm$adj_pval[rna_glm$adj_pval != 0])

rna_join <- rbind(rna_devil, rna_glm, rna_nebula)
rna_join <- rna_join %>%
  dplyr::mutate(
    isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
    DEtype = if_else(!isDE, "n.s.", if_else(lfc > 0, "Up-reg", "Down-reg")))


p_volcanos <- rna_join %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(col = DEtype), size = 2.0, alpha = 0.2) +
  scale_color_manual(values = de_colors) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)), 
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
  labs(x = expression(Log[2] ~ FC),
       y = expression(-log[10] ~ Pvalue),
       col = "DE type") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p_volcanos

ggsave("plot/revision/subsampled/volcanos_interaction_subsampled.png", dpi = 400, width = 16.0, height = 5.0, plot = p_volcanos)


### Compare in full lfc_age vs lfc_inter, pval_age vs pval_inter ###





