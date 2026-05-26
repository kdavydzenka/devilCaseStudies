
rm(list = ls())
require(tidyverse)
require(patchwork)
library(ggplot2)

dir.create("source_data", recursive = T)

# Figure 2 ####
# Use MacaqueBrain and HSC
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/runtime.RDS")
panelA_runtime = pA@data %>% dplyr::select(n_genes, n_cells, model_name, y)
# panelA_runtime %>% 
#   ggplot(mapping = aes(x = factor(n_genes), y = y, fill = model_name)) + 
#   facet_wrap(~n_cells, scales = "free_y") + 
#   scale_y_sqrt() + 
#   geom_col(position = "dodge")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory.RDS")
panelA_memory = pB@data %>% dplyr::select(n_genes, n_cells, model_name, y)
panelA_memory %>% 
  ggplot(mapping = aes(x = factor(n_genes), y = y, fill = model_name)) + 
  facet_wrap(~n_cells, scales = "free_y") + 
  scale_y_sqrt() + 
  geom_col(position = "dodge")

p15mill = readRDS("timing_scaling/img/RDS/MacaqueBrain/large.RDS")
panelB_data = p15mill@data %>% dplyr::select(model_name, n_genes, n_cells, time)
panelB_data %>%
  ggplot(mapping = aes(x = n_cells, y = time, col = model_name)) +
  geom_point() +
  geom_line() +
  scale_y_sqrt() +
  scale_x_log10()

# panel d
pD = readRDS("de_analysis/nullpower/figures/RDS/main/qq_plot.rds")
panelD_data = pD@data %>% dplyr::select(name, is.pb, x, y)
panelD_data %>% 
  ggplot(mapping = aes(x = x, y = y, colour = name)) +
  geom_line(linewidth = 1.2) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = 2,
    linewidth = 0.5,
    colour = "black"
  ) +
  facet_wrap(~is.pb) +
  scale_y_continuous(limits = c(0, 5))

# panel e
pE = readRDS("de_analysis/nullpower/figures/RDS/main/power_curve.rds")
panelE_data = pE@data %>% dplyr::select(FDR, TPR, method, splitval)
panelE_data %>%
  ggplot(mapping = aes(x = FDR, y = TPR, colour = method, group = method)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~splitval)

# panel f
p = readRDS("de_analysis/nullpower/figures/RDS/main/MCC_boxplot.rds")
panelF_data = p@data %>% dplyr::select(name, MCC, patients, author)
panelF_data %>%
  ggplot(mapping = aes(x = name, y = MCC, col = name)) +
  geom_boxplot() +
  ggh4x::facet_nested(~is.pb + patients) +
  coord_flip()

# Save to xlsx file with one sheet per panel (split panel a in memory and runtime)
require(openxlsx)

source_data_fig2 <- list(
  panelA_runtime = panelA_runtime,
  panelA_memory  = panelA_memory,
  panelB         = panelB_data,
  panelD         = panelD_data,
  panelE         = panelE_data,
  panelF         = panelF_data
)

write.xlsx(
  source_data_fig2,
  file = "source_data/Figure2_source_data.xlsx",
  overwrite = TRUE
)

rm(list = ls())
gc()

# Figure 3 ####
pA = readRDS("cell_types_analysis/plot_figure/main_RDS/panelA.RDS")
pA_data <- pA$layers[[1]]$data %>%
  dplyr::select(umap1, umap2, legend)
pA_data %>%
  ggplot(mapping = aes(x=umap1, y=umap2, col = legend)) +
  geom_point()

pB = readRDS("cell_types_analysis/plot_figure/main_RDS/panelB.RDS")
pB_data = pB@data %>% dplyr::select(name, pval, adj_pval, lfc, is_marker)
pB_data = pB_data %>% dplyr::filter(lfc > -5)
pB_data %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval), col = is_marker)) +
  geom_point()

pD = readRDS("cell_types_analysis/plot_figure/main_RDS/panelD.RDS")
panelD_data <- dplyr::bind_rows(
  pD$layers[[1]]$data %>%
    dplyr::select(annot, pred, score) %>%
    dplyr::mutate(layer = "all_scores", is_correct = NA),
  
  pD$layers[[2]]$data %>%
    dplyr::select(annot, pred, score, is_correct) %>%
    dplyr::mutate(layer = "max_prediction"),
  
  pD$layers[[3]]$data %>%
    dplyr::select(annot, pred, score, is_correct) %>%
    dplyr::mutate(layer = "max_prediction")
)
panelD_data %>%
  ggplot(mapping = aes(x = annot, y = pred, size = score, colour = is_correct)) +
  geom_point()

pE = readRDS("cell_types_analysis/plot_figure/main_RDS/panelE.RDS")
panelE_data = pE@data %>% dplyr::select(model, n_markers, mean_acc, sd_acc)
panelE_data %>%
  ggplot(mapping = aes(x = n_markers, y = mean_acc, ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc, col = model)) +
  geom_pointrange() +
  geom_line() +
  scale_x_continuous(transform = "log10")

pF = readRDS("cell_types_analysis/plot_figure/main_RDS/panelF.RDS")
panelF_data = pF@data %>% dplyr::select(cluster, n, model)
panelF_data %>%
  ggplot(mapping = aes(x = cluster, y = n, fill = model)) +
  geom_col(position = "dodge")

pG = readRDS("cell_types_analysis/plot_figure/main_RDS/panelG.RDS")
panelG_data = pG@data %>% dplyr::select(cluster, model, fold_change)
panelG_data %>% ggplot(mapping = aes(x = as.numeric(cluster), y = fold_change, col = model)) + geom_line() + geom_point()

# Save to xlsx file with one sheet per panel
require(openxlsx)

source_data_fig3 <- list(
  panelA = pA_data,
  panelB = pB_data,
  panelD = panelD_data,
  panelE = panelE_data,
  panelF = panelF_data,
  panelG = panelG_data
)

write.xlsx(
  source_data_fig3,
  file = "source_data/Figure3_source_data.xlsx",
  overwrite = TRUE
)

rm(list = ls())
gc()

# Figure 4 #####
pA = readRDS("cell_types_analysis/marker_finder/results/umap_data.rds")
panelA_data = pA
panelA_data %>% ggplot(mapping = aes(x = UMAP_1, y = UMAP_2, colour = cell_type)) +
  geom_point()

pB = readRDS("cell_types_analysis/marker_finder/img/auc_bar_plot.RDS")
panelB_and_C_data = pB@data
panelB_and_C_data %>% ggplot(mapping = aes(x = cell_type, y = auc, fill = model)) + geom_col(position = "dodge") + coord_flip()
panelB_and_C_data %>% ggplot(mapping = aes(x = model, y = auc, fill = model)) + geom_boxplot()

# Save to xlsx file with one sheet per panel
require(openxlsx)

source_data_fig4 <- list(
  panelA       = panelA_data,
  panelB_and_C = panelB_and_C_data
)

write.xlsx(
  source_data_fig4,
  file = "source_data/Figure4_source_data.xlsx",
  overwrite = TRUE
)

rm(list = ls())
gc()

# Figure 5 #####
pA = readRDS("muscle_case_study_revision/figures/umap_all.rds")
panelA_data = pA@data %>% dplyr::select(umap_1, umap_2, cell_type, age_pop)
panelA_data %>% 
  ggplot(mapping = aes(x = umap_1, y = umap_2, colour = age_pop)) +
  geom_point()

pB = readRDS("muscle_case_study_revision/figures/umap_zoom.rds")
panelB_data = pB@data %>% dplyr::select(umap_1, umap_2, cell_type, age_pop)
panelB_data %>%
  ggplot(mapping = aes(x = umap_1, y = umap_2, colour = age_pop)) +
  geom_point()

pC = readRDS("muscle_case_study_revision/figures/volcano_devil.RDS")
panelC_data = pC@data
panelC_data %>% ggplot(mapping = aes(x = lfc, y = -log10(padj), col = category)) +
  geom_point() +
  facet_wrap(~test_type)

pD = readRDS("muscle_case_study_revision/figures/go_plot_devil.rds")
panelD_data = pD@data
panelD_data %>% 
  ggplot(mapping = aes(x = Description, y = -log10(p.adjust), size = abs(enrichmentScore))) +
  geom_point() +
  facet_wrap(~condition, space = "free_y") +
  coord_flip()

pE = readRDS("muscle_case_study_revision/figures/semantic_plot_overlay.rds")
panelE_data = pE@data
panelE_data %>% 
  ggplot(mapping = aes(x = dim1, y = dim2, size = size, col = method)) +
  geom_point() +
  facet_wrap(~condition)

# Save to xlsx file with one sheet per panel
require(openxlsx)

source_data_fig5 <- list(
  panelA = panelA_data,
  panelB = panelB_data,
  panelC = panelC_data,
  panelD = panelD_data,
  panelE = panelE_data
)

write.xlsx(
  source_data_fig5,
  file = "source_data/Figure5_source_data.xlsx",
  overwrite = TRUE
)

rm(list = ls())
gc()

# Figure 6 ####

pB = readRDS("muscle_case_study_revision/results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/gsea_tissue_plot.RDS")
panelB_data = pB@data %>% dplyr::select(method, Tissue, n, f)
panelB_data %>% ggplot(mapping = aes(x = Tissue, y = f, fill = method)) +
  geom_col(position = "dodge")

pC = readRDS("muscle_case_study_revision/results/MuscleRNA/per_contrast_vector_analysis/sub_v_full_comparison/age_type1/venn_bar_plot_gsea.RDS")
panelC_data = pC@data
panelC_data %>% 
  ggplot(mapping = aes(x = method, y = prop, fill = name)) +
  geom_col()
  
pD = readRDS("muscle_case_study_revision/results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/venn_bar_plot.RDS")
panelD_data = pD@data
panelD_data %>% 
  ggplot(mapping = aes(x = "", y = value, fill = name)) +
  geom_col()

umaps = readRDS("muscle_case_study_revision/results/MuscleRNA/per_contrast_vector_analysis/full/age_type1/umaps.RDS")
pE = umaps$`glmGamPoi private`
panelE_data = pE@data
panelE_data %>% ggplot(mapping = aes(x=umap_1, y=umap_2, colour = age_cluster)) + geom_point()

pF = umaps$`glmGamPoi and devil`
panelF_data = pF@data
panelF_data %>% ggplot(mapping = aes(x=umap_1, y=umap_2, colour = age_cluster)) + geom_point()

# Save to xlsx file with one sheet per panel
require(openxlsx)

source_data_fig6 <- list(
  panelB = panelB_data,
  panelC = panelC_data,
  panelD = panelD_data,
  panelE = panelE_data,
  panelF = panelF_data
)

write.xlsx(
  source_data_fig6,
  file = "source_data/Figure6_source_data.xlsx",
  overwrite = TRUE
)

rm(list = ls())
gc()
