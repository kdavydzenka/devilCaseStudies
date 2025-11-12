rm(list = ls())
require(patchwork)
require(tidyverse)

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)", "flash")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)", "flash")

method_colors = c(
  "glmGamPoi (fixed)" = "#A22E29",
  "edgeR" = "#7D629E",
  "edgeR (Pb)" = "#7D629E",
  "limma" = "#B96461",
  "limma (Pb)" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "glmGamPoi" = "#EAB578",
  "Nebula" =  'steelblue', #"#B0C4DE",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "Devil (base)" = "#099668",
  "Devil (mixed)" = "#099668",
  "Devil" = "#099668",
  "devil" = "#099668",
  "flash" = "black"
)

res = readRDS("nullpower/final_res/results.rds") %>%
  dplyr::filter((is.pb == TRUE & name %in% method_patientwise) | (is.pb == FALSE & name %in% method_cellwise))

all_timing = lapply(list.files("nullpower/timing_results/"), function(p) {
  readRDS(file.path("nullpower/timing_results/", p)) %>% dplyr::mutate(dataset = p)
}) %>% do.call("bind_rows", .) %>%
  dplyr::filter((is.pb == TRUE & algo %in% method_patientwise) | (is.pb == FALSE & algo %in% method_cellwise))

all_MCC_boxplots = res %>%
  dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
  dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
  ggplot(mapping = aes(x=name, y=MCC, col=name)) +
  geom_boxplot() +
  facet_grid(author~is.pb) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(color = "") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
all_MCC_boxplots

all_time_boxplots = all_timing %>%
  dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
  dplyr::mutate(name = algo) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
  dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
  ggplot(mapping = aes(x=name, y=timings, col=name)) +
  geom_boxplot() +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(color = "") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(transform = "log10") +
  coord_flip() +
  facet_wrap(~author) +
  labs(y = "Time (s)")

failure_rate_plot = all_timing %>%
  dplyr::filter(is.pb) %>%
  dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
  dplyr::mutate(name = algo) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
  dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
  dplyr::mutate(is_bad = is.na(timings)) %>%
  dplyr::group_by(name, is.pb, author) %>%
  dplyr::summarise(`Failure rate` = sum(is_bad) / n()) %>%
  ggplot(mapping = aes(x=name, y=`Failure rate`, fill=name, col=name)) +
  geom_col(col="black") +
  facet_grid(author~is.pb) +
  scale_fill_manual(values = method_colors) +
  labs(color = "") +
  theme_bw() +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )

design = "
A
A
A
B
B
C
C"

final_plot = free(all_MCC_boxplots) + free(all_time_boxplots) + free(failure_rate_plot) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a")

ggsave("figures/all_models_comparison.pdf", width = 11, height = 16, units = "in", dpi = 600, plot = final_plot)


method_patientwise = c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)", "flash")
res$name %>% unique()
res %>% 
  dplyr::filter(is.pb == T, name %in% method_patientwise) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "Devil", name)) %>% 
  ggplot(aes(MCC, colour = name)) +
  stat_ecdf(size = 1) +
  theme_bw() +
  facet_wrap(~author) +
  scale_color_manual(values = method_colors)


method_cellwise = c("Nebula", "Devil (base)", "limma", "glmGamPoi (cell)", "flash")
res$name %>% unique()
res %>% 
  dplyr::filter(is.pb == F, name %in% method_cellwise) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "Devil", name)) %>% 
  ggplot(aes(MCC, colour = name)) +
  stat_ecdf(size = 1) +
  theme_bw() +
  facet_wrap(~author) +
  scale_color_manual(values = method_colors)

