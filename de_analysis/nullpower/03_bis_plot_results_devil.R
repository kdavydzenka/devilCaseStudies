
rm(list = ls())
library(ggrepel)
library(tidyr)
library(magrittr)
library(dplyr)

MY_PALETTE = c(
  "devil (MOM-NB)" = "#099668",
  "devil (MLE-NB)" = "#2A4747"
)

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

METHOD_LEVELS = names(MY_PALETTE)

DEVIL_NAME_MAPPING = c("devil (old overdisp)" = "devil (MLE-NB)",
                       "devil (MOM overdisp)" = "devil (MOM-NB)")

methods = names(MY_PALETTE)

author_mapping = list(
  "bca" = "Reed",
  "hsc" = "Suo", 
  "yazar" = "Yazar",
  "kumar" = "Kumar"
)

# method_patientwise_supp = method_cellwise_supp = METHOD_LEVELS
# method_patientwise_main = method_cellwise_main = METHOD_LEVELS

MAX_GENE = 1000

plot_timing_across_authors = function(a, methods, ratio = NULL) {
  res = readRDS("final_res/results.rds")
  res$author = lapply(res$author, function(a) {
    author_mapping[a]
  }) %>% unlist()
  # if (!is.null(a)) {
  #   res = res %>% dplyr::filter(author == a)
  # }
  
  df = res %>%
    #dplyr::filter(is.pb == FALSE, name %in% methods) %>%
    dplyr::filter(name %in% methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))
  
  if (!is.null(ratio)) {
    df = df %>%
      dplyr::group_by(idx, author) %>%
      dplyr::mutate(Time = Time[name == ratio] / Time)
  }
  
  p = df %>%
    dplyr::mutate(name = factor(name, levels = METHOD_LEVELS)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = author, y=Time, col=name)) +
    geom_boxplot() +
    # geom_point() +
    theme_bw() +
    labs(x = "", col = "") +
    # scale_y_continuous(transform = "log10") +
    scale_color_manual(values = MY_PALETTE)
  
  if (!is.null(ratio)) {
    p = p + labs(y = paste0("Speedup (vs. ", ratio, ")")) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_y_continuous(transform = "log10")
  } else {
    p = p + labs(y = "Time (seconds)")
  }
  #p + coord_flip() + theme(legend.position = "none")
  p + theme(legend.position = "none")
}

get_idx_result = function(author, idx) {
  fp = paste0("results/",author,"_",idx,".rds")
  if (file.exists(fp)) return(readRDS(fp))
  dplyr::tibble()
}

get_result = function(author, is.pb, n_patients = 4, ngenes = 50, cell_index = 1, i = 2, stop_on_error = T) {
  
  param_grid = readRDS(paste0("data/",author,"_param_grid.rds"))
  this = param_grid %>%
    dplyr::mutate(idx = row_number()) %>%
    dplyr::mutate(ng = as.integer(MAX_GENE * prob_de)) %>%
    dplyr::rename(is_pb = is.pb) %>%
    dplyr::filter(n.sample == n_patients, ng == ngenes, int.ct == cell_index, iter == i, is_pb == is.pb)
  
  if (nrow(this) != 1) {
    if (stop_on_error) {
      stop("found either zero or too many (>2) paths")
    } else {
      return(dplyr::tibble())
    }
  }
  
  idx = this$idx
  get_idx_result(author, idx) %>% cbind(this)
}

get_all_results = function(author) {
  param_grid = readRDS(paste0("data/",author,"_param_grid.rds"))
  param_grid$idx = 1:nrow(param_grid)
  df = lapply(param_grid$idx, function(i) {
    get_idx_result(author, i) %>% dplyr::mutate(idx = i)
  }) %>% do.call("bind_rows", .)
  
  for (i in 1:length(names(DEVIL_NAME_MAPPING))) {
    mask = df$name == names(DEVIL_NAME_MAPPING)[i]
    df$name[mask] = DEVIL_NAME_MAPPING[i]
  }
  
  dplyr::left_join(df, param_grid, by = "idx")
}

plot_MCCs_boxplots_across_authors = function(methods) {
  
  res = readRDS("final_res/results.rds")
  
  df = res %>%
    dplyr::mutate(
      author = unlist(lapply(author, function(a) author_mapping[a])),
      name = dplyr::case_when(
        grepl("glmGamPoi", name) ~ "glmGamPoi",
        grepl("Nebula", name) ~ "NEBULA",
        TRUE ~ name
      )
    ) %>%
    dplyr::filter(name %in% methods) %>%
    dplyr::mutate(
      analysis = ifelse(is.pb, "Patient-wise", "Cell-wise"),
      name = factor(name, levels = METHOD_LEVELS)
    )
  
  df %>% dplyr::group_by(name, author, analysis, patients) %>% 
    dplyr::summarise(n = n())
  
  df %>%
    ggplot(aes(x = author, y = MCC, color = name)) +
    geom_boxplot(outlier.shape = NA, width = .7) +
    # geom_jitter(width = .15, alpha = .35, size = 1) +
    ggh4x::facet_nested(~analysis + paste0(patients, " patients")) +
    scale_color_manual(values = MY_PALETTE) +
    theme_bw() +
    labs(
      x = "Dataset",
      y = "MCC",
      color = "Algorithm"
    ) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggpubr::stat_compare_means(
      aes(group = name),
      method = "wilcox.test",
      label = "p.format", show.legend = F, label.y = 1.02
    )
}

# Helper: rename the single facet strip on a ggplot
rename_single_facet <- function(p, new_label) {
  if (inherits(p$facet, "FacetWrap")) {
    fac_var <- names(p$facet$params$facets)
    if (length(fac_var) == 0) fac_var <- names(p$facet$params$rows)
    stopifnot(length(fac_var) == 1)
    lab_map <- setNames(list(c("overall" = new_label)), fac_var)
    p + facet_wrap(reformulate(fac_var),
                   labeller = do.call(labeller, lab_map))
  } else if (inherits(p$facet, "FacetGrid")) {
    # grid: could be in rows or cols; handle both
    rows <- names(p$facet$params$rows)
    cols <- names(p$facet$params$cols)
    lab_args <- list()
    if (length(rows) == 1) lab_args[[rows]] <- c("overall" = new_label)
    if (length(cols) == 1) lab_args[[cols]] <- c("overall" = new_label)
    p + facet_grid(
      rows = if (length(rows) == 1) reformulate(rows) else NULL,
      cols = if (length(cols) == 1) reformulate(cols) else NULL,
      labeller = do.call(labeller, lab_args)
    )
  } else {
    # if there is no facet, just add a title as a fallback
    p + ggtitle(new_label)
  }
}

dir.create("figures/RDS/devil_comparison", recursive = T)

# SUPP FIGURE ####
MCC_box = plot_MCCs_boxplots_across_authors(methods) + MY_THEME
ptiming_ratio = plot_timing_across_authors(NULL, methods, ratio = "devil (MLE-NB)") + MY_THEME

dir.create(paste0("figures/RDS/devil_comparison/"), recursive = T)
saveRDS(MCC_box, paste0("figures/RDS/devil_comparison/MCC_box.rds"))
saveRDS(ptiming_ratio, paste0("figures/RDS/devil_comparison/ptiming_ratio.rds"))



