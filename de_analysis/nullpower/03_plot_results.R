
rm(list = ls())
library(ggrepel)
library(tidyr)
library(magrittr)
library(dplyr)

# NAME_MAPPING = c(
#   "edgeR..Pb." = "edgeR (Pb)",
#   "edgeR..cell." = "edgeR (cell)",
#   "MAST..cell." = "MAST (cell)",
#   "Seurat..cell." = "Seurat (cell)",
#   "limma..Pb." = "limma (Pb)",
#   "limma..cell." = "limma (cell)",
#   "glmGamPoi..Pb." = "glmGamPoi (Pb)",
#   "glmGamPoi..cell." = "glmGamPoi (cell)",
#   "NEBULA" = "Nebula",
#   "devil..base." = "Devil (base)",
#   "devil..mixed." = "Devil (mixed)",
#   "devil..sf.base." = "DevilSF (base)",
#   "devil..sf.mixed." = "DevilSF (mixed)"
# )

DEVIL_NAME_MAPPING = c("devil" = "devil (Poisson)",
                       "devil (new overdisp)" = "devil (I-NB)",
                       "devil (old overdisp)" = "devil (MLE-NB)",
                       "devil (MOM overdisp)" = "devil (MOM-NB)")

MY_PALETTE = c(
  "devil (MOM-NB)" = "#099668",
  "devil (MLE-NB)" = "#2A4747",
  # "devil (I-NB)" = "#6F9283",
  "devil (Poisson)" = "#92AA83",
  "Nebula" = "steelblue",
  "NEBULA" = "steelblue",
  "edgeR" = "#7D629E",
  "edgeR (Pb)" = "#7D629E",
  "limma" = "#B96461",
  "limma (Pb)" = "#B96461",
  "limma (cell)" = "#A02D2C",
  "glmGamPoi (cell)" = "#EAB578",
  "glmGamPoi" = "#EAB578",
  "limmaDupCorr (cell)" = "#8B0000",
  "limmaDupCorr" = "#8B0000",
  "Seurat (cell)" = "lightgray",
  "Seurat (Wilcox)" = "lightgray",
  "Seurat" = "lightgray",
  "MAST (cell)" = "#D8BFD8",
  "MAST" = "#D8BFD8"
)

MY_PALETTE %>% names()

METHOD_LEVELS = c("Devil (base)", "Devil (mixed)", "Devil", "devil (new)", "devil", "devil (overdisp)", "NEBULA", "Nebula", "glmGamPoi (cell)", "glmGamPoi",
                  "MAST (cell)", "MAST", "Seurat (cell)", "Seurat", "Seurat (Wilcox)","limma (cell)", "limmaDupCorr (cell)", "limmaDupCorr",
                  "limma", "limma (Pb)", "edgeR", "edgeR (Pb)")

METHOD_LEVELS = names(MY_PALETTE)

method_patientwise_supp = method_cellwise_supp = c("edgeR (Pb)", "MAST (cell)", "Seurat (cell)", "limma (Pb)", "limma (cell)",
                         "glmGamPoi (cell)", "NEBULA", "limmaDupCorr (cell)", "devil (MLE-NB)", "devil (MOM-NB)")

method_patientwise_main = method_cellwise_main = c("edgeR (Pb)", "MAST (cell)", "limma (Pb)",
                                                   "glmGamPoi (cell)", "NEBULA", "devil (MOM-NB)")

MAX_GENE = 1000

plot_timing = function(a, methods, ratio = NULL) {
  res = readRDS("final_res/results.rds")
  if (!is.null(a)) {
    res = res %>% dplyr::filter(author == a)
  }

  df = res %>%
    #dplyr::filter(is.pb == FALSE, name %in% methods) %>%
    dplyr::filter(name %in% methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  if (!is.null(ratio)) {
    df = df %>%
      dplyr::group_by(idx) %>%
      dplyr::mutate(Time = Time[name == ratio] / Time)
  }

  p = df %>%
    dplyr::mutate(name = factor(name, levels = METHOD_LEVELS)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = name, y=Time, col=name)) +
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
  p + coord_flip() + theme(legend.position = "none")
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

plot_MCCs_boxplots = function(a, cellwise_methods, patientwise_methods) {
  res = readRDS("final_res/results.rds")
  if (!is.null(a)) {
    res = res %>% dplyr::filter(author == a)
  }

  df_cellwise = res %>%
    #na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% cellwise_methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_patientwise = res %>%
    #na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% patientwise_methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_all = dplyr::bind_rows(df_patientwise, df_cellwise)
  df_all %>%
    dplyr::mutate(name = factor(name, levels = METHOD_LEVELS)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = name, y=MCC, col=name)) +
    geom_boxplot() +
    coord_flip() +
    ggh4x::facet_nested(~is.pb+paste0(patients, " patients")) +
    scale_color_manual(values = MY_PALETTE) +
    theme_bw() +
    labs(x = "", col = "") +
    theme(legend.position = "none")
}

plot_FDRs_boxplots = function(a, cellwise_methods, patientwise_methods) {
  res = readRDS("final_res/results.rds")
  if (!is.null(a)) {
    res = res %>% dplyr::filter(author == a)
  }

  df_cellwise = res %>%
    #na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% cellwise_methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_patientwise = res %>%
    #na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% patientwise_methods) %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_all = dplyr::bind_rows(df_patientwise, df_cellwise)
  df_all %>%
    dplyr::mutate(name = factor(name, levels = METHOD_LEVELS)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = name, y=FDR, col=name)) +
    geom_boxplot() +
    coord_flip() +
    ggh4x::facet_nested(~is.pb+paste0(patients, " patients")) +
    scale_color_manual(values = MY_PALETTE) +
    theme_bw() +
    labs(x = "", col = "") +
    theme(legend.position = "none")
}

plot_qq_null_pvals = function(author,
                              cellwise_methods, patientwise_methods, n.points,
                              n.patients = NULL, cell.index = NULL) {

  df = get_all_results(author)

  if (!is.null(n.patients)) df = df %>% dplyr::filter(n.sample == n.patients)
  if (!is.null(cell.index)) df = df %>% dplyr::filter(int.ct == cell.index)

  df = dplyr::bind_rows(
    df %>% dplyr::filter(is.pb, name %in% patientwise_methods),
    df %>% dplyr::filter(!is.pb, name %in% cellwise_methods)
  )

  df = df %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::filter(!is_de)

  if (!is.null(n.points)) {
    df = df %>%
      group_by(name) %>%
      dplyr::sample_n(n.points)
  }

  p <- df %>%
    group_by(name, idx) %>%
    arrange(p_val, .by_group = TRUE) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1),
      x = -log10(expected),
      y = -log10(p_val)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = x, y = y, colour = name)) +
    geom_point(alpha = 0.05, size = 0.1, show.legend = FALSE) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5, colour = "black") +
    scale_color_manual(values = MY_PALETTE) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 3), expand = FALSE) +
    facet_wrap(~is.pb, ncol = 1) +
    labs(
      x = expression(Expected~-log[10](italic(p))),
      y = expression(Observed~-log[10](italic(p)))
    ) +
    labs(color = "Model")
  p
}

plot_power_curve = function(author,
                            methods, is.pwise = TRUE,
                            n.patients = NULL, cell.index = NULL) {

  df = get_all_results(author)
  if (!is.null(n.patients)) df = df %>% dplyr::filter(n.sample == n.patients)
  if (!is.null(cell.index)) df = df %>% dplyr::filter(int.ct == cell.index)

  df = df %>% na.omit() %>% dplyr::distinct()
  df = df %>% dplyr::group_by(name, idx) %>%
    dplyr::mutate(adj_pval = p.adjust(p_val, "BH"))

  df = df %>% dplyr::filter(name %in% methods)
  df = df %>% dplyr::filter(is.pb == is.pwise)

  df = df %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  # TPR vs FDR curve
  # Prep cobra object
  pval = df %>%
    dplyr::select(gene, name, p_val, idx) %>%
    tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))

  padj = df %>%
    #dplyr::mutate(adj_pval = ifelse(abs(lfc) < .5, 1, adj_pval)) %>%
    dplyr::select(gene, name, adj_pval, idx) %>%
    tidyr::pivot_wider(values_from = adj_pval, names_from = name) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))

  truth = df %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, is_de, idx) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()

  rownames(pval) = paste0("Gene", 1:nrow(pval))
  rownames(padj) = paste0("Gene", 1:nrow(pval))
  rownames(truth) = paste0("Gene", 1:nrow(pval))
  truth$feature = rownames(truth)
  pval = pval %>% dplyr::select(!idx)
  truth = truth %>% dplyr::select(!idx)
  padj = padj %>% dplyr::select(!idx)

  library(iCOBRA)
  cobradata = iCOBRA::COBRAData(pval, truth = truth, padj = padj)
  cobraperf <- calculate_performance(cobradata, binary_truth = "is_de")
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2")
  plot_fdrtprcurve(cobraplot, pointsize = 0, linewidth = 1.2) +
    scale_color_manual(values = MY_PALETTE) +
    scale_fill_manual(values = MY_PALETTE) +
    theme_bw()
}

plot_power_curve_faceted = function(author,
                                    cellwise_methods, patientwise_methods,
                                    n.points,
                                    n.patients = NULL, cell.index = NULL) {

  df = get_all_results(author)
  df = df %>% na.omit()
  if (!is.null(n.patients)) df = df %>% dplyr::filter(n.sample == n.patients)
  if (!is.null(cell.index)) df = df %>% dplyr::filter(int.ct == cell.index)

  df = df %>% na.omit() %>% dplyr::distinct()
  df = df %>% dplyr::group_by(name, idx) %>%
    dplyr::mutate(adj_pval = p.adjust(p_val, "BH"))

  df = dplyr::bind_rows(
    df %>% dplyr::filter(is.pb, name %in% patientwise_methods),
    df %>% dplyr::filter(!is.pb, name %in% cellwise_methods)
  )

  df = df %>%
    # dplyr::mutate(name = if_else(grepl("devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  # TPR vs FDR curve
  # Prep cobra object
  pval = df %>%
    dplyr::select(gene, name, p_val, idx) %>%
    tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))

  padj = df %>%
    #dplyr::mutate(adj_pval = ifelse(abs(lfc) < .5, 1, adj_pval)) %>%
    dplyr::select(gene, name, adj_pval, idx) %>%
    tidyr::pivot_wider(values_from = adj_pval, names_from = name) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))

  truth = df %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, is_de, idx, is.pb) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()

  rownames(pval) = paste0("Gene", 1:nrow(pval))
  rownames(padj) = paste0("Gene", 1:nrow(pval))
  rownames(truth) = paste0("Gene", 1:nrow(pval))
  truth$feature = rownames(truth)
  pval = pval %>% dplyr::select(!idx)
  truth = truth %>% dplyr::select(!idx)
  padj = padj %>% dplyr::select(!idx)

  library(iCOBRA)

  if (!is.null(n.points)) {
    idxs = sample(1:nrow(truth), size = n.points)
    pval = pval[idxs,]
    padj = padj[idxs,]
    truth = truth[idxs,]
  }

  cobradata = iCOBRA::COBRAData(pval, truth = truth, padj = padj)
  cobraperf <- calculate_performance(cobradata, binary_truth = "is_de", splv = "is.pb")
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", facetted = TRUE)
  cobraplot@fdrtprcurve = cobraplot@fdrtprcurve %>%
    dplyr::filter(splitval != "overall") %>%
    dplyr::mutate(splitval = ifelse(grepl("Cell-wise", splitval), "Cell-wise", "Patient-wise"))

  plot_fdrtprcurve(cobraplot, plottype = "curve") +
    scale_color_manual(values = MY_PALETTE) +
    scale_fill_manual(values = MY_PALETTE) +
    theme_bw() +
    labs(color = "Model")
}

plot_ks <- function(author,
                    cellwise_methods, patientwise_methods, base_method,
                    n.patients = NULL, cell.index = NULL) {

  res = readRDS("final_res/results.rds")
  if (!is.null(n.patients)) res <- dplyr::filter(res, patients == n.patients)
  if (!is.null(cell.index)) res <- dplyr::filter(res, ct.index == cell.index)

  a <- author

  # CELL-WISE (is.pb == FALSE)
  r_false <- res %>%
    dplyr::filter(author == a) %>%
    dplyr::filter(is.pb == FALSE, name %in% cellwise_methods)


  cellwise_methods = cellwise_methods[cellwise_methods %in% unique(r_false$name)]
  ks_false_pvals <- lapply(cellwise_methods[cellwise_methods != base_method], function(m) {
    pval <- stats::ks.test(
      #dplyr::filter(r_false, name == "devil (base)")$MCC,
      dplyr::filter(r_false, name == base_method)$MCC,
      dplyr::filter(r_false, name == m)$MCC
    )$p.value
    pval_text <- if (pval <= .001) "< .001" else paste0("= ", sprintf("%.3f", pval))
    dplyr::tibble(m = m, pval = pval_text, pval_numeric = pval)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      name = dplyr::case_when(
        grepl("glmGamPoi", m) ~ "glmGamPoi",
        grepl("Nebula", m) ~ "NEBULA",
        TRUE ~ m
      )
    )

  # Canonicalize names for plotting (as in your code)
  r_false_plot <- r_false %>%
    dplyr::mutate(
      name = dplyr::case_when(
        # grepl("devil", name) ~ "devil",
        grepl("glmGamPoi", name) ~ "glmGamPoi",
        grepl("Nebula", name) ~ "NEBULA",
        TRUE ~ name
      ),
      facet = "Cell-wise"
    )

  # Anchor points and annotations for Cell-wise
  curve_false <- r_false %>%
    dplyr::mutate(
      name = dplyr::case_when(
        # grepl("devil", name) ~ "devil",
        grepl("glmGamPoi", name) ~ "glmGamPoi",
        grepl("Nebula", name) ~ "NEBULA",
        TRUE ~ name
      )
    ) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      x = stats::quantile(MCC, 0.75, na.rm = TRUE),
      y = 0.75,
      .groups = "drop"
    ) %>%
    dplyr::left_join(ks_false_pvals, by = "name") %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::mutate(
      label = paste0("p ", pval),
      facet = "Cell-wise"
    )


  # PATIENT-WISE (is.pb == TRUE)
  r_true <- res %>%
    dplyr::filter(author == a) %>%
    dplyr::filter(is.pb == TRUE, name %in% patientwise_methods)

  # KS p-values vs "devil (mixed)" EXACTLY like before
  # ks_true_pvals <- lapply(patientwise_methods[patientwise_methods != "devil (mixed)"], function(m) {
  patientwise_methods = patientwise_methods[patientwise_methods %in% unique(r_true$name)]
  ks_true_pvals <- lapply(patientwise_methods[patientwise_methods != base_method], function(m) {
    pval <- stats::ks.test(
      #dplyr::filter(r_true, name == "devil (mixed)")$MCC,
      dplyr::filter(r_true, name == base_method)$MCC,
      dplyr::filter(r_true, name == m)$MCC
    )$p.value
    pval_text <- if (pval <= .001) "< .001" else paste0("= ", sprintf("%.3f", pval))
    dplyr::tibble(m = m, pval = pval_text, pval_numeric = pval)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      name = dplyr::case_when(
        # grepl("devil", m) ~ "devil",
        grepl("glmGamPoi", m) ~ "glmGamPoi",
        grepl("Nebula", m) ~ "NEBULA",
        TRUE ~ m
      )
    )

  # Canonicalize names for plotting (as in your code)
  r_true_plot <- r_true %>%
    dplyr::mutate(
      name = dplyr::case_when(
        # grepl("devil", name) ~ "devil",
        grepl("glmGamPoi", name) ~ "glmGamPoi",
        grepl("Nebula", name) ~ "NEBULA",
        TRUE ~ name
      ),
      facet = "Patient-wise"
    )

  # Anchor points and annotations for Patient-wise
  curve_true <- r_true %>%
    dplyr::mutate(
      name = dplyr::case_when(
        # grepl("devil", name) ~ "devil",
        grepl("glmGamPoi", name) ~ "glmGamPoi",
        grepl("Nebula", name) ~ "NEBULA",
        TRUE ~ name
      )
    ) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      x = stats::quantile(MCC, 0.75, na.rm = TRUE),
      y = 0.75,
      .groups = "drop"
    ) %>%
    dplyr::left_join(ks_true_pvals, by = "name") %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::mutate(
      label = paste0("p ", pval),
      facet = "Patient-wise"
    )

  plot_df <- dplyr::bind_rows(r_false_plot, r_true_plot)
  ann_df  <- dplyr::bind_rows(curve_false, curve_true)

  p <- ggplot(plot_df, aes(x = MCC, color = name)) +
    stat_ecdf(geom = "line", linewidth = 1) +
    facet_wrap(~ facet, ncol = 2) +
    scale_color_manual(values = MY_PALETTE, name = "Model") +
    labs(x = "MCC", y = "Empirical CDF") +
    theme_bw() +
    theme(legend.position = "left")

  if (nrow(ann_df) > 0) {

    p <- p + ggrepel::geom_label_repel(
      data = ann_df,
      aes(x = x, y = y, label = label, color = name),
      inherit.aes = FALSE,
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      direction = "both",
      show.legend = FALSE
    )
  }

  p
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



# PANELS FOR MAIN AND EXT ####
set.seed(12345)

qq_plot_hsc = plot_qq_null_pvals(author = "hsc", n.points = 10000,
                                 cellwise_methods = method_cellwise_main,
                                 patientwise_methods = method_patientwise_main,
                                 n.patients = 20,
                                 cell.index = NULL)

power_curve = plot_power_curve_faceted(author = "hsc",
                                       cellwise_methods = method_cellwise_main,
                                       patientwise_methods = method_patientwise_main,
                                       n.points = 10000,
                                       n.patients = 20,
                                       cell.index = NULL)

MCC_boxplot = plot_MCCs_boxplots("hsc", method_cellwise_main, method_patientwise_main)
FDR_boxplot = plot_FDRs_boxplots("hsc", method_cellwise_main, method_patientwise_main) +
  geom_hline(yintercept = .05, linetype = "dashed")

ecfd_ks_plot = plot_ks(author = "hsc",
                       cellwise_methods = method_cellwise_main,
                       patientwise_methods = method_patientwise_main,
                       base_method = "devil (MOM-NB)")

ptiming_ratio = plot_timing("hsc", methods = c("NEBULA", "glmGamPoi (cell)", "devil (MOM-NB)"), ratio = "NEBULA")

dir.create("figures/RDS/main", recursive = T)
saveRDS(qq_plot_hsc, "figures/RDS/main/qq_plot.rds")
saveRDS(power_curve, "figures/RDS/main/power_curve.rds")
saveRDS(MCC_boxplot, "figures/RDS/main/MCC_boxplot.rds")
saveRDS(ecfd_ks_plot, "figures/RDS/main/ecfd_ks_plot.rds")
saveRDS(ptiming_ratio, "figures/RDS/main/ptiming_ratio.rds")

# SUPP FIGURE ####
a = "yazar"
for (a in c("hsc", "bca", "yazar", "kumar")) {
  set.seed(12345)

  MCC_box = plot_MCCs_boxplots(a, method_cellwise_supp, method_patientwise_supp)
  FDR_box = plot_FDRs_boxplots(a, method_cellwise_supp, method_patientwise_supp) +
    geom_hline(yintercept = .05, linetype = "dashed")

  qq20 = plot_qq_null_pvals(a, method_cellwise_supp, method_patientwise_supp,
                            n.points = 10000, n.patients = 20, cell.index = NULL)
  qq4 = plot_qq_null_pvals(a, method_cellwise_supp, method_patientwise_supp,
                           n.points = 10000, n.patients = 4, cell.index = NULL)

  power_curve_4 = plot_power_curve_faceted(author = a,
                                           cellwise_methods = method_cellwise_supp,
                                           patientwise_methods = method_patientwise_supp,
                                           n.patients = 4, n.points = 10000,
                                           cell.index = NULL)

  power_curve_20 = plot_power_curve_faceted(author = a,
                                            cellwise_methods = method_cellwise_supp,
                                            patientwise_methods = method_patientwise_supp,
                                            n.patients = 20, n.points = 10000,
                                            cell.index = NULL)

  ecfd_ks_plot = plot_ks(author = a,
                         cellwise_methods = method_cellwise_supp,
                         patientwise_methods = method_patientwise_supp,
                         base_method = "devil (MOM-NB)")

  ptiming = plot_timing(a, method_cellwise_supp, ratio = NULL) +
    scale_y_continuous(transform = "log10")
  ptiming_ratio = plot_timing(a, method_cellwise_supp, ratio = "devil (MOM-NB)")

  dir.create(paste0("figures/RDS/",a), recursive = T)
  saveRDS(MCC_box, paste0("figures/RDS/",a,"/MCC_box.rds"))

  saveRDS(qq20, paste0("figures/RDS/",a,"/qq20.rds"))
  saveRDS(qq4, paste0("figures/RDS/",a,"/qq4.rds"))

  saveRDS(power_curve_4, paste0("figures/RDS/",a,"/power_curve_4.rds"))
  saveRDS(power_curve_20, paste0("figures/RDS/",a,"/power_curve_20.rds"))
  saveRDS(ecfd_ks_plot, paste0("figures/RDS/",a,"/ecfd_ks_plot.rds"))
  saveRDS(ptiming, paste0("figures/RDS/",a,"/ptiming.rds"))
  saveRDS(ptiming_ratio, paste0("figures/RDS/",a,"/ptiming_ratio.rds"))
}
