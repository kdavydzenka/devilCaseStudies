
NAME_MAPPING = c(
  "edgeR..Pb." = "edgeR (Pb)",
  "edgeR..cell." = "edgeR (cell)",
  "MAST..cell." = "MAST (cell)",
  "Seurat..cell." = "Seurat (cell)",
  "limma..Pb." = "limma (Pb)",
  "limma..cell." = "limma (cell)",
  "glmGamPoi..Pb." = "glmGamPoi (Pb)",
  "glmGamPoi..cell." = "glmGamPoi (cell)",
  "NEBULA" = "Nebula",
  "devil..base." = "Devil (base)",
  "devil..mixed." = "Devil (mixed)",
  "devil..sf.base." = "DevilSF (base)",
  "devil..sf.mixed." = "DevilSF (mixed)"
)

method_cellwise = c("edgeR (Pb)", "MAST (cell)", "Seurat (cell)", "limma (Pb)", 
                    "glmGamPoi (cell)", "NEBULA", "devil (base)")
method_patientwise = c("edgeR (Pb)", "MAST (cell)", "Seurat (cell)", "limma (Pb)", 
                       "glmGamPoi (cell)", "NEBULA", "devil (mixed)", "devil (sf.mixed)")


MAX_GENE = 1000

get_idx_result = function(author, idx) {
  fp = paste0("nullpower/results/",author,"_",idx,".rds")
  if (file.exists(fp)) return(readRDS(fp))
  dplyr::tibble()
}

get_result = function(author, is.pb, n_patients = 4, ngenes = 50, cell_index = 1, i = 2, stop_on_error = T) {
  
  param_grid = readRDS(paste0("nullpower/data/",author,"_param_grid.rds"))
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
  param_grid = readRDS(paste0("nullpower/data/",author,"_param_grid.rds"))
  param_grid$idx = 1:nrow(param_grid)
  df = lapply(param_grid$idx, function(i) {
    get_idx_result(author, i) %>% dplyr::mutate(idx = i)
  }) %>% do.call("bind_rows", .)
  dplyr::left_join(df, param_grid, by = "idx")
}


is.cellwise = F
n.patients = 20
author = "yazar"
cell.index = 1

plot_null_pvals = function(author, is.cellwise, n.patients, cell.index) {
  df = get_all_results(author)
  df = df %>% dplyr::filter(is.pb == !is.cellwise, n.sample == n.patients, int.ct == cell.index)
  
  if (!is.cellwise) df = df %>% dplyr::filter(name %in% method_patientwise)
  if (is.cellwise) df = df %>% dplyr::filter(name %in% method_cellwise)
  
  df %>% 
    dplyr::filter(!is_de) %>% 
    group_by(name, idx) %>%
    arrange(p_val, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1)
    ) %>% 
    ggplot(mapping = aes(x = -log10(expected), y =-log10(p_val), col = name)) +
    # ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
    geom_smooth(method = "gam") +
    # geom_point(size = .5) +
    theme_bw() +
    ggsci::scale_color_nejm() +
    #scale_x_continuous(limits = c(0, 5)) +
    # scale_y_continuous(limits = c(0, 5)) +
    geom_abline(slope = 1, intercept = 0)
  
  # TPR vs FDR curve
  df = df %>% na.omit() %>% dplyr::distinct()
  df = df %>% dplyr::group_by(name, idx) %>% 
    dplyr::mutate(adj_pval = p.adjust(p_val, "BH"))
  
  # Prep cobra object
  pval = df %>%
    dplyr::select(gene, name, p_val, idx) %>%
    tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(!gene) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))
  
  padj = df %>%
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
  #cobradata <- calculate_adjp(cobradata, method = "BH")
  cobraperf <- calculate_performance(cobradata, binary_truth = "is_de")
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2")
  plot_fdrtprcurve(cobraplot) +
    ggsci::scale_color_bmj()
  plot_roc(cobraplot) +
    ggsci::scale_color_bmj()
  plot_fdrnbrcurve(cobraplot) +
    ggsci::scale_color_bmj()
}


df %>% 
  ungroup() %>% 
  mutate(
    predicted = adj_pval <= 0.05 & lfc > .5,
    TP = as.numeric(is_de & predicted),    # True Positive
    TN = as.numeric(!is_de & !predicted),  # True Negative
    FP = as.numeric(!is_de & predicted),   # False Positive
    FN = as.numeric(is_de & !predicted)    # False Negative
  ) %>%
  dplyr::group_by(name) %>% 
  summarise(
    TP = sum(TP),
    TN = sum(TN),
    FP = sum(FP),
    FN = sum(FN),
    numerator = (TP * TN) - (FP * FN),
    denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    MCC = ifelse(denominator == 0, 0, numerator / denominator), 
    TPR = TP / (TP + FN),
    FDR = FP / (TP + FP)
  ) %>%
  dplyr::select(name, MCC, TPR, FDR)

df %>% 
  ungroup() %>% 
  mutate(
    predicted = adj_pval <= 0.05,
    TP = as.numeric(is_de & predicted),    # True Positive
    TN = as.numeric(!is_de & !predicted),  # True Negative
    FP = as.numeric(!is_de & predicted),   # False Positive
    FN = as.numeric(is_de & !predicted)    # False Negative
  ) %>% 
  dplyr::filter(FP == 1) %>% 
  ggplot(mapping = aes(x = name, y = lfc)) +
  geom_boxplot()
  


df %>% 
  dplyr::filter(name %in% c("NEBULA", "devill (mixed)")) %>% 
  dplyr::select(lfc, gene, name) %>%
  tidyr::pivot_wider(values_from = lfc, names_from = name) %>% 
  ggplot(mapping = aes())
