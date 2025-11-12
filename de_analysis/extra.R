
author = "hsc"
is.pb = F
n_patients = 20
ngenes = 50
cell_index = 1
iter = 5

df = lapply(1:6, function(ct) {
  lapply(1:5, function(iter) {
    get_result(author, is.pb, n_patients, ngenes, ct, iter)
  }) %>% do.call("bind_rows", .)
}) %>% do.call("bind_rows", .)

#df = get_result(author, is.pb, n_patients, ngenes, cell_index, iter)

df
thresh_col = "adj_pval"

tools = unique(d$name)
tpr_fdr_df = lapply(tools, function(tool) {
  dt = df %>% dplyr::filter(name == tool) %>% 
    dplyr::select(X, status, name, adj_pval, p_val)
  
  thresh = sort(unlist(c(dt[,thresh_col]))) %>% 
    unique()
  
  lapply(thresh, function(t) {
    dt$pred = dt[,thresh_col] <= t
    dt$gt = (dt$status == 1)
    P = sum(dt$gt)
    N = sum(!dt$gt)
    TP = sum(dt$pred & dt$gt)
    FP = sum(dt$pred & !dt$gt)
    FDR = FP / (TP + FP)
    TPR = TP / P
    dplyr::tibble(thresh = t, FDR = FDR, TPR = TPR)  
  }) %>% do.call("bind_rows", .) %>% dplyr::mutate(name = tool)  
}) %>% do.call("bind_rows", .)

eps <- 1e-6
tpr_fdr_df %>%
  ggplot(aes(x = FDR, y = TPR, col = name)) +
  geom_line()



# Build PR points per tool by cumulative counts as threshold moves
pr_points <- d %>%
  group_by(name) %>%
  mutate(
    tp = cumsum(status == 1),
    fp = cumsum(status == 0),
    TPR = tp / (tp + fp),
    FDR    = tp / P
  ) %>%
  # Optional: add an initial anchor at recall=0, precision=1 for nicer paths
  bind_rows(
    d %>%
      distinct(name) %>%
      mutate(tp = 0L, fp = 0L, precision = 1, recall = 0, adj_pval = NA_real_)
  ) %>%
  arrange(name, recall, adj_pval) %>%
  ungroup()

# Plot
ggplot(pr_points, aes(x = recall, y = precision, color = name)) +
  geom_path(linewidth = 1) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = prevalence, linetype = "dashed") +
  #scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Precision–Recall curves by tool",
    x = "Recall (TP / all positives)",
    y = "Precision (TP / (TP+FP))",
    color = "Tool",
    subtitle = sprintf("Dashed line = class prevalence = %.3f", prevalence)
  ) +
  theme_minimal(base_size = 12)


library(dplyr)
library(ggplot2)
library(tidyr)