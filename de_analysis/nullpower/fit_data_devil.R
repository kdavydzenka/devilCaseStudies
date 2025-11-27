
rm(list = ls())

source("utils/utils.R")
source("utils/devil.R")
source("utils/glmGamPoi.R")
library(Seurat)
library(ggplot2)
library(tidyverse)

list.func <- list(
  glmgp.cell.mult,  
  devil.fit, 
  devil.overdisp.new,
  devil.overdisp.old,
  devil.overdisp.MOM
)

method_names = c(
  "glmGamPoi",
  "devil", "devil (new overdisp)", 
  "devil (old overdisp)", "devil (MOM overdisp)", "devil (SC overdisp)"
  #, "devil (base)", "devil (mixed)", "devil (pure)"
)

set.seed(12345)
# --- constants ---------------------------------------------------------------
MAX_GENE <- 1000

# --- args & env --------------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) >= 1)
# author <- args[1]
# stopifnot(author %in% c("bca","yazar","hsc","kumar"))

# author = "yazar"
# idx = 112
author = "hsc"
idx = 39
idx = 1

# idx  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
param_grid = readRDS(paste0("data/",author,"_param_grid.rds")) %>%
  dplyr::mutate(i = row_number())
this = param_grid[idx,]

data = readRDS(paste0("data/",author,"_",idx,".rds"))

message('start regression')

cnt.select = data$cnt
col.data.select = data$meta
rm(data)
dim(cnt.select)
list.result.method <- list()
timings <- c()
cnt.input <- cnt.select %>% as.matrix()
rownames(cnt.input) = paste0("Gene", 1:nrow(cnt.input))

s <- Sys.time()
df.result = lapply(1:length(list.func), function(int.test) {
  print(int.test)
  
  is.pb <<- this$is.pb
  df_res = list.func[[int.test]](
    cnt.input, 
    col.data.select
  ) %>% as.tibble()
  df_res$gene = rownames(cnt.input)
  df_res$is_de = FALSE
  df_res$is_de[1:as.integer(this$prob_de * MAX_GENE)] = T
  
  df_res %>%
    dplyr::select(Estimate, `Pr(>|t|)`, Time, gene, is_de) %>% 
    dplyr::rename(lfc = Estimate, p_val = `Pr(>|t|)`) %>% 
    dplyr::mutate(name = method_names[int.test])
}) %>% do.call("bind_rows", .)
e = Sys.time()

max_gene <- as.integer(MAX_GENE * this$prob_de)
n.gene <- max_gene

timing_df = df.result %>% 
  dplyr::select(name, Time) %>% 
  dplyr::bind_cols(this) %>% distinct() %>% 
  dplyr::mutate(n.cells=dim(cnt.select)[2])

tibb = df.result

df = tibb %>%
  dplyr::group_by(name) %>%
  na.omit() %>% 
  dplyr::mutate(p.adj = p.adjust(p_val, method = "BH"))

df %>%
  dplyr::filter(is_de == FALSE) %>%
  ggplot(mapping = aes(x = p_val)) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black") +
  facet_wrap(~name, scales = "free_y") +
  labs(title = "P-value Distribution (Non-DE genes)",
       subtitle = "Should be approximately uniform",
       x = "P-value", y = "Count") +
  theme_bw()

df %>%
  dplyr::filter(is_de == FALSE) %>%
  group_by(name) %>%
  arrange(p_val, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    expected = rank / (n() + 1)
  ) %>%
  ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

df %>%
  dplyr::filter(is_de == FALSE) %>%
  group_by(name) %>%
  arrange(p_val, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    expected = rank / (n() + 1)
  ) %>%
  ggplot(mapping = aes(x = -log10(expected), y = -log10(p_val), col = name)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  #ggsci::scale_color_bmj() +
  scale_x_continuous(limits = c(0, 4)) +
  scale_y_continuous(limits = c(0, 4)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "-log10(Expected p-value)",
       y = "-log10(Observed p-value)",
       title = "Q-Q Plot: P-value Uniformity (Non-DE genes)",
       color = "Method") +
  theme(legend.position = "bottom")


df %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(
    calls = p.adj <= .05,
    TP = sum(calls & is_de),
    FP = sum(calls & !is_de),
    FN = sum(!calls & is_de),
    precision = TP / (TP + FP + 1e-9),
    recall    = TP / (TP + FN + 1e-9),
    FDR       = FP / (TP + FP + 1e-9),
    N = length(calls),
    S = (TP + FN) / N,
    P = (TP + FP) / N,
    num = TP / N - S * P,
    den = sqrt(P * S * (1 - S) * (1 - P)),
    MCC = num / den,
    F1 = 2 * TP / (2 * TP + FP + FN)
  ) %>% dplyr::select(name, TP, FP, FN, MCC, F1, FDR) %>% distinct() %>%
  dplyr::arrange(-MCC)

df$X = df$gene
df$adj_pval = df$p.adj
df$status = df$is_de

pval = df %>%
  dplyr::select(X, name, p_val) %>%
  tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
  dplyr::arrange(X) %>%
  dplyr::select(!X) %>%
  as.data.frame()
rownames(pval) = paste0("Gene", 1:nrow(pval))

padj = df %>%
  dplyr::select(X, name, adj_pval) %>%
  tidyr::pivot_wider(values_from = adj_pval, names_from = name) %>%
  dplyr::arrange(X) %>%
  dplyr::select(!X) %>%
  as.data.frame()
rownames(pval) = paste0("Gene", 1:nrow(pval))

truth = df %>%
  ungroup() %>% 
  dplyr::select(X, status) %>%
  dplyr::distinct() %>%
  dplyr::arrange(X) %>%
  dplyr::select(!X) %>%
  as.data.frame()

df %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(TP = sum(adj_pval <= .05 & status == 1))

df %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(FP = sum(adj_pval <= .05 & status == 0))

rownames(pval) = paste0("Gene", 1:nrow(pval))
rownames(padj) = paste0("Gene", 1:nrow(pval))
rownames(truth) = paste0("Gene", 1:nrow(pval))
truth$feature = rownames(truth)

library(iCOBRA)
cobradata = iCOBRA::COBRAData(pval, truth = truth, padj = padj)
#cobradata <- calculate_adjp(cobradata, method = "BH")
cobraperf <- calculate_performance(cobradata, binary_truth = "status")
cobraplot = iCOBRA::prepare_data_for_plot(cobraperf, colorscheme = "Dark2")

auroc_per_method <- cobraperf@roc %>%
  arrange(method, FPR) %>% 
  group_by(method) %>%
  summarise(
    AUROC = sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2),
    .groups = "drop"
  ) %>% 
  dplyr::arrange(-AUROC)

plot_fdrtprcurve(cobraplot) +
  ggsci::scale_color_bmj()

plot_roc(cobraplot)

d_tpr_fdr <- df %>%
  group_by(name) %>%
  arrange(adj_pval, .by_group = TRUE) %>%
  mutate(
    is_de = as.integer(is_de),
    P     = sum(is_de),
    TP    = cumsum(is_de),
    FP    = cumsum(1L - is_de),
    TPR   = ifelse(P > 0, TP / P, NA_real_),
    denom = TP + FP,
    FDR   = ifelse(denom > 0, FP / denom, NA_real_)
  ) %>%
  # collapse duplicates so each unique threshold (p_adj) has the *last* cumulative stats
  group_by(name, adj_pval) %>%
  summarise(TP = last(TP), FP = last(FP), P = last(P),
            TPR = last(TPR), FDR = last(FDR), .groups = "drop") %>%
  dplyr::rename(p = adj_pval)

