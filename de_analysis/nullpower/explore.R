
rm(list = ls())
source("utils/utils.R")
source("utils/edgeR.R")
source("utils/limma.R")
# source("utils/flash.R")
source("utils/glmGamPoi.R")
source("utils/nebula.R")
source("utils/devil.R")
source("utils/SeuratWilcox.R")
source("utils/MAST.R")
source("utils/deSeq2.R")
# source("utils/devil_robust_sandwich.R")
library(Seurat)
library(ggplot2)
library(tidyverse)

a = "hsc"

res = readRDS("final_res/results.rds")
res %>% 
  dplyr::filter(author == a) %>% 
  ggplot(mapping = aes(x = name, y = MCC, col = name)) +
  geom_boxplot() +
  facet_grid(ct.index~is.pb+patients)

res %>% dplyr::filter(grepl("devil", name), is.pb, patients == 20, author == a) %>% 
  dplyr::select(idx, name, MCC) %>% 
  tidyr::pivot_wider(values_from = MCC, names_from = name) %>%
  dplyr::mutate(delta = `devil (overdisp)` - devil) %>% view()

idx = 115

list.func <- list(
  #deseq2.mult,
  #edger.mult,
  devil.fit, 
  devil.overdisp#,
  #glmgp.cell.mult,
  #nebula.mult
)

method_names = c("deSeq2","edgeR", "devil", "devil.disp", "glmGamPoi", "NEBULA")
method_names = c("devil", "devil.disp")
MAX_GENE = 1000
author = "hsc"
#idx = 47

data = readRDS(paste0("data/",author,"_",idx,".rds"))
params = readRDS(paste0("data/",author,"_param_grid.rds"))[idx,]

cnt.select = data$cnt
col.data.select = data$meta

list.result.method <- list()
timings <- c()
cnt.input <- cnt.select %>% as.matrix()

count = cnt.input
df = col.data.select
is.pb = params$is.pb

df$tx_cell = as.factor(df$tx_cell)
design_matrix <- model.matrix(~1+tx_cell, data = df)
clusters = as.factor(paste0(df$id))

if (is.pb) {
  sf = "psinorm"
} else {
  sf = NULL
}


fit_w_disp <- devil::fit_devil(count, design_matrix, 
                        size_factors=sf, 
                        verbose=F, parallel.cores=1, 
                        init_overdispersion = NULL, overdispersion = T)

fit_wo_disp <- devil::fit_devil(count, design_matrix, 
                               size_factors=sf, 
                               verbose=F, parallel.cores=1, 
                               init_overdispersion = 100, overdispersion = F)

# fit_w_disp$overdispersion = 1 / fit_w_disp$overdispersion
test_w_disp <- devil::test_de(fit_w_disp, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
test_wo_disp <- devil::test_de(fit_wo_disp, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)

# fit_wo_disp$size_factors = fit_w_disp$size_factors = rep(1, length(fit_w_disp$size_factors))
# test_w_disp_sf <- devil::test_de(fit_w_disp, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
# test_wo_disp_sf <- devil::test_de(fit_wo_disp, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)

is_de = 1:nrow(count) <= as.integer(params$prob_de * MAX_GENE)

#test_w_disp_sf %>% 
#test_w_disp %>% 
test_wo_disp %>% 
  dplyr::mutate(is_de = is_de) %>% 
  dplyr::mutate(
    calls = adj_pval <= .05,
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
  ) %>% dplyr::select(MCC, FDR, TP, FP) %>% dplyr::distinct()


plot(test_w_disp$pval, test_wo_disp$pval)
plot(-log10(test_w_disp$adj_pval), -log10(test_wo_disp$adj_pval))


