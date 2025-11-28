
rm(list = ls())
require(tidyverse)
require(patchwork)

# method_cellwise = c("glmGamPoi (cell)", "Devil (base)", "limma (cell)", "Nebula", "Seurat (Cell)", "MAST (cell)", 
#                     "edgeR (cell)", "edgeR (PB)", "limma (Pb)")
# method_patientwise = c("glmGamPoi (cell)", "Devil (mixed)", "limma (cell)", "Nebula", "Seurat (Cell)", 
#                        "MAST (cell)", "edgeR (cell)", "edgeR (PB)", "limma (Pb)")

# FPR test ####
beta <- 0.5
MAX_GENE = 1000
res <- dplyr::tibble()

authors = c("bca", "hsc", "yazar", "kumar")

a = "hsc"
idx = 327

for (a in authors) {
  grid = readRDS(paste0("data/",a,"_param_grid.rds"))
  for (idx in 1:nrow(grid)) {
    this = grid[idx,]
    is.pb = this$is.pb
    
    fp = paste0("results/", a, "_", idx, ".rds")
    fp_dup = paste0("results/", a, "_", idx, "_dupCorr.rds")
    
    if (!file.exists(fp)) {
      print(paste0("Skipping ", a, " - ", idx))
      next
    }
    
    print(paste0(a, " - ", idx))
    
    d = readRDS(fp)
    dtime = readRDS(paste0("timing_results/", a, "_", idx, ".rds")) %>%
      dplyr::select(Time, n.cells, name)
    
    if (file.exists(fp_dup)) {
      d = dplyr::bind_rows(d, readRDS(fp_dup))
      dtime = dplyr::bind_rows(
        dtime, 
        readRDS(paste0("timing_results/", a, "_", idx, "_dupCorr.rds")) %>%
          dplyr::select(Time, n.cells, name)
        )
    }
    
    r <- d %>% 
      na.omit() %>% 
      dplyr::select(lfc, p_val, gene, is_de, name) %>% 
      dplyr::group_by(name) %>% 
      dplyr::mutate(p_val = ifelse(is.na(p_val), 1, p_val)) %>% 
      dplyr::mutate(padj = p.adjust(p_val, "BH")) %>%
      mutate(
        predicted = padj <= 0.05,
        TP = as.numeric(is_de & predicted),    # True Positive
        TN = as.numeric(!is_de & !predicted),  # True Negative
        FP = as.numeric(!is_de & predicted),   # False Positive
        FN = as.numeric(is_de & !predicted),    # False Negative
      ) %>%
      summarise(
        TP = sum(TP),
        TN = sum(TN),
        FP = sum(FP),
        FN = sum(FN),
        numerator = (TP * TN) - (FP * FN),
        denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
        MCC = ifelse(denominator == 0, 0, numerator / denominator), 
        TPR = TP / (TP + FN),
        FDR = FP / (TP + FP), 
        FPR = FP / sum(!is_de)
      ) %>%
      dplyr::select(name, MCC, TPR, FDR, FPR) %>% 
      dplyr::arrange(-MCC)
    
    r = r %>% dplyr::left_join(dtime, by = "name")
    
    r = cbind(r, this) %>% 
      dplyr::mutate(ngenes = as.integer(prob_de * MAX_GENE), author = a, idx = idx) %>% 
      dplyr::rename(ct.index = int.ct, i.iter = iter, patients = n.sample) %>% 
      dplyr::arrange(-MCC)
    
    res <- dplyr::bind_rows(res, r)
  }
}

ggplot(res, aes(x = name, y = MCC, col = name)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(author~is.pb+patients)

ggplot(res, aes(x = name, y = FDR, col = name)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(author~is.pb+patients)

res %>% 
  dplyr::filter(is.pb) %>% 
  ggplot(mapping = aes(x = n.cells, y = MCC, col = paste0(patients, is.pb))) +
  geom_point() +
  facet_grid(author~name, scales = "free") +
  theme_bw() +
  scale_x_continuous(transform = "log10")

res %>% 
  dplyr::filter(is.pb) %>% 
  ggplot(mapping = aes(x = n.cells, y = FDR, col = paste0(patients, is.pb))) +
  geom_point() +
  facet_grid(author~name, scales = "free") +
  theme_bw() +
  scale_x_continuous(transform = "log10")

saveRDS(res, "final_res/results.rds")
