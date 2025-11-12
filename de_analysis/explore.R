

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

MAX_GENE = 1000

get_result = function(author, is.pb, n_patients = 4, ngenes = 50, cell_index = 1, i = 2, stop_on_error = T) {
  
  param_grid = readRDS("nullpower/data/param_grid.rds")
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
  
  # Read DE res
  if (is.pb) {
    null_folder = "nullpower/null_subject/"
    pow_folder = "nullpower/pow_subject/"
  } else {
    null_folder = "nullpower/null_cell/"
    pow_folder = "nullpower/pow_cell//"
  }
  
  idx = this$idx
  
  if (!file.exists(file.path(null_folder, paste0(author, "_", idx, ".csv")))) {
    if (stop_on_error) {
      stop("file not found")  
    } else {
      return(dplyr::tibble())
    }
  }
  
  d_non_de = read.delim(file.path(null_folder, paste0(author, "_", idx, ".csv")), sep = ",") %>% 
    tidyr::pivot_longer(!X) %>% 
    dplyr::rename(gene = X, p_val = value) %>% 
    dplyr::mutate(is_de = FALSE)
  d_de = read.delim(file.path(pow_folder, paste0(author, "_", idx, ".csv")), sep = ",") %>% 
    tidyr::pivot_longer(!X) %>% 
    dplyr::rename(gene = X, p_val = value) %>% 
    dplyr::mutate(is_de = TRUE)
  
  # Merge results and compute adj.pvalues
  dplyr::bind_rows(d_non_de, d_de) %>%
    dplyr::group_by(name) %>% 
    dplyr::mutate(adj_pval = p.adjust(p_val, method = "BH")) %>% 
    dplyr::group_by(name) %>% 
    dplyr::mutate(author = author, is.pb = is.pb, n_patients = n_patients, ngenes = ngenes, cell_index = cell_index, iter = i) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(idx = idx)
}

is.pb = T
author = "kumar"
ct.index = 2
n.patients = 20
iters = c(1:5)
ngenes = c(5, 25, 50)
ct.indexes = 2

df = lapply(iters, function(i) {
  lapply(ngenes, function(ng) {
    lapply(ct.indexes, function(ct) {
      lapply(n.patients, function(np) {
        print(paste(is.pb, np, ng, ct, i, sep = " - "))
        get_result(author, is.pb = is.pb, n_patients = np, ngenes = ng, cell_index = ct, i = i, stop_on_error = F)
      }) %>% do.call("bind_rows", .)
    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .)
}) %>% do.call("bind_rows", .)

df$name = NAME_MAPPING[df$name]
method_patientwise = c(method_patientwise, "DevilSF (mixed)")

if (is.pb) df = df %>% dplyr::filter(name %in% method_patientwise)
if (!is.pb) df = df %>% dplyr::filter(name %in% method_cellwise)

idxs = df$idx %>% unique()

idx = idxs[5]

df_data = lapply(idxs, function(idx) {
  x = readRDS(paste0("nullpower/data/", author, "_", idx, ".rds"))
  
  cell_per_types = x$meta %>% 
    dplyr::select(id, tx_cell) %>% table() %>% colSums()
  
  m = x$meta %>% 
    dplyr::mutate(donor_id = as.character(donor_id)) %>% 
    dplyr::select(donor_id, tx_cell) %>% table()
  mzeros = (m == 0) %>% colSums()
  
  print(idx)
  # print(x$meta %>% 
  #         dplyr::mutate(donor_id = as.character(donor_id)) %>% 
  #         dplyr::select(donor_id, tx_cell) %>% table())
  
  m = x$meta %>% 
    dplyr::mutate(donor_id = as.character(donor_id)) %>% 
    dplyr::select(donor_id, id, tx_cell) %>% table()
  
  print(rowSums(m[,,1] > 0) %>% table())
  print(rowSums(m[,,2] > 0) %>% table())
  
  unique_donors = x$meta$donor_id %>% unique() %>% length()
  unique_id = x$meta$id %>% unique() %>% length()
  
  dplyr::tibble(idx = idx, cell0 = cell_per_types[1], cell1 = cell_per_types[2], unique_donors = unique_donors, unique_id, zero0 = mzeros[1], zero1 = mzeros[2])  
}) %>% do.call("bind_rows", .)



df %>% 
  dplyr::mutate(pred = adj_pval <= .05) %>% 
  dplyr::mutate(TP = pred & is_de, FP = pred & !is_de) %>% 
  dplyr::group_by(name, idx) %>% 
  dplyr::summarise(TP = sum(TP), FP = sum(FP), P = sum(is_de), N = sum(!is_de)) %>% 
  dplyr::mutate(TPR = TP / P, FDR = FP / (TP + FP)) %>% 
  dplyr::left_join(df_data) %>% view()






df_long = lapply(idxs, function(idx) {
  x = readRDS(paste0("nullpower/data/", author, "_", idx, ".rds"))
  
  cell_per_types = x$meta %>% 
    dplyr::select(id, tx_cell) %>% table() %>% colSums()
  
  m = x$meta %>% 
    dplyr::mutate(donor_id = as.character(donor_id)) %>% 
    dplyr::select(donor_id, tx_cell) %>% table()
  mzeros = (m == 0) %>% colSums()
  
  x$meta %>% 
    dplyr::mutate(donor_id = as.character(donor_id)) %>% 
    dplyr::select(donor_id, tx_cell) %>% 
    dplyr::group_by(donor_id, tx_cell) %>% 
    dplyr::summarise(n = n()) %>% 
    tidyr::pivot_wider(values_from = n, names_from = tx_cell, values_fill = 0) %>% 
    dplyr::mutate(time = idx)
}) %>% do.call(bind_rows, .)


X <- df_long %>%
  ungroup() %>% 
  select(`0`,`1`) %>%
  as.matrix()

pca <- prcomp(X, scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col = as.factor(df_long$time))
