
rm(list = ls())

source("utils/utils.R")
# source("utils/edgeR.R")
# source("utils/limma.R")
# # source("utils/flash.R")
# source("utils/glmGamPoi.R")
# source("utils/nebula.R")
# source("utils/devil.R")
# source("utils/SeuratWilcox.R")
# source("utils/MAST.R")
# source("utils/deSeq2.R")
library(Seurat)


set.seed(123456)
# --- constants ---------------------------------------------------------------
MAX_GENE <- 1000
exp.cut  <- 0.5
sample.prob = .5
n.sim    <- 5           # keep only iterations inside each task

# --- args & env --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
author <- args[1]
stopifnot(author %in% c("bca","yazar","hsc","kumar"))

if (!(author %in% c('bca', 'yazar', 'hsc', 'kumar'))){
  stop('author must be bca or yazar or hsc or kumar')
}

if (author == "bca") {
  seurat.obj <- readRDS("datasets/bca.seurat.rds")
} else if (author == "yazar") {
  seurat.obj <- readRDS('datasets/yazar.seurat.rds')
} else if (author == 'hsc') {
  seurat.obj <- readRDS('datasets/hsc.rds')
} else if (author == "kumar") {
  seurat.obj <- readRDS('datasets/kumar.rds')
} else {
  stop('author must be bca or yazar')
}

# prepare data
seurat.obj@meta.data
  
# cols <- c('donor_id', 'cell_type')
col.data <- seurat.obj@meta.data
col.data$cell_id <- rownames(col.data)
#cnt <- GetAssayData(object = seurat.obj, slot = "counts")
if (author %in% c("bca", "hsc", "kumar")) {
  cnt <- seurat.obj@assays$RNA$counts
} else if (author %in% c("yazar")) {
  cnt <- seurat.obj@assays$RNA$data
} else  {
  stop("author must be bca or yazar")
}

head(col.data)

# sort cell type by numbers, select top 6
ct.used <- col.data %>%
  group_by(cell_type) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  top_n(6) %>%
  pull(cell_type)

# select donors with more than 50 cells per selected cell types
donor.cells <- col.data %>%
  group_by(donor_id, cell_type, .drop=FALSE) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from=cell_type, values_from=n) %>%
  select(ct.used) %>%
  filter(if_all(ct.used,~.>100))  

donor.cells$cell.sum <- rowSums(donor.cells[,2:7])
donor.used <- donor.cells %>% arrange(-cell.sum) %>% pull(donor_id) %>% as.character()
if (length(donor.used) > 20) {
  donor.used = donor.used[1:20]
}

col.data %>%
  group_by(cell_type) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  top_n(6)

cell.used <- col.data$donor_id %in% donor.used
cnt.used <- cnt[rowMeans(cnt) > 0.1, cell.used]
col.data.used <- col.data[cell.used,]

# Simulations ####
# --- build parameter grid for the ARRAY ---------------------------------------
IS.PB    <- c(TRUE, FALSE)
PROB.DE  <- c(.025, .05, .1)
N.SAMPLE <- c(4, 20)
INT.CT   <- seq_along(ct.used)  # 1..6 by construction
ITERS = 1:n.sim

grid <- expand.grid(
  is.pb    = IS.PB,
  int.ct   = INT.CT,
  prob_de  = PROB.DE,
  n.sample = N.SAMPLE,
  iter = ITERS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
saveRDS(grid, paste0("data/",author,"_param_grid.rds"))

for (idx in 1:nrow(grid)) {
  print(idx)
  
  this = grid[idx,]
  is.pb = this$is.pb
  int.ct = this$int.ct
  prob_de = this$prob_de
  n.sample = this$n.sample
  iter = this$iter
  
  # 1) Focus on one cell type ----------------------------------------------------
  cell.type <- ct.used[[int.ct]]
  bool.ct   <- col.data.used$cell_type == cell.type
  col.data.ct <- col.data.used[bool.ct, , drop = FALSE]
  cnt.ct      <- cnt.used[, bool.ct, drop = FALSE]
  
  # Guard: require enough cells to continue
  # if (ncol(cnt.ct) < 100) next
  
  # 2) Gene filtering: mean + prevalence + HVG -----------------------------------
  # Keep genes with enough mean AND expressed in >=10% of cells (by CPM ~ 0.1)
  libsize <- colSums(cnt.ct)
  cpm     <- edgeR::cpm(cnt.ct)
  preval  <- rowMeans(cpm > 0.1)
  keep_g  <- rowMeans(cnt.ct) > 0.2 & preval >= 0.10
  if (!any(keep_g)) next
  cand    <- which(keep_g)
  
  # Choose HVGs among candidates (still respects a mean floor)
  sel_genes <- sample(cand, MAX_GENE, replace = T)
  #if (length(sel_genes) < 50) next  # safety
  cnt.ct.bm <- cnt.ct[sel_genes, , drop = FALSE]
  
  # 3) Sample donors (prefer without replacement), then balance per donor --------
  all_donors <- unique(col.data.ct$donor_id)
  if (length(all_donors) >= n.sample) {
    donors.sel <- sample(all_donors, size = n.sample, replace = FALSE)
  } else {
    donors.sel <- sample(all_donors, size = n.sample, replace = TRUE)
  }
  
  # Assemble per-donor matrices
  col.data.donor <- lapply(donors.sel, function(id) col.data.ct[col.data.ct$donor_id == id, , drop = FALSE])
  cnt.donor      <- lapply(donors.sel, function(id) cnt.ct.bm[, col.data.ct$donor_id == id, drop = FALSE])
  
  SAFE_COLS <- intersect(c("donor_id", "cell_type", "cell_id"), colnames(col.data.ct))
  
  if (length(SAFE_COLS) == 0) {
    stop("No safe columns found; ensure donor_id/cell_type/cell_id exist in meta.data.")
  }
  
  # Build per-donor parts using SAFE_COLS only; iterate over selection positions
  col.data.parts <- vector("list", length(donors.sel))
  cnt.parts      <- vector("list", length(donors.sel))
  
  for (k in seq_along(donors.sel)) {
    idk <- donors.sel[k]
    rows_k <- which(col.data.ct$donor_id == idk)
    if (length(rows_k) == 0) next
    
    # meta slice with only safe columns
    dfk <- as.data.frame(col.data.ct[rows_k, SAFE_COLS, drop = FALSE])
    dfk$id <- k  # position in the *sampled* donor list (handles replacement)
    rownames(dfk) <- rownames(col.data.ct)[rows_k]
    col.data.parts[[k]] <- dfk
    
    # counts slice matching those cells
    cnt.parts[[k]] <- cnt.ct.bm[, rows_k, drop = FALSE]
  }
  
  # Drop empties (in case a sampled donor had 0 cells after earlier filters)
  nz_idx <- which(vapply(col.data.parts, function(x) !is.null(x) && nrow(x) > 0, logical(1)))
  if (length(nz_idx) == 0) next
  
  col.data.select <- data.table::rbindlist(col.data.parts[nz_idx], use.names = TRUE, fill = TRUE)
  cnt.select      <- do.call(cbind,     cnt.parts[nz_idx])
  
  # ensure colnames(cnt.select) align to rownames(meta)
  colnames(cnt.select) <- rownames(col.data.select)
  
  cnt.select      <- do.call(cbind, cnt.donor)
  colnames(cnt.select) <- rownames(col.data.select)
  
  # Balance: keep up to a target per donor (gentle) to avoid dominance
  per_donor <- col.data.select %>% 
    dplyr::count(id, name = "n")
  target <- max(40L, floor(stats::quantile(per_donor$n, 0.2)))  # >= 40 or 20th pct
  
  idx_keep <- logical(nrow(col.data.select))
  offset <- 0L
  for (k in seq_len(nrow(per_donor))) {
    nk  <- per_donor$n[k]
    sel <- sample.int(nk, size = min(nk, target))
    idx_keep[(offset + sel)] <- TRUE
    offset <- offset + nk
  }
  col.data.select <- col.data.select[idx_keep, , drop = FALSE]
  cnt.select      <- cnt.select[, idx_keep, drop = FALSE]
  
  # Drop zero-count cells (rare after balancing)
  nz <- colSums(cnt.select) > 0
  col.data.select <- col.data.select[nz, , drop = FALSE]
  cnt.select      <- cnt.select[, nz, drop = FALSE]
  
  # Recompute donor counts after balancing
  per_donor <- col.data.select %>% 
    dplyr::count(id, name = "n")
  
  # 4) Treatment assignment (PB vs cell-wise) ------------------------------------
  if (is.pb) {
    # balanced one label per donor
    U <- sort(unique(col.data.select$id))
    n_tx <- floor(length(U) / 2)
    donor_tx <- dplyr::tibble(id = U, tx = 0L)
    donor_tx$tx[sample.int(length(U), n_tx)] <- 1L
    col.data.select <- col.data.select |>
      dplyr::left_join(donor_tx, by = "id") |>
      dplyr::mutate(tx_cell = as.integer(tx))
    col.data.select$tx <- NULL
  } else {
    # per-cell random assignment
    col.data.select$tx_cell <- stats::rbinom(nrow(col.data.select), size = 1L, prob = 0.5)
  }
  
  # for power simulation, cut expression to half (+ force int)
  max_gene <- as.integer(MAX_GENE * prob_de)
  n.gene <- max_gene
  idx.de <- 1:max_gene
  idx.tx <- which(col.data.select$tx_cell == 1)
  downCells(cnt.select, idx.de, idx.tx, exp.cut)
  
  # Filter out zero-count cells after downCells modification
  cell.sums <- colSums(cnt.select)
  bool.nonzero.cells <- cell.sums > 0
  if(sum(bool.nonzero.cells) < ncol(cnt.select)) {
    cnt.select <- cnt.select[, bool.nonzero.cells]
    col.data.select <- col.data.select[bool.nonzero.cells, ]
  }
  
  data = list(cnt = cnt.select, meta = col.data.select)
  
  saveRDS(data, file = paste0("data/", author, "_", idx, ".rds"))
}
