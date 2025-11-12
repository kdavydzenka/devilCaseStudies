
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
  filter(if_all(ct.used,~.>50))  

donor.cells$cell.sum <- rowSums(donor.cells[,2:7])
donor.used <- donor.cells %>% arrange(-cell.sum) %>% pull(donor_id)
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
  
  cell.type <- ct.used[[int.ct]]
  bool.ct <- col.data.used$cell_type == cell.type
  col.data.ct <- col.data.used[bool.ct,]
  
  # select genes with sufficient mean (nblmm doesn't work for small mean)
  cnt.ct <- cnt.used[,bool.ct]
  cnt.ct.bm <- cnt.ct[rowMeans(cnt.ct) > 0.1,]
  cnt.ct.bm <- cnt.ct.bm[sample(1:nrow(cnt.ct.bm), MAX_GENE),]
  
  # Filter out cells (columns) with zero total counts AFTER gene filtering
  cell.sums <- colSums(cnt.ct.bm)
  bool.nonzero.cells <- cell.sums > 0
  cnt.ct.bm <- cnt.ct.bm[, bool.nonzero.cells]
  col.data.ct <- col.data.ct[bool.nonzero.cells, ]
  
  # boolean index of cells
  bool.cell.donor <- lapply(
    unique(col.data.ct$donor_id),
    function(id.donor){
      return(col.data.ct$donor_id == id.donor)
    }
  )
  col.data.donor <- lapply(
    bool.cell.donor,
    function(bool.donor){
      return(col.data.ct[bool.donor,])
    }
  )
  cnt.donor <- lapply(
    bool.cell.donor,
    function(bool.donor){
      return(cnt.ct.bm[,bool.donor])
    }
  )
  
  list.result.null <- list()
  list.result.pow <- list()
  
  # sample individuals
  ind.select <- sample.int(
    n=length(unique(col.data.ct$donor_id)),
    size=n.sample,
    replace=TRUE
  )
  
  # construct col.data & cnt matrix
  col.data.select <- data.table::rbindlist(col.data.donor[ind.select], idcol="id")
  cnt.select <- do.call(cbind, cnt.donor[ind.select])
  colnames(cnt.select) <- rownames(col.data.select)
  
  # select cells randomly
  # n.per.donor <- col.data.select %>% group_by(id, donor_id) %>% summarise(n=n()) %>% pull(n)
  # cell.select <- as.logical(unlist(sapply(n.per.donor, rbinom, size=1, p=min(n.per.donor) / n.per.donor)))
  n.per.donor <- col.data.select %>% group_by(id, donor_id) %>% summarise(n=n()) %>% pull(n)
  cell.select <- as.logical(unlist(sapply(n.per.donor, rbinom, size=1, p=sample.prob)))
  
  col.data.select <- col.data.select[cell.select,]
  cnt.select <- cnt.select[,cell.select]
  
  # assign treatment
  col.data.select$tx_cell <- rbinom(n=nrow(col.data.select), size=1, p=0.5)
  if (is.pb){
    # assign treatment label
    n.tx <- as.integer(n.sample/2)
    urn <- c(rep(1,n.tx), rep(0,n.sample-n.tx))
    tx.ind <- sample(x=urn, size=n.sample, replace=FALSE)
    
    # assign tx to cells
    cell.per.ind <- col.data.select %>%
      group_by(id) %>%
      summarise(n=n())
    
    col.data.select$tx_cell <- rep(tx.ind, times=cell.per.ind$n)
  }
  
  col.data.select %>% dplyr::select(id, tx_cell) %>% table()
  
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

