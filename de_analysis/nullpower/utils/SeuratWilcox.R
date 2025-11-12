
library(Seurat)

seurat_wilcox.mult <- function(count, df){
  
  clusters = as.factor(df$id)
  df$tx_cell <- factor(df$tx_cell)
  if (length(levels(df$tx_cell)) != 2) stop("tx_cell must have exactly two levels.")
  lev <- levels(df$tx_cell)
  
  # Orient counts to genes x cells
  mat <- if (ncol(count) == nrow(df)) count else t(count)
  
  so <- CreateSeuratObject(counts = mat, meta.data = df)
  so$tx_cell <- factor(so$tx_cell, levels = lev)
  Idents(so) <- "tx_cell"
  so <- NormalizeData(so, verbose = FALSE)
  
  s <- Sys.time()
  mrk <- FindMarkers(
    so,
    ident.1 = lev[2],
    ident.2 = lev[1],
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  est  <- mrk$avg_log2FC
  pval <- mrk$p_val
  # reconstruct z-like stat and SE to match your edgeR output shape
  tval <- qnorm(1 - pval/2) * sign(est)
  se   <- est / tval
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(est),
    `Std. Error` = as.numeric(se),
    `t value`    = as.numeric(tval),
    `Pr(>|t|)`   = as.numeric(pval),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  rownames(out) <- rownames(mrk)
  out
}
