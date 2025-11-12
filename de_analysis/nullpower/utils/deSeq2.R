
library(DESeq2)

deseq2.mult <- function(count, df){
  
  # tx_cell factor with exactly two levels
  df$tx_cell <- factor(df$tx_cell)
  if (length(levels(df$tx_cell)) != 2) stop("tx_cell must have exactly two levels.")
  lev <- levels(df$tx_cell)
  
  # Pseudobulk by id, tx_cell (same as your edgeR path)
  sce.obj <- SingleCellExperiment(list(counts = if (ncol(count) == nrow(df)) count else t(count)),
                                  colData = df)
  sce.pb  <- glmGamPoi::pseudobulk(sce.obj, group_by = vars(id, tx_cell), verbose = FALSE)
  
  
  #dds <- dds[rowSums(counts(dds)) > 0, ]
  
  s <- Sys.time()
  dds <- DESeqDataSetFromMatrix(countData = counts(sce.pb),
                                colData   = as.data.frame(colData(sce.pb)),
                                design    = ~ tx_cell)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("tx_cell", lev[2], lev[1]), tidy = FALSE)
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  beta <- res$log2FoldChange
  se   <- res$lfcSE
  stat <- res$stat
  pval <- res$pvalue
  
  # fill se if NA but beta & stat available
  fill <- is.na(se) & !is.na(beta) & !is.na(stat) & stat != 0
  se[fill] <- beta[fill] / stat[fill]
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(beta),
    `Std. Error` = as.numeric(se),
    `t value`    = as.numeric(stat),
    `Pr(>|t|)`   = as.numeric(pval),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  rownames(out) <- rownames(res)
  out
}
