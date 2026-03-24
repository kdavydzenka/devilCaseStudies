
library(MAST)
library(Matrix)
library(edgeR)

run_MAST_DE <- function(counts, df, tx_col = "tx_cell") {
  
  library(MAST)
  library(dplyr)
  
  # checks
  stopifnot(ncol(counts) == nrow(df))
  stopifnot(tx_col %in% colnames(df))
  
  # expression matrix (cells x genes)
  expr <- as.matrix(counts)
  
  # ensure binary factor
  df[[tx_col]] <- factor(df[[tx_col]])
  
  # feature metadata
  fdata <- data.frame(primerid = rownames(expr))
  
  # build SingleCellAssay
  sca <- FromMatrix(
    exprsArray = log1p(expr),
    cData = df,
    fData = fdata
  )
  
  # cell detection rate
  cdr <- colSums(expr > 0)
  colData(sca)$cdr <- scale(cdr)
  
  # fit hurdle model
  s <- Sys.time()
  formula_str <- paste0("~ cdr + ", tx_col)
  #formula_str <- paste0("~ ", tx_col)
  zlmFit <- zlm(as.formula(formula_str), sca)
  
  # likelihood ratio test
  lrt <- lrTest(zlmFit, "tx_cell")
  e = Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  # p-values
  pvals <- lrt %>%
    as.data.frame() %>%
    dplyr::select("hurdle.Pr(>Chisq)")
  
  out <- tibble::tibble(
    `Estimate`   = NA,
    `Std. Error` = NA,
    `t value`    = NA,
    `Pr(>|t|)`   = pvals[,1],
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  out
}


mast.mult <- function(count, df){
  
  # Ensure tx_cell is a 2-level factor; use df directly (no helpers)
  tx <- factor(df$tx_cell)
  if (length(levels(tx)) != 2) stop("tx_cell must have exactly two levels.")
  lev <- levels(tx)
  
  # Orient counts so columns = cells (must match nrow(df))
  mat <- if (ncol(count) == nrow(df)) count else t(count)
  
  # CPM per cell then log2(1+CPM) for MAST
  lib_size <- colSums(mat)
  cpm <- sweep(mat, 2, pmax(1, lib_size), "/") * 1e6
  logcpm1 <- log2(cpm + 1)
  
  sca <- FromMatrix(exprsArray = as.matrix(logcpm1),
                    cData = data.frame(tx_cell = tx, row.names = colnames(logcpm1)))
  
  # cellular detection rate
  colData(sca)$cngeneson <- scale(colSums(assay(sca) > 0))
  
  s <- Sys.time()
  z <- zlm(~ tx_cell + cngeneson, sca = sca)
  sm <- summary(z, doLRT = "tx_cell1")
  summaryDt <- sm$datatable
  fcHurdle <- merge(summaryDt[contrast=='tx_cell1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='tx_cell1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(fcHurdle$coef),
    `Std. Error` = NA,
    `t value`    = NA,
    `Pr(>|t|)`   = as.numeric(fcHurdle$`Pr(>Chisq)`),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  out
}
