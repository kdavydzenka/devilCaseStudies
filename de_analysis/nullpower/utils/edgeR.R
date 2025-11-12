
require(glmGamPoi)
require(SingleCellExperiment)
require(edgeR)

edger.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
  s <- Sys.time()
  edger.obj <- edgeR::DGEList(counts(sce.pb))
  edger.obj <- edgeR::estimateDisp(edger.obj, design)
  fit <- edgeR::glmQLFit(y=edger.obj, design=design)
  test <- edgeR::glmTreat(fit, coef=2)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  beta <- test$coefficients[,2]
  pval <- test$table[,'PValue']
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

edger_cell.mult <- function(count, df){
  suppressPackageStartupMessages(library(edgeR))
  df$tx_cell <- factor(df$tx_cell)
  if (length(levels(df$tx_cell)) != 2) stop("tx_cell must have exactly two levels.")
  
  mat <- if (ncol(count) == nrow(df)) count else t(count)
  design <- model.matrix(~ 1 + tx_cell, data = df)
  
  s = Sys.time()
  y <- DGEList(mat)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  test <- glmQLFTest(fit, coef = 2)  # no lfc thresholding; plain QL F-test
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  tab  <- test$table
  beta <- tab$logFC
  # For 1 df numerator tests, F = t^2; recover signed t using logFC sign
  tval <- sign(beta) * sqrt(tab$F)
  pval <- tab$PValue
  se   <- beta / tval
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(beta),
    `Std. Error` = as.numeric(se),
    `t value`    = as.numeric(tval),
    `Pr(>|t|)`   = as.numeric(pval),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  rownames(out) <- rownames(tab)
  out
}
