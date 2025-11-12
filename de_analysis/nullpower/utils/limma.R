
limma.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
  s <- Sys.time()
  edger.obj <- edgeR::DGEList(counts(sce.pb))
  v <- limma::voom(edger.obj, design)
  vfit <- limma::lmFit(v, design)
  efit <- limma::eBayes(vfit)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()

  beta <- efit$coefficients[,2] * log(2)
  pval <- efit$p.value[,2]
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

limma_voom_cell.mult <- function(count, df){
  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
  })
  df$tx_cell <- factor(df$tx_cell)
  if (length(levels(df$tx_cell)) != 2) stop("tx_cell must have exactly two levels.")
  lev <- levels(df$tx_cell)
  
  mat <- if (ncol(count) == nrow(df)) count else t(count)
  design <- model.matrix(~ 1 + tx_cell, data=df)
  
  s <- Sys.time()
  y <- DGEList(mat)
  y <- calcNormFactors(y)  # TMM
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  # Map to requested columns
  beta <- tt$logFC
  tval <- tt$t
  se   <- beta / tval
  pval <- tt$P.Value
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(beta),
    `Std. Error` = as.numeric(se),
    `t value`    = as.numeric(tval),
    `Pr(>|t|)`   = as.numeric(pval),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  rownames(out) <- rownames(tt)
  out
}

limma.dupCorr.mult <- function(count, df){
  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
  })
  
  # ---- checks ----
  if (!("tx_cell" %in% names(df))) stop("df must contain 'tx_cell'.")
  if (!("id" %in% names(df))) stop("df must contain 'patient' for duplicateCorrelation blocking.")
  df$tx_cell <- factor(df$tx_cell)
  df$id <- factor(df$id)
  if (length(levels(df$tx_cell)) != 2) stop("tx_cell must have exactly two levels.")
  
  # orient matrix: genes x cells
  mat <- if (ncol(count) == nrow(df)) count else t(count)
  
  # design: intercept + treatment
  design <- model.matrix(~ 1 + tx_cell, data = df)
  
  s <- Sys.time()
  # edgeR container + TMM
  y <- DGEList(mat)
  y <- calcNormFactors(y)  # TMM
  # voom transform
  v <- voom(y, design, plot = FALSE)
  # estimate intra-patient correlation
  dupcor <- duplicateCorrelation(v, design = design, block = df$id)
  rho_hat <- dupcor$consensus.correlation
  
  # fit with blocking & correlation
  fit <- lmFit(v, design, block = df$id, correlation = rho_hat)
  fit <- eBayes(fit)  # keep default to match your original
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  # results (coef 2 = tx_cell level 2 vs level 1)
  
  
  beta <- tt$logFC
  tval <- tt$t
  se   <- beta / tval
  pval <- tt$P.Value
  
  out <- tibble::tibble(
    `Estimate`   = as.numeric(beta),
    `Std. Error` = as.numeric(se),
    `t value`    = as.numeric(tval),
    `Pr(>|t|)`   = as.numeric(pval),
    `Time`       = delta_time
  )
  out <- as.matrix(out)
  rownames(out) <- rownames(tt)
  
  # (optional) store rho as an attribute for bookkeeping
  attr(out, "dupCor_consensus") <- rho_hat
  
  out
}
