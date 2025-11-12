# raw count or pseudobulk data as input processed=rawcount
run_nebula <- function(processed, cellinfo, cov = TRUE, Det = FALSE, former.meth = "") {
  stopifnot(is.matrix(processed) || is.data.frame(processed))
  suppressPackageStartupMessages({
    library(nebula)
    library(dplyr)
  })
  
  count_df <- as.matrix(processed)
  storage.mode(count_df) <- "integer"
  
  rownames(cellinfo) <- cellinfo$Cell
  cellinfo <- cellinfo[colnames(count_df), ]
  cellinfo$Group <- factor(cellinfo$Group)
  cellinfo$Batch <- factor(cellinfo$Batch)
  
  # Detection-rate covariate (per-cell % detected), centered & scaled
  cdr <- scale(colMeans(count_df > 0))
  
  # Build predictor matrix with a single primary effect "is_cell"
  if(Det){
    if(cov){
      design <- model.matrix(~Group+cdr+Batch, data = cellinfo)
    }else{
      design <- model.matrix(~Group+cdr, data = cellinfo)
    }
  }else{
    if(cov){
      design <- model.matrix(~Group+Batch, data = cellinfo)
    }else{
      design <- model.matrix(~Group, data = cellinfo)
    }
  }
  pred_mat = design
  
  # Subject/cluster id: use Batch as the grouping unit (fallback to single cluster)
  id_vec <- if (!is.null(cellinfo$Batch)) as.integer(cellinfo$Batch) else rep(1L, nrow(pred_mat))
  
  fit <- nebula::nebula(
    count_df,
    id   = id_vec,
    pred = pred_mat,
    mincp = 0,
    cpc = 0,
    ncore = 1
  )
  
  # Prefer columns tied to the "is_cell" predictor
  sumdf <- fit$summary
  # Typical names: logFC_is_cell, p_is_cell; sometimes *_is_cellTRUE
  lfc_col = "logFC_GroupGroup2"
  p_col = "p_GroupGroup2"
  if (is.na(lfc_col) || is.na(p_col)) stop("Could not locate NEBULA logFC/p columns.")
  
  lfc  <- sumdf[[lfc_col]]
  pval <- sumdf[[p_col]]
  gene <- sumdf[["gene"]]
  padj <- p.adjust(pval, method = "BH")
  
  res <- dplyr::tibble(gene = gene, pval = pval, padj = padj, lfc = lfc,
                       cell_type = levels(cellinfo$Group)[2])
  
  res_name <- paste0(
    ifelse(former.meth == "", "", paste0(former.meth, "+")),
    "nebula",
    ifelse(Det, "_Detrate", ""),
    ifelse(cov, "_Cov", "")
  )
  save(res, cellinfo, file = paste0("./", res_name, ".rda"))
  return(res_name)
}
