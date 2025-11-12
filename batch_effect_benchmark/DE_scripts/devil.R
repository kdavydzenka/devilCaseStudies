# raw count or pseudobulk data as input processed=rawcount
run_devil <- function(processed, cellinfo, cov = TRUE, Det = FALSE, former.meth = "") {
  stopifnot(is.matrix(processed) || is.data.frame(processed))
  suppressPackageStartupMessages({
    library(edgeR)   # for TMM size factors
    library(devil)
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
  
  # Build design with explicit "is_cell" (the target Group effect)
  # is_cell = 1 for the second level of Group, 0 for the reference
  # (You can change the reference with relevel(cellinfo$Group, ref="..."))
  # is_cell <- as.numeric(cellinfo$Group == levels(cellinfo$Group)[2])
  # design_df <- data.frame(
  #   is_cell = is_cell
  # )
  
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
  
  
  # Fit DEVIL
  fit <- devil::fit_devil(
    count_df,
    design,
    overdispersion = TRUE,
    size_factors = NULL,
    max_iter = 500,
    parallel.cores = 1,
    verbose = TRUE
  )
  
  # Build a contrast vector that tests the "is_cell" coefficient
  # coef_names <- colnames(design_mat)
  # k <- length(coef_names)
  contrast <- c(0,1,0)
  # target_idx <- which(coef_names == "is_cell")
  # if (length(target_idx) != 1) stop("Couldn't find 'is_cell' in design.")
  # contrast[target_idx] <- 1
  
  test.res <- devil::test_de(fit, contrast, max_lfc = 20)
  colnames(test.res) <- c("gene", "pval", "padj", "lfc")
  res <- dplyr::mutate(as.data.frame(test.res), cell_type = levels(cellinfo$Group)[2])
  
  res_name <- paste0(
    ifelse(former.meth == "", "", paste0(former.meth, "+")),
    "devil",
    ifelse(Det, "_Detrate", ""),
    ifelse(cov, "_Cov", "")
  )
  save(res, cellinfo, file = paste0("./", res_name, ".rda"))
  return(res_name)
}
