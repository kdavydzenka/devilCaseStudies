# raw count as input: count = rawcount (genes x cells)
run_pseudobulk <- function(count, cellinfo, former.meth = "") {
  stopifnot(is.matrix(count) || is.data.frame(count))
  count <- as.matrix(count)
  storage.mode(count) <- "integer"
  
  rownames(cellinfo) <- cellinfo$Cell
  cellinfo <- cellinfo[colnames(count), , drop = FALSE]
  cellinfo$Group <- factor(cellinfo$Group)
  cellinfo$Batch <- factor(cellinfo$Batch)
  
  # Sample key = Group_Batch
  sample_key <- interaction(cellinfo$Group, cellinfo$Batch, drop = TRUE, sep = "_")
  
  # Aggregate to pseudobulk: sum counts across cells within each (Group,Batch)
  # count is genes x cells, so aggregate over columns -> transpose twice
  pb <- t(rowsum(t(count), group = sample_key))  # genes x samples
  colnames(pb) <- levels(sample_key)
  
  # Build pseudobulk-level metadata, preserving Group and Batch
  parts <- do.call(rbind, strsplit(colnames(pb), "_", fixed = TRUE))
  cellinfo_pb <- data.frame(
    Cell  = colnames(pb),
    Group = factor(parts[, 1], levels = levels(cellinfo$Group)),
    Batch = factor(parts[, 2], levels = levels(cellinfo$Batch)),
    row.names = colnames(pb)
  )
  
  res <- pb
  processed <- pb
  
  res_name <- paste0(ifelse(former.meth == "", "", paste0(former.meth, "+")), "pseudobulk")
  save(res, processed, cellinfo = cellinfo_pb, file = paste0("./", res_name, ".rda"))
  return(res_name)
}
