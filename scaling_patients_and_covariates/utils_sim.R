
# install.packages("BiocManager"); BiocManager::install("splatter")
library(splatter)
library(SingleCellExperiment)

simulate_splatter_simple <- function(
    n_cells,
    n_patients,
    n_covariates,
    n_genes = 2000,
    de_prob = 0.10,           # fraction of DE genes between groups
    de_facLoc = 0.1,          # logFC location
    de_facScale = 0.4,        # logFC scale
    batch_facLoc = 0.0,       # patient (batch) effect location (log-scale)
    batch_facScale = 0.2,     # patient (batch) effect scale
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- batches (patients): split cells roughly evenly ---
  batch_cells <- rep(floor(n_cells / n_patients), n_patients)
  batch_cells[seq_len(n_cells - sum(batch_cells))] <- batch_cells[seq_len(n_cells - sum(batch_cells))] + 1
  
  # --- groups from binary covariates (all combinations) ---
  if (n_covariates <= 0) {
    group_mat <- matrix(0, nrow = 1, ncol = 0)
    groups <- 1L
  } else {
    # all binary combinations for K covariates -> 2^K groups
    groups <- 2^n_covariates
    # (cap at a sane number if K is big)
    if (groups > 32) stop("This simple wrapper caps at 32 groups; reduce n_covariates.")
    group_mat <- as.matrix(expand.grid(rep(list(0:1), n_covariates)))
    colnames(group_mat) <- paste0("cov", seq_len(n_covariates))
  }
  
  group_prob <- rep(1 / groups, groups)  # evenly spread groups
  # Build parameters
  p <- newSplatParams(
    nGenes = n_genes,
    batchCells = batch_cells,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale
  )
  
  # Simulate with group-level DE
  sim <- splatSimulateGroups(
    params = p,
    group.prob = group_prob,
    de.prob = de_prob,
    de.downProb = 0.5,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    verbose = FALSE
  )
  
  # Map Splatter's group factor (1..G) to binary covariate columns
  g <- as.integer(colData(sim)$Group)
  if (n_covariates > 0) {
    # recycle the 2^K design rows across groups in order
    cov_df <- as.data.frame(group_mat[g, , drop = FALSE])
  } else {
    cov_df <- data.frame()
  }
  
  # Patient labels from batches
  patient <- factor(colData(sim)$Batch, labels = paste0("P", seq_len(n_patients)))
  
  meta <- data.frame(
    cell = colnames(sim),
    patient = patient,
    group = paste0("G", g),
    cov_df,
    row.names = colnames(sim),
    check.names = FALSE
  )
  
  list(
    counts = as.matrix(counts(sim)),  # genes x cells
    meta = meta,
    sce = sim,                        # keep the SCE for downstream Splatter helpers
    truth = list(
      groups = levels(colData(sim)$Group),
      patients = levels(patient),
      params = p
    )
  )
}