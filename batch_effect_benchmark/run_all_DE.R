#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(fs)
})

# ---- helpers ---------------------------------------------------------------

run_in_dir <- function(fun, out_dir, ...) {
  old <- getwd(); on.exit(setwd(old), add = TRUE); fs::dir_create(out_dir); setwd(out_dir)
  tryCatch(fun(...), error = function(e){ message("  ✗ ", conditionMessage(e)); NA_character_ })
}

read_counts_matrix <- function(path_counts) {
  # Robust reader: supports first column as gene names or already rownames
  df <- readr::read_tsv(path_counts, col_types = cols())
  # If first column is genes and not numeric, treat it as rownames
  first_is_gene <- !is.numeric(df[[1]]) && !is.logical(df[[1]])
  if (first_is_gene) {
    rn <- df[[1]]
    df <- df[, -1, drop = FALSE]
    mat <- as.matrix(df)
    rownames(mat) <- rn
  } else {
    mat <- as.matrix(df)
    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0("g", seq_len(nrow(mat)))
    }
  }
  storage.mode(mat) <- "integer"
  mat = mat[,!is.na(colMeans(mat))]
  mat
}

read_cellinfo_df <- function(path_cellinfo) {
  ci <- read.delim(path_cellinfo, sep = "\t")
  # Expect columns: Cell, Group, Batch (create if missing)
  if (!"Cell"  %in% names(ci))  stop("cellinfo.txt must contain a 'Cell' column.")
  if (!"Group" %in% names(ci))  stop("cellinfo.txt must contain a 'Group' column.")
  if (!"Batch" %in% names(ci))  ci$Batch <- "Batch1"
  ci
}

safe_source <- function(path) {
  if (fs::file_exists(path)) source(path, local = TRUE)
}

run_one_tool <- function(fun, ..., out_dir) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  fs::dir_create(out_dir)
  setwd(out_dir)
  tryCatch(fun(...),
           error = function(e) {
             message("  ✗ Error: ", conditionMessage(e))
             return(NA_character_)
           })
}

# ---- main ------------------------------------------------------------------

project_root <- fs::path_abs(".")
data_root    <- fs::path(project_root, "DataPrep")
scripts_root <- fs::path(project_root, "DE_scripts")
results_root <- fs::path(project_root, "results")
fs::dir_create(results_root)

# Source tool wrappers
source(fs::path(scripts_root, "pseudobulk.R"))
source(fs::path(scripts_root, "pseudobulk.R"))
source(fs::path(scripts_root, "edgeR.R"))
source(fs::path(scripts_root, "devil.R"))
source(fs::path(scripts_root, "nebula.R"))
source(fs::path(scripts_root, "wilcox.R"))  # optional

# Check required functions exist
need_funs <- c("run_pseudobulk", "run_edgeR", "run_devil", "run_nebula")
missing <- need_funs[!vapply(need_funs, exists, logical(1))]
if (length(missing)) stop("Missing required functions in DE_scripts/: ", paste(missing, collapse = ", "))

# Discover simulation folders: DataPrep/simul*
sim_dirs <- fs::dir_ls(data_root, type = "directory", regexp = "simul")
if (!length(sim_dirs)) stop("No simulation folders found under DataPrep/ (expected 'simul*').")

cat("Found", length(sim_dirs), "simulation folders.\n")

sim_path = sim_dirs[1]
for (sim_path in sim_dirs) {
  sim_name <- fs::path_file(sim_path)
  cat("\n=== ", sim_name, " ===\n", sep = "")
  counts_path   <- fs::path(sim_path, "counts.txt")
  cellinfo_path <- fs::path(sim_path, "cellinfo.txt")
  if (!fs::file_exists(counts_path) || !fs::file_exists(cellinfo_path)) { cat("  Skipping (missing files)\n"); next }
  
  counts   <- read_counts_matrix(counts_path)
  cellinfo <- read_cellinfo_df(cellinfo_path)
  
  # align
  common_cells <- intersect(colnames(counts), cellinfo$Cell)
  if (!length(common_cells)) { cat("  Skipping (no overlapping cells)\n"); next }
  counts <- counts[, common_cells, drop = FALSE]
  rownames(cellinfo) <- cellinfo$Cell
  cellinfo <- cellinfo[common_cells, , drop = FALSE]
  cellinfo$Group <- factor(cellinfo$Group); cellinfo$Batch <- factor(cellinfo$Batch)
  
  sim_out <- fs::path(results_root, sim_name); fs::dir_create(sim_out)
  
  ## 1) PSEUDOBULK (first)
  cat("  • Pseudobulk ... ")
  pb_dir <- fs::path(sim_out, "pseudobulk")
  pb_tag <- run_in_dir(run_pseudobulk, pb_dir,
                       count = counts, cellinfo = cellinfo, former.meth = sim_name)
  if (is.na(pb_tag)) { cat("FAILED\n"); next } else { cat("OK (", pb_tag, ")\n", sep="") }
  
  # load pseudobulk objects saved by run_pseudobulk(): processed, cellinfo
  pb_rda <- fs::path(pb_dir, paste0(pb_tag, ".rda"))
  if (!fs::file_exists(pb_rda)) { cat("  ✗ Missing pseudobulk .rda at ", pb_rda, "\n", sep=""); next }
  load(pb_rda)  # provides: processed (genes x samples), cellinfo (samples metadata)
  pseudobulk_mat <- processed
  pseudobulk_info <- cellinfo_pb
  
  ## 2) edgeR ON PSEUDOBULK
  cat("  • edgeR (pseudobulk) ... ")
  ed_dir <- fs::path(sim_out, "edgeR_pseudobulk")
  ed_tag <- run_in_dir(run_edgeR, ed_dir,
                       processed = pseudobulk_mat,
                       cellinfo  = pseudobulk_info,
                       cov = TRUE, Det = FALSE,  # typically no CDR at pseudobulk level
                       former.meth = sim_name)
  cat(ifelse(is.na(ed_tag), "FAILED\n", paste0("OK (", ed_tag, ")\n")))
  
  ## 3) DEVIL on raw counts (single-cell, cluster-robust)
  cat("  • DEVIL (raw) ... ")
  dv_tag <- run_in_dir(run_devil, fs::path(sim_out, "devil_raw"),
                       processed = counts, cellinfo = cellinfo,
                       cov = TRUE, Det = FALSE, former.meth = sim_name)
  cat(ifelse(is.na(dv_tag), "FAILED\n", paste0("OK (", dv_tag, ")\n")))
  
  ## 4) NEBULA on raw counts (single-cell, mixed model)
  cat("  • NEBULA (raw) ... ")
  nb_tag <- run_in_dir(run_nebula, fs::path(sim_out, "nebula_raw"),
                       processed = counts, cellinfo = cellinfo,
                       cov = TRUE, Det = TRUE, former.meth = sim_name)
  cat(ifelse(is.na(nb_tag), "FAILED\n", paste0("OK (", nb_tag, ")\n")))
}

cat("\nAll done.\n")
