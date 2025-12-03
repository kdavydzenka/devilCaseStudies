#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(fs)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
})

# ----------------------- CONFIG ---------------------------------------------

DATAPREP_ROOT <- "DataPrep"
RESULTS_ROOT  <- "results"

# Methods you currently have available from your pipeline
# pattern is a glob/pcre that matches the saved .rda inside the method folder
METHOD_SPECS <- list(
  Pseudobulk_edgeR = list(
    folder = "edgeR_pseudobulk",
    pattern = "*\\.rda$",
    # Direction rule for 'UP'/'DOWN' filtering of logFC sign
    # Many bulk-ish methods in your original script used <0 for "UP".
    up_sign = "negative"  # "positive" or "negative"
  ),
  DEVIL = list(
    folder = "devil_raw",
    pattern = "*\\.rda$",
    up_sign = "positive"  # adjust if your contrast was reversed
  ),
  NEBULA = list(
    folder = "nebula_raw",
    pattern = "*\\.rda$",
    up_sign = "positive"  # adjust if needed
  )
)

ALPHA <- 0.05

# -------------------- HELPERS (IO & NORMALIZATION) --------------------------

# Read counts (genes x cells). If first column looks like gene ids, use as rownames.
read_counts_matrix <- function(path_counts) {
  df <- readr::read_tsv(path_counts, col_types = cols())
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
  mat
}

# Load ground truth (Up/Down/All) from a simulation dir in DataPrep
load_ground_truth <- function(sim_dir) {
  true_up   <- readr::read_tsv(fs::path(sim_dir, "true_up_genes.txt"), col_names = FALSE, show_col_types = FALSE)[[1]]
  true_down <- readr::read_tsv(fs::path(sim_dir, "true_down_genes.txt"), col_names = FALSE, show_col_types = FALSE)[[1]]
  counts    <- read_counts_matrix(fs::path(sim_dir, "counts.txt"))
  all_genes <- rownames(counts)
  list(Up = true_up, Down = true_down, All = all_genes)
}

# Find and load your saved method result (.rda with object 'res')
# Normalize to data.frame with rownames = gene and columns: pval, adj.pval, log2FC
load_method_res <- function(sim_name, method) {
  spec <- METHOD_SPECS[[method]]
  mdir <- fs::path(RESULTS_ROOT, sim_name, spec$folder)
  if (!fs::dir_exists(mdir)) return(NULL)
  
  # look for "<sim>+<method>*.rda"
  files <- list.files(mdir, full.names = T)
  if (!length(files)) return(NULL)
  
  # pick the newest in case multiple
  f <- files[1]
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  #load(f)
  
  if (!exists("res", envir = env)) return(NULL)
  df <- get("res", envir = env)
  
  # normalize columns
  # accepted column names:
  # gene|Gene, pval|pvalue|PValue, padj|adjpvalue|FDR, lfc|logFC|log2FC
  nm <- tolower(colnames(df))
  colnames(df)[match(nm, nm)] <- nm
  
  # add gene column if not present
  if (!("gene" %in% nm)) {
    # if rownames are genes
    if (!is.null(rownames(df))) {
      df$gene <- rownames(df)
    } else {
      stop("No 'gene' column and no rownames found in: ", f)
    }
  }
  
  # create canonical columns
  df$pval    <- df[[if ("pval" %in% nm) "pval" else if ("pvalue" %in% nm) "pvalue" else if ("pvalue" %in% names(df)) "pvalue" else if ("p.value" %in% nm) "p.value" else "pval"]]
  if ("padj" %in% nm) {
    df$adj.pval <- df$padj
  } else if ("adjpvalue" %in% nm) {
    df$adj.pval <- df$adjpvalue
  } else if ("fdr" %in% nm) {
    df$adj.pval <- df$fdr
  } else {
    # fall back if missing
    df$adj.pval <- p.adjust(df$pval, method = "BH")
  }
  
  if ("log2fc" %in% nm) {
    df$log2FC <- df$log2fc
  } else if ("logfc" %in% nm) {
    df$log2FC <- df$logfc
  } else if ("lfc" %in% nm) {
    df$log2FC <- df$lfc
  } else {
    stop("No logFC/lfc/log2FC column in: ", f)
  }
  
  rownames(df) <- df$gene
  df[, c("pval", "adj.pval", "log2FC"), drop = FALSE]
}

# sign filter per method + 'select' (UP/DOWN)
filter_by_direction <- function(tab, method, select) {
  rule <- METHOD_SPECS[[method]]$up_sign
  if (select == "UP") {
    if (rule == "positive")  tab <- tab[tab$log2FC > 0, , drop = FALSE]
    if (rule == "negative")  tab <- tab[tab$log2FC < 0, , drop = FALSE]
  } else if (select == "DOWN") {
    if (rule == "positive")  tab <- tab[tab$log2FC < 0, , drop = FALSE]
    if (rule == "negative")  tab <- tab[tab$log2FC > 0, , drop = FALSE]
  } else {
    stop("select must be 'UP' or 'DOWN'")
  }
  tab
}

# Compute F-score (beta = 0.5 like your first script)
fscore_beta <- function(PPV, TPR, beta = 0.5) {
  b2 <- beta^2
  (1 + b2) * (PPV * TPR) / (b2 * PPV + TPR)
}

# -------------------- MAIN F-SCORE FUNCTION ---------------------------------

main_Fscore <- function(select = c("UP", "DOWN")) {
  select <- match.arg(select)
  
  simul_dirs <- fs::dir_ls(DATAPREP_ROOT, type = "directory", regexp = "simul")
  if (!length(simul_dirs)) stop("No simulation folders found under ", DATAPREP_ROOT)
  
  # collect per-simulation per-method rows
  rows <- list()
  
  sim_dir = simul_dirs[1]
  for (sim_dir in simul_dirs) {
    sim_name <- fs::path_file(sim_dir)
    GT <- load_ground_truth(sim_dir)  # Up/Down/All
    
    method = "DEVIL"
    for (method in names(METHOD_SPECS)) {
      res <- load_method_res(sim_name, method)
      if (is.null(res)) next
      
      # direction & cutoff
      res <- filter_by_direction(res, method, select)
      res_sig <- res[res$adj.pval <= ALPHA, , drop = FALSE]
      
      N  <- GT$All
      GTset <- if (select == "UP") GT$Up else GT$Down
      pred  <- rownames(res_sig)
      
      TP <- sum(GTset %in% pred)
      FP <- length(pred) - TP
      FN <- length(GTset) - TP
      TN <- length(N) - TP - FP - FN
      
      TPR <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      PPV <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      TNR <- if ((TN + FP) > 0) TN / (TN + FP) else 0
      F   <- fscore_beta(PPV, TPR, beta = 0.5)
      
      rows[[length(rows) + 1]] <- data.frame(
        simu = sim_name, method = method,
        TP = TP, FN = FN, FP = FP, TN = TN,
        precision = PPV, recall = TPR, specificity = TNR,
        Fscore = F,
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) stop("No results found in ", RESULTS_ROOT)
  
  # Write a wide table like your original (optional)
  wide_conf <- out %>%
    select(simu, method, TP, FN, FP, precision, recall, specificity, Fscore) %>%
    group_by(simu, method) %>% slice(1) %>% ungroup() %>%
    arrange(simu, method)
  
  write.table(wide_conf,
              paste0("summary_confusion_matrix_", select, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Return only the F-scores (for plotting)
  fscores <- out %>% select(simu, method, Fscore)
  return(fscores)
}

# -------------------- RUN & PLOT --------------------------------------------

Fscore_up   <- main_Fscore(select = "UP")
Fscore_down <- main_Fscore(select = "DOWN")

Fall <- rbind(
  transform(Fscore_up,   direction = "UP"),
  transform(Fscore_down, direction = "DOWN")
)

Fall$method <- factor(Fall$method, levels = names(METHOD_SPECS))

p <- ggplot(Fall, aes(x = method, y = Fscore, color = method)) +
  geom_boxplot(outlier.size = 0.8) +
  coord_flip() +
  labs(y = "F-score (β = 0.5)", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

pdf("fscore_boxplot.pdf", width = 8, height = 6)
print(p)
dev.off()
