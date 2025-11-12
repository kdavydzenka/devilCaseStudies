require(Rcpp)
require(Matrix)
require(data.table)
require(SingleCellExperiment)
require(dplyr)
require(tidyr)
require(magrittr)
require(Rcpp)

selectCol <- function(mat, j.col){
  x.col.dense <- rep(0,nrow(mat))
  p.begin <- mat@p[j.col]+1
  p.end <- mat@p[j.col+1]
  i.col <- mat@i[p.begin:p.end]+1 # i counts from 0
  x.col <- mat@x[p.begin:p.end]
  x.col.dense[i.col] <- x.col
  return(x.col.dense)
}

selectCols <- function(mat, j.cols){
  return(sapply(j.cols, selectCol, mat=mat))
}


src <-
  "
#include <Rcpp.h>

// [[Rcpp::export]]
void vec_down_sample(
    Rcpp::NumericVector data,
    const Rcpp::LogicalVector which,
    int begin,
    int end,
    double prob
    ){
    for(int i=begin; i<end; i++){
        if(which[i]){
            data[i] = data[i] + R::rbinom(data[i], prob);
          }
        }
    }
"
sourceCpp(code = src)

# # Pure R version of vec_down_sample (no Rcpp needed)
# vec_down_sample <- function(data, sel, begin, end, prob, end_inclusive = FALSE) {
#   stopifnot(is.numeric(data), is.logical(sel), length(sel) >= length(data))
#   stopifnot(is.numeric(begin), is.numeric(end), length(begin) == 1, length(end) == 1)
#   stopifnot(is.numeric(prob), length(prob) == 1, prob >= 0, prob <= 1)
#   
#   n <- length(data)
#   # Emulate the C++ for (i = begin; i < end; i++)   (end is exclusive by default)
#   i1 <- max(1L, as.integer(begin))
#   i2 <- as.integer(end) - if (end_inclusive) 0L else 1L
#   i2 <- min(n, i2)
#   if (i1 > i2) return(data)
#   
#   idx <- i1:i2
#   # pick positions where sel is TRUE in the range
#   pos <- idx[sel[idx]]
#   if (!length(pos)) return(data)
#   
#   # Binomial sizes must be non-negative integers
#   sizes <- pmax(0L, as.integer(floor(data[pos])))
#   inc   <- stats::rbinom(length(pos), size = sizes, prob = prob)
#   
#   data[pos] <- data[pos] + inc
#   data
# }

downCells <- function(spmat, i.rows, j.cols, p){
  data <- spmat@x
  i.bool <- spmat@i %in% (i.rows-1) # spmat@i begins from 0, i.rows begins from 1
  for (j.col in j.cols){
    begin <- spmat@p[j.col]
    end <- spmat@p[j.col+1]
    vec_down_sample(data, i.bool, begin, end, p)
  }
}


qc_filter_cells <- function(cnt, coldata, min_umi = 500, min_genes = 200) {
  stopifnot(ncol(cnt) == nrow(coldata))
  libsz   <- colSums(cnt)
  detgn   <- Matrix::colSums(cnt > 0)
  keep    <- (libsz >= min_umi) & (detgn >= min_genes)
  list(
    cnt = cnt[, keep, drop = FALSE],
    coldata = coldata[keep, , drop = FALSE]
  )
}

# ------- Gene filter: CPM-based, then cap to MAX_GENE -------------------------
select_genes_cpm <- function(cnt, max_gene = 5000, cpm_thresh = 1) {
  libsz <- colSums(cnt); libsz[libsz == 0] <- 1
  cpm   <- t(t(cnt) / libsz) * 1e6
  keep  <- which(rowMeans(cpm) > cpm_thresh)
  if (length(keep) == 0) stop("No genes pass CPM threshold.")
  sel   <- if (length(keep) > max_gene) sample(keep, max_gene) else keep
  cnt[sel, , drop = FALSE]
}

# ------- Equalize cells per donor (after QC) ----------------------------------
equalize_cells_per_donor <- function(coldata, target = c("min","fixed"), fixed_n = NULL) {
  target <- match.arg(target)
  cd <- coldata %>% mutate(.row = row_number())
  # cells per donor
  tab <- cd %>% count(donor_id, name = "n")
  n_target <- if (target == "min") min(tab$n) else fixed_n
  if (is.null(n_target) || n_target < 1) stop("Invalid target per-donor cell count.")
  # sample per donor
  idx_keep <- cd %>%
    group_by(donor_id) %>%
    reframe(.row = sample(.row, size = min(n_target, n()), replace = n() < n_target)) %>%
    pull(.row)
  idx_keep
}

# ------- Assign treatment at donor level (is.pb = TRUE) -----------------------
assign_treatment_subject <- function(coldata) {
  donors <- unique(coldata$donor_id)
  if (length(donors) < 2) stop("Need at least 2 donors for subject-level assignment.")
  n_tx <- floor(length(donors) / 2)
  # ensure both groups non-empty
  tx_donors <- sample(donors, n_tx)
  tx_map <- setNames(integer(length(donors)), donors)
  tx_map[tx_donors] <- 1L
  as.integer(tx_map[coldata$donor_id])
}

# ------- Assign treatment within donor (is.pb = FALSE) ------------------------
assign_treatment_within_donor <- function(coldata, p = 0.5) {
  # Ensure each donor has both labels; flip the smallest group if needed
  cd <- coldata %>% mutate(tx = NA_integer_)
  cd <- cd %>%
    group_by(donor_id) %>%
    group_modify(~{
      n <- nrow(.x)
      # at least 2 cells per donor are required to split
      if (n < 2) stop("Each donor needs >= 2 cells for within-donor assignment.")
      tx <- rbinom(n, 1, p)
      if (length(unique(tx)) == 1L) {
        # force one opposite label
        flip_idx <- sample.int(n, 1)
        tx[flip_idx] <- 1L - tx[flip_idx]
      }
      .x$tx <- tx
      .x
    }) %>% ungroup()
  cd$tx
}

# ------- Sanity checks --------------------------------------------------------
check_subject_case <- function(coldata) {
  # donors per arm
  donors_by_arm <- coldata %>% distinct(donor_id, tx_cell) %>% count(tx_cell, name = "n_donors")
  if (!all(c(0,1) %in% donors_by_arm$tx_cell)) stop("Both donor arms must exist.")
  # per-donor purity of label (no mixing)
  mix <- coldata %>% group_by(donor_id) %>% summarise(nlab = n_distinct(tx_cell), .groups="drop")
  if (any(mix$nlab != 1L)) stop("Subject-level case has mixed labels within a donor.")
}

check_cell_case <- function(coldata) {
  # each donor has both labels
  ok <- coldata %>% group_by(donor_id) %>% summarise(nlab = n_distinct(tx_cell), .groups="drop")
  if (any(ok$nlab < 2L)) stop("Cell-level case requires both labels within each donor.")
}
