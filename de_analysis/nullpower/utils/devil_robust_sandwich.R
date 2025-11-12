
# --- Compile the C++ helpers (adjust path if needed) --------------------------
Rcpp::sourceCpp("utils/robust_sandwich.cpp")

# --- Helpers ------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# numerically stable solve(A) %*% M %*% solve(A)
# tries Cholesky; if that fails, adds a small ridge; if still fails, uses ginv
solve_A_M_A <- function(A, M, ridge0 = 1e-8, max_tries = 5) {
  ok <- FALSE
  r  <- ridge0
  for (k in 0:max_tries) {
    A_try <- if (k == 0) A else A + diag(r, nrow(A))
    cholA <- tryCatch(chol(A_try), error = function(e) NULL)
    if (!is.null(cholA)) {
      Ainv <- chol2inv(cholA)
      return(Ainv %*% M %*% Ainv)
    }
    r <- r * 10
  }
  # last resort
  Ainv <- MASS::ginv(A)
  Ainv %*% M %*% Ainv
}

# working weights for NB2 (matches your C++)
nb2_working <- function(y, X, beta, size_factors, overdispersion) {
  eta <- as.vector(X %*% beta)
  mu  <- size_factors * exp(eta)
  w   <- mu / (1 + mu / overdispersion)
  list(mu = mu, w = w)
}

# Effective df for one contrast using cluster influence (handles unbalanced clusters).
# df_eff ≈ G_eff - 1 with G_eff = (sum a_g)^2 / sum a_g^2 where
# a_g = c' A^{-1} (X_g' W_g X_g) A^{-1} c
effective_df_contrast <- function(X, contrast, clusters, w) {
  contrast <- as.numeric(contrast)
  Whalf <- sqrt(pmax(w, 0))
  Xw    <- X * Whalf
  # A = X' W X
  A <- crossprod(Xw, Xw)
  
  # Robust inversion of A
  Ainv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  
  split_idx <- split(seq_len(nrow(X)), as.integer(as.factor(clusters)))
  a_g <- vapply(split_idx, function(idx) {
    Xg <- X[idx, , drop = FALSE]
    wg <- pmax(w[idx], 0)
    Xgw <- Xg * sqrt(wg)
    Sg <- crossprod(Xgw, Xgw)  # X_g' W_g X_g
    as.numeric(t(contrast) %*% Ainv %*% Sg %*% Ainv %*% contrast)
  }, numeric(1))
  
  sum_a  <- sum(a_g)
  sum_a2 <- sum(a_g^2)
  
  if (!is.finite(sum_a) || !is.finite(sum_a2) || sum_a2 <= .Machine$double.eps) return(1)
  
  G_eff <- (sum_a^2) / sum_a2
  df_eff <- max(G_eff - 1, 1)
  df_eff
}

# --- 1) CR2-like sandwich VCOV for one gene -----------------------------------
# V = A^{-1} * (cfac * B_CR) * A^{-1}
# A = expected Hessian (X'WX), B_CR = clustered meat from your C++.
# cfac = small-sample HC1-style finite-sample factor based on G, N, K.
compute_sandwich <- function(design_matrix, y, beta, overdispersion, size_factors, clusters) {
  A <- compute_hessian_expected(beta, overdispersion, y, design_matrix, size_factors)  # X' W X
  B <- compute_clustered_meat(design_matrix, y, beta, overdispersion, size_factors, clusters)
  
  N <- nrow(design_matrix); K <- ncol(design_matrix)
  G <- dplyr::n_distinct(clusters)
  
  # Finite-sample correction (HC1-style). This is conservative; df will further adjust via Satterthwaite.
  cfac <- if (G > 1 && N > K + 1) (G/(G-1)) * ((N-1)/(N-K)) else 1
  
  # Stable A^{-1} M A^{-1}
  solve_A_M_A(A, cfac * B)
}

# --- 2) Differential test (one contrast) with CR2 + Satterthwaite df ----------
test_de_v2 <- function(devil.fit, contrast, pval_adjust_method = "BH", max_lfc = 10, clusters = NULL, parallel.cores=NULL) {
  
  # Detect cores to use
  max.cores <- parallel::detectCores()
  if (is.null(parallel.cores)) {
    n.cores = max.cores
  } else {
    if (parallel.cores > max.cores) {
      message(paste0("Requested ", parallel.cores, " cores, but only ", max.cores, " available."))
    }
    n.cores = min(max.cores, parallel.cores)
  }
  
  # Extract necessary information
  ngenes <- nrow(devil.fit$input_matrix)
  nsamples <- nrow(devil.fit$design_matrix)
  contrast <- as.array(contrast)
  
  # Calculate log fold changes
  lfcs <- (devil.fit$beta %*% contrast) %>%
    unlist() %>%
    unname() %>%
    c()
  
  # Calculate p-values in parallel
  if (!is.null(clusters)) {
    
    if (!is.numeric(clusters)) {
      message("Converting clusters to numeric factors")
      clusters = as.numeric(as.factor(clusters))
    }
    
    p_values <- parallel::mclapply(1:nrow(devil.fit$input_matrix), function(gene_idx) {
      mu_test <- lfcs[gene_idx]
      
      H <- compute_sandwich(
        devil.fit$design_matrix,
        devil.fit$input_matrix[gene_idx,],
        devil.fit$beta[gene_idx,], devil.fit$overdispersion[gene_idx],
        devil.fit$size_factors,
        clusters
      )

      total_variance <- t(contrast) %*% H %*% contrast
      p_val_clust = 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = F)
      
      H <- devil:::compute_hessian(devil.fit$beta[gene_idx,], 1 / devil.fit$overdispersion[gene_idx], devil.fit$input_matrix[gene_idx,], devil.fit$design_matrix, devil.fit$size_factors)
      total_variance <- t(contrast) %*% H %*% contrast
      p_val_indep = 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = F)
      
      max(p_val_clust, p_val_indep)
    }, mc.cores = n.cores) %>% unlist()
    
  } else {
    p_values <- parallel::mclapply(1:ngenes, function(gene_idx) {
      mu_test <- lfcs[gene_idx]
      H <- compute_hessian(devil.fit$beta[gene_idx,], 1 / devil.fit$overdispersion[gene_idx], devil.fit$input_matrix[gene_idx,], devil.fit$design_matrix, devil.fit$size_factors)
      total_variance <- t(contrast) %*% H %*% contrast
      #1 - stats::pchisq(mu_test^2 / total_variance, df = 1)
      2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = F)
    }, mc.cores = n.cores) %>% unlist()
  }
  
  # Create tibble with results
  result_df <- dplyr::tibble(
    name = rownames(devil.fit$beta),
    pval = p_values,
    adj_pval = stats::p.adjust(p_values, method = pval_adjust_method),
    lfc = lfcs / log(2)
  )
  
  # Filter results based on max_lfc
  result_df <- result_df %>%
    dplyr::mutate(lfc = ifelse(.data$lfc >= max_lfc, max_lfc, .data$lfc)) %>%
    dplyr::mutate(lfc = ifelse(.data$lfc <= -max_lfc, -max_lfc, .data$lfc))
  
  # if (sum(is.na(result_df))) {
  #   message('Warning: the results for some genes are unrealiable (i.e. NaN)\n This might be due to gene very lowly expressed or not expressed at all for some conditions')
  # }
  # if (sum(is.na(result_df))) {
  #   na_genes_idxs <- which(is.na(result_df$pval))
  #   dm <- as.matrix(devil.fit$design_matrix[,contrast != 0])
  #   cell_idx <- dm != 0
  #   dm <- as.matrix(dm[cell_idx,])
  #   tmp <- lapply(na_genes_idxs, function(gene_idx) {
  #     beta0 <- init_beta(t(devil.fit$input_matrix[gene_idx,cell_idx]), design_matrix = dm, offset_matrix = devil.fit$offset_matrix[gene_idx,cell_idx])
  #     new_beta <- beta_fit(devil.fit$input_matrix[gene_idx,cell_idx], X = dm, mu_beta = beta0, off = devil.fit$offset_matrix[gene_idx,cell_idx], k = 1 / devil.fit$overdispersion[gene_idx], max_iter = 500, eps = 1e-3)
  #     new_beta <- new_beta$mu_beta
  #     mu_test <- sum(new_beta %*% contrast)
  #     if (!is.null(clusters)) {
  #
  #       H <- compute_sandwich(
  #         devil.fit$design_matrix,
  #         devil.fit$input_matrix[gene_idx,],
  #         new_beta, devil.fit$overdispersion[gene_idx],
  #         devil.fit$size_factors,
  #         clusters
  #       )
  #
  #       total_variance <- t(contrast) %*% H %*% contrast
  #       new_pval <- 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = F)
  #     } else {
  #       H <- compute_hessian(new_beta, 1 / devil.fit$overdispersion[gene_idx], devil.fit$input_matrix[gene_idx,], devil.fit$design_matrix, devil.fit$size_factors)
  #       total_variance <- t(contrast) %*% H %*% contrast
  #       new_pval <- 2 * stats::pt(abs(mu_test) / sqrt(total_variance), df = nsamples - 2, lower.tail = F)
  #     }
  #     result_df$pval[gene_idx] <<- new_pval
  #   })
  #
  #   #message('Warning: the results for some genes are unrealiable (i.e. NaN)\n This might be due to gene very lowly expressed or not expressed at all for some conditions')
  # }
  
  result_df$adj_pval = stats::p.adjust(result_df$pval, method = pval_adjust_method)
  
  return(result_df)
}


new_test_de <- function(devil.fit,
                        contrast,
                        pval_adjust_method = "BH",
                        max_lfc = 10,
                        clusters,
                        parallel.cores = NULL) {
  
  stopifnot(is.matrix(devil.fit$design_matrix))
  stopifnot(is.matrix(devil.fit$beta))
  stopifnot(is.matrix(devil.fit$input_matrix))
  
  X <- devil.fit$design_matrix
  p <- ncol(X)
  
  contrast <- as.numeric(contrast)
  if (length(contrast) != p) stop("contrast must have length = ncol(X)")
  
  cvec <- matrix(contrast, ncol = 1)
  lfcs <- drop(devil.fit$beta %*% cvec)  # linear contrast on beta (per gene)
  
  max.cores <- parallel::detectCores()
  n.cores <- if (is.null(parallel.cores)) max.cores else min(max.cores, parallel.cores)
  
  cl_int <- as.integer(as.factor(clusters))
  ngenes <- nrow(devil.fit$input_matrix)
  
  g = which(rownames(devil.fit$input_matrix) == "ENSG00000106070")
  
  p_values <- parallel::mclapply(
    X = seq_len(ngenes),
    mc.cores = n.cores,
    FUN = function(g) {
      y_g    <- devil.fit$input_matrix[g, ]
      beta_g <- devil.fit$beta[g, ]
      phi_g  <- devil.fit$overdispersion[g]
      sf_g   <- devil.fit$size_factors %||% rep(1, length(y_g))
      
      # working weights for df_eff (depend on beta_g, phi_g)
      ww <- nb2_working(y = y_g, X = X, beta = beta_g,
                        size_factors = sf_g, overdispersion = phi_g)$w
      
      # CR2-like VCOV for the whole beta vector
      V <- devil:::compute_sandwich(design_matrix = X, y = y_g, beta = beta_g,
                                    overdispersion = 1 / phi_g, size_factors = sf_g,
                                    clusters = cl_int)
      
      H = devil:::compute_hessian(beta = beta_g, overdispersion = 1 / phi_g, y = y_g, design_matrix = X, size_factors = sf_g)
      
      # contrast variance
      v_c <- as.numeric(t(contrast) %*% V %*% contrast)
      if (!is.finite(v_c) || v_c <= .Machine$double.eps) {
        return(1.0)  # ultra-conservative fallback
      }
      
      # Satterthwaite/Bell–McCaffrey effective df for this gene/contrast
      df_eff <- effective_df_contrast(X, contrast, cl_int, ww)
      df_eff = max(cl_int) - 1
      
      # Two-sided t p-value with df_eff
      mu_hat <- lfcs[g]
      tstat  <- mu_hat / sqrt(v_c)
      2 * stats::pt(abs(tstat), df = df_eff, lower.tail = FALSE)
    }
  ) %>% unlist(use.names = FALSE)
  
  tibble::tibble(
    name     = rownames(devil.fit$beta),
    pval     = p_values,
    adj_pval = stats::p.adjust(p_values, method = pval_adjust_method),
    lfc      = (lfcs / log(2)) |> pmin(max_lfc) |> pmax(-max_lfc)
  )
}

# --- 3) Convenience wrapper to fit and report per-gene coefficients ------------
new_devil <- function(count, df) {
  df <- tibble::as_tibble(df)
  df$tx_cell <- as.factor(df$tx_cell)
  X <- model.matrix(~ 1 + tx_cell, data = df)
  clusters <- df$id
  
  # Fit DEVIL
  s <- Sys.time()
  fit <- devil::fit_devil(
    count, X,
    size_factors = "psinorm",
    verbose = FALSE,
    parallel.cores = 1,
    init_overdispersion = 100
  )
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  # Test the coefficient on tx_cell
  p <- ncol(X)
  contrast <- rep(0, p); contrast[2] <- 1
  
  test <- test_de_v2(
    devil.fit = fit,
    contrast = contrast,
    pval_adjust_method = "BH",
    max_lfc = 100,
    clusters = clusters,
    parallel.cores = 1
  )
  
  beta <- fit$beta[, 2]
  pval <- test$pval
  
  # Back out t-stats via normal approx is unstable with small df; instead recompute from p + sign
  # Here we use normal only to get a signed "t-like" value for the summary table.
  tval <- stats::qnorm(pmax(1e-300, 1 - pval / 2)) * sign(beta)
  se   <- beta / pmax(abs(tval), .Machine$double.eps)
  
  res <- tibble::tibble(
    `Estimate`   = beta,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = pval,
    `Time`       = delta_time
  )
  as.matrix(res)
}

new_devil_2 <- function(count, df) {
  df <- tibble::as_tibble(df)
  df$tx_cell <- as.factor(df$tx_cell)
  X <- model.matrix(~ 1 + tx_cell, data = df)
  clusters <- df$id
  
  # Fit DEVIL
  s <- Sys.time()
  fit <- devil::fit_devil(
    count, X,
    size_factors = "psinorm",
    verbose = FALSE,
    parallel.cores = 1,
    init_overdispersion = 100
  )
  e <- Sys.time()
  delta_time <- as.numeric(difftime(e, s, units = "secs"))
  
  # Test the coefficient on tx_cell
  p <- ncol(X)
  contrast <- rep(0, p); contrast[2] <- 1
  
  test <- test_de_v2(
    devil.fit = fit,
    contrast = contrast,
    pval_adjust_method = "BH",
    max_lfc = 100,
    clusters = clusters,
    parallel.cores = 1
  )
  
  beta <- fit$beta[, 2]
  pval <- test$pval
  
  # Back out t-stats via normal approx is unstable with small df; instead recompute from p + sign
  # Here we use normal only to get a signed "t-like" value for the summary table.
  tval <- stats::qnorm(pmax(1e-300, 1 - pval / 2)) * sign(beta)
  se   <- beta / pmax(abs(tval), .Machine$double.eps)
  
  res <- tibble::tibble(
    `Estimate`   = beta,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = pval,
    `Time`       = delta_time
  )
  as.matrix(res)
}
