
library(devil)

devil.pure <- function(count, df) {
  log_library_size = log(colSums(count))
  
  df$tx_cell = as.factor(df$tx_cell)
  df$lls = log_library_size
  design_matrix <- model.matrix(~lls+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=NULL, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = F)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}


devil.fit <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  if (is.pb) {
    sf = "psinorm"
  } else {
    sf = NULL
  }
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=sf, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = F, 
                          init_beta_rough = TRUE)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$input_parameters$parallel = FALSE
  
  fit$overdispersion = rep(0, length(fit$overdispersion))
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.overdisp.new <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  if (is.pb) {
    sf = "psinorm"
  } else {
    sf = NULL
  }
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=sf, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = "new", 
                          init_beta_rough = TRUE)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.overdisp.old <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  if (is.pb) {
    sf = "psinorm"
  } else {
    sf = NULL
  }
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=sf, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = "old", 
                          init_beta_rough = TRUE)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.overdisp.MOM <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  if (is.pb) {
    sf = "psinorm"
  } else {
    sf = NULL
  }
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=sf, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = "MOM", 
                          init_beta_rough = TRUE)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  # 
  # estimate_rough_dispersion = function(count_matrix, design_matrix, beta_matrix, sf) {
  #   inner_estimate_rough_dispersion = function(design_matrix, ys, beta, sf) {
  #     n = length(ys)
  #     p = length(beta)
  #     
  #     mus = as.numeric(sf * exp(design_matrix %*% beta))
  #     
  #     num = sum((ys - mus)^2 - mus)
  #     den = sum(mus^2)
  #     
  #     corr = n / (n - p)
  #     
  #     theta_hat = corr * num / den
  #     max(theta_hat, 0)  # truncate to non-negative
  #   }  
  #   
  #   
  #   ngenes = nrow(count_matrix)
  #   lapply(1:ngenes, function(gene_idx) {
  #     inner_estimate_rough_dispersion(design_matrix = design_matrix, 
  #                                     ys = count_matrix[gene_idx,], 
  #                                     beta = beta_matrix[gene_idx,], 
  #                                     sf = sf)
  #   }) %>% unlist()
  # }
  # 
  # rough_disps = estimate_rough_dispersion(fit$input_matrix, 
  #                                         design_matrix = fit$design_matrix, 
  #                                         beta_matrix = fit$beta, 
  #                                         sf = fit$size_factors)
  # 
  # # fit$overdispersion = devil:::estimate_dispersion(count, fit$offset_vector)
  # fit$overdispersion = rough_disps
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.overdisp.SC <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  if (is.pb) {
    sf = "psinorm"
  } else {
    sf = NULL
  }
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=sf, 
                          verbose=F, parallel.cores=1, 
                          init_overdispersion = NULL, overdispersion = FALSE)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  estimate_dispersion_sc <- function(count_matrix,
                                     design_matrix,
                                     beta_matrix,
                                     sf,
                                     do_newton_for_outliers = TRUE,
                                     newton_steps = 1L,
                                     log_resid_thresh = 1,     # |log θ - trend| > 1 ⇒ refine
                                     shrinkage_strength = 0.3, # 0 = no shrink, 1 = full trend
                                     trend_span = 0.5) {
    # count_matrix: G x N
    # design_matrix: N x p
    # beta_matrix: G x p
    # sf: length N
    
    G <- nrow(count_matrix)
    N <- ncol(count_matrix)
    p <- ncol(design_matrix)
    
    if (nrow(beta_matrix) != G) {
      stop("beta_matrix must have one row per gene.")
    }
    if (nrow(design_matrix) != N) {
      stop("design_matrix must have one row per cell (columns of count_matrix).")
    }
    if (length(sf) != N) {
      stop("sf must have length equal to number of cells (columns of count_matrix).")
    }
    
    # -----------------------------
    # 1) Compute μ for all genes
    # -----------------------------
    # eta_mat: N x G (cells x genes)
    eta_mat <- design_matrix %*% t(beta_matrix)        # big BLAS call
    mu_mat  <- sweep(exp(eta_mat), 1, sf, `*`)         # multiply each row by sf
    # Reorient to G x N to match count_matrix
    mu_mat  <- t(mu_mat)                               # genes x cells
    
    # -----------------------------
    # 2) Rough GLM-aware MoM θ0
    # -----------------------------
    Y   <- count_matrix
    mu  <- mu_mat
    rm(eta_mat)  # free some memory
    
    num <- rowSums((Y - mu)^2 - mu)
    den <- rowSums(mu^2)
    
    df_corr <- N / (N - p)   # df correction
    theta0  <- df_corr * num / den
    
    # truncate negatives
    theta0[theta0 < 0 | !is.finite(theta0)] <- 0
    
    # -----------------------------
    # 3) Mean–dispersion trend (initial, using θ0)
    # -----------------------------
    eps <- 1e-8
    mu_bar <- rowMeans(mu)   # mean μ per gene
    
    valid <- (theta0 > 0) & is.finite(theta0) & (mu_bar > 0)
    log_mu  <- log(mu_bar[valid])
    log_th0 <- log(theta0[valid])
    
    if (sum(valid) < 10) {
      warning("Too few valid genes for trend; skipping shrinkage & Newton logic.")
      theta_refined <- theta0
      theta_trend   <- rep(NA_real_, G)
      theta_shrunk  <- theta_refined
      return(list(
        theta_raw    = theta_refined,
        theta_trend  = theta_trend,
        theta_shrunk = theta_shrunk
      ))
    }
    
    trend_fit0 <- loess(log_th0 ~ log_mu, span = trend_span)
    trend_pred0 <- rep(NA_real_, G)
    trend_pred0[valid] <- predict(trend_fit0,
                                  newdata = data.frame(log_mu = log(mu_bar[valid])))
    
    # -----------------------------
    # 4) Newton refinement for outliers (optional)
    # -----------------------------
    theta_refined <- theta0
    
    if (isTRUE(do_newton_for_outliers)) {
      # define "problematic" genes: big deviation from trend or zero θ
      log_theta0_full <- log(theta0 + eps)
      dev <- abs(log_theta0_full - trend_pred0)
      dev[!is.finite(dev)] <- Inf
      
      is_outlier <- (theta0 == 0) | (!is.finite(theta0)) |
        (dev > log_resid_thresh)
      
      out_idx <- which(is_outlier)
      
      if (length(out_idx) > 0L) {
        for (g in out_idx) {
          ys  <- Y[g, ]
          mu_g <- mu[g, ]
          th  <- theta0[g]
          if (!is.finite(th) || th < 0) th <- 0
          
          # Perform 1 or more Newton steps on θ
          for (step in seq_len(newton_steps)) {
            # Score U(θ) and info I(θ) for NB2-style update
            denom <- (1 + th * mu_g)^2
            U <- sum(((ys - mu_g)^2 - mu_g) / denom)
            I <- sum(mu_g^2 / denom)
            
            if (!is.finite(U) || !is.finite(I) || I <= 0) break
            
            th_new <- th + U / I
            if (!is.finite(th_new) || th_new < 0) {
              th_new <- 0
              break
            }
            th <- th_new
          }
          theta_refined[g] <- th
        }
      }
      
      # non-outliers keep θ0
    }
    
    # -----------------------------
    # 5) Refit trend on refined θ
    # -----------------------------
    valid2 <- (theta_refined > 0) & is.finite(theta_refined) & (mu_bar > 0)
    log_th_ref <- log(theta_refined[valid2])
    log_mu2    <- log(mu_bar[valid2])
    
    trend_fit <- loess(log_th_ref ~ log_mu2, span = trend_span)
    trend_pred <- rep(NA_real_, G)
    trend_pred[valid2] <- predict(trend_fit,
                                  newdata = data.frame(log_mu2 = log(mu_bar[valid2])))
    
    theta_trend <- exp(trend_pred)
    
    # -----------------------------
    # 6) Shrinkage: combine θ_refined with trend
    # -----------------------------
    ss <- shrinkage_strength
    ss <- max(min(ss, 1), 0)  # clamp to [0,1]
    
    log_theta_ref   <- log(theta_refined + eps)
    log_theta_trend <- log(theta_trend + eps)
    
    log_theta_shrunk <- (1 - ss) * log_theta_ref + ss * log_theta_trend
    # Where trend is NA, fall back to refined
    log_theta_shrunk[is.na(log_theta_trend)] <- log_theta_ref[is.na(log_theta_trend)]
    
    theta_shrunk <- exp(log_theta_shrunk)
    
    list(
      theta_raw    = theta_refined,  # per-gene, MoM + possible Newton
      theta_trend  = theta_trend,    # smooth mean–dispersion curve
      theta_shrunk = theta_shrunk    # shrunken per-gene dispersions
    )
  }
  
  rough_disps = estimate_dispersion_sc(count_matrix = fit$input_matrix, 
                                    design_matrix = fit$design_matrix, 
                                    beta_matrix = fit$beta, 
                                    sf = fit$size_factors, 
                                    do_newton_for_outliers = TRUE, 
                                    newton_steps = 1, shrinkage_strength = .3, trend_span = 1)
  
  # fit$overdispersion = devil:::estimate_dispersion(count, fit$offset_vector)
  fit$overdispersion = rough_disps$theta_shrunk
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters, parallel.cores = 1)
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.base <- function(count, df){
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  
  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, size_factors=NULL, verbose=T, parallel.cores=1, init_overdispersion = 100)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$beta %>% sum()
  fit$overdispersion %>% sum()
  
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)))
  
  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.mixed <- function(count, df) {
  df$tx_cell = as.factor(df$tx_cell)
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(paste0(df$id))
  
  s <- Sys.time()
  #fit <- devil::fit_devil(count, design_matrix, size_factors="psinorm", verbose=F, parallel.cores=1, init_overdispersion = 100)
  fit <- devil::fit_devil(count, design_matrix, size_factors=NULL, verbose=F, parallel.cores=1, init_overdispersion = 100)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit$input_parameters$parallel = FALSE
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters)

  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
