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
                          init_overdispersion = NULL, overdispersion = F)
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

devil.overdisp <- function(count, df) {
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
                          init_overdispersion = NULL, overdispersion = T)
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
