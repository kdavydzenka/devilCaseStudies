require(glmGamPoi)
require(SingleCellExperiment)

glmgp.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  s <- Sys.time()
  fit <- glmGamPoi::glm_gp(sce.pb, design=~1+tx_cell, size_factors = 1)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  test <- glmGamPoi::test_de(fit,contrast = c(0,1))

  beta <- fit$Beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

glmgp.cell.mult <- function(count, df){
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  
  if (is.pb) {
    sf = devil:::calculate_sf(count, method = "psinorm")
  } else {
    sf = FALSE
  }

  s <- Sys.time()
  fit <- glmGamPoi::glm_gp(count, design_matrix, on_disk=FALSE, size_factors=sf)
  #fit <- glmGamPoi::glm_gp(count, design_matrix, , size_factors = 1)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- glmGamPoi::test_de(fit, contrast = c(0,1))

  beta <- fit$Beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

# glmGam v2 ####
glmgp.cell.fixed <- function(count, df){
  #design_matrix <- model.matrix(~1+tx_cell, data = df)
  df$tx_cell = factor(df$tx_cell, levels = c(0,1))
  df$donor_id = factor(df$donor_id, levels = unique(df$donor_id))
  design_matrix <- model.matrix(~1+tx_cell+donor_id, data = df)
  design_matrix[,1] <- 1
  
  result <- tryCatch(
    {
      s <- Sys.time()
      fit <- glmGamPoi::glm_gp(count, design_matrix, on_disk=FALSE, size_factors=FALSE)
      delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
      e <- Sys.time()
      
      contrast = rep(0, ncol(design_matrix))
      contrast[2] = 1
      test <- glmGamPoi::test_de(fit, contrast = contrast)
      
      beta <- fit$Beta[,2]
      pval <- test$pval
      tval <- qnorm(1-pval/2) * sign(beta)
      se <- beta/tval
      
      result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
      result$delta_time <- delta_time
      colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
      result
    }
    ,
    error = function(e) {
      print(e)
      matrix(NA, nrow = dim(count)[1], ncol = 5)
    }
  )
  
  return(result %>% as.matrix())
}
