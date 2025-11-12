
library(FLASHMM)

flash.base <- function(count, df){
  X <- model.matrix(~1+tx_cell, data = df)
  Y <- log(count + 1)
  Z <- model.matrix(~ 0 + id, data = df)
  d <- ncol(Z)
  
  result <- tryCatch(
    {
      s <- Sys.time()
      fit <- lmmfit(Y, X, Z, d = d)
      e <- Sys.time()
      delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
      
      test <- lmmtest(fit, contrast=as.array(c(0,1)))
      
      beta = test[,1]
      se = beta * 0
      tval = test[,2]
      pval = test[,3]
      
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