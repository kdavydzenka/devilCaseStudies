
nebula.mult <- function(count, df){
  pred <- model.matrix(~~1+tx_cell, data=df)

  data_g = nebula::group_cell(count=count, id=df$id, pred=pred)

  if (is.pb) {
    sf = devil:::calculate_sf(count, method = "psinorm")
  } else {
    sf = rep(1, ncol(count))
  }
  # offset=Matrix::colSums(count)
  # offset = devil:::calculate_sf(count, method = "psinorm")

  s <- Sys.time()
  sid <- df$id
  fit.nebula <- nebula::nebula(
    count,
    id = df$id,
    pred = pred,
    cpc=0, 
    offset = sf,
    mincp=0,
    ncore=1
  )
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()

  fit.result <- fit.nebula$summary
  rownames(fit.result) <- fit.result$gene
  
  result <- fit.result %>%
    mutate(
      Estimate=logFC_tx_cell / log(2),
      'Std. Error'=se_tx_cell,
      't value'=logFC_tx_cell/se_tx_cell,
      'Pr(>|t|)'=p_tx_cell
    ) %>%
    select(Estimate, 'Std. Error', 't value', 'Pr(>|t|)')

  result <- dplyr::as_tibble(result)
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
