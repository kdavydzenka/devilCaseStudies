
# ---------- run methods (illustrative baselines) ----------
run_edgeR_with_batch <- function(Y, meta){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=Y), colData=meta)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(patient, condition, batch),
    verbose=FALSE
  )
  
  design <- model.matrix(~batch+condition, data=colData(sce.pb))
  s <- Sys.time()
  # edger.obj <- edgeR::DGEList(counts(sce.pb))
  edger.obj <- edgeR::DGEList(as.matrix(round(counts(sce.pb))))
  edger.obj <- edgeR::estimateDisp(edger.obj, design)
  fit <- edgeR::glmQLFit(y=edger.obj, design=design)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- edgeR::glmTreat(fit, coef=ncol(design))
  
  beta <- test$coefficients[,ncol(test$coefficients)]
  pval <- test$table[,'PValue']
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  tibble(gene = rownames(Y),
         pval = pval,
         p_adj = p.adjust(pval, "BH"),
         lfc = beta,
         method = "edgeR+batch")
}

run_edgeR_wo_batch <- function(Y, meta){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=Y), colData=meta)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(patient, condition),
    verbose=FALSE
  )
  
  design <- model.matrix(~condition, data=colData(sce.pb))
  s <- Sys.time()
  edger.obj <- edgeR::DGEList(as.matrix(round(counts(sce.pb))))
  # edger.obj <- edgeR::DGEList(counts(sce.pb))
  edger.obj <- edgeR::estimateDisp(edger.obj, design)
  fit <- edgeR::glmQLFit(y=edger.obj, design=design)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- edgeR::glmTreat(fit, coef=2)
  
  beta <- test$coefficients[,2]
  pval <- test$table[,'PValue']
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval
  
  tibble(gene = rownames(Y),
         pval = pval,
         p_adj = p.adjust(pval, "BH"),
         lfc = beta,
         method = "edgeR+nobatch")
}

run_nebula_w_batch <- function(Y, meta){
  pred <- model.matrix(~ batch + condition, data = meta)
  s <- Sys.time()
  sid <- meta$patient
  fit.nebula <- nebula::nebula(
    Y,
    id = meta$patient,
    pred = pred,
    cpc=0,
    mincp=0,
    ncore=1
  )
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit.result <- fit.nebula$summary
  rownames(fit.result) <- fit.result$gene
  
  tibble(gene = rownames(Y),
         pval = fit.result$p_conditionB,
         p_adj = p.adjust(fit.result$p_conditionB, "BH"),
         lfc = fit.result$logFC_conditionB,
         method = "NEBULA+batch")
}

run_nebula_wo_batch <- function(Y, meta){
  pred <- model.matrix(~ condition, data = meta)
  s <- Sys.time()
  sid <- meta$patient
  fit.nebula <- nebula::nebula(
    Y,
    id = meta$patient,
    pred = pred,
    cpc=0,
    mincp=0,
    ncore=1
  )
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  
  fit.result <- fit.nebula$summary
  rownames(fit.result) <- fit.result$gene
  
  tibble(gene = rownames(Y),
         pval = fit.result$p_conditionB,
         p_adj = p.adjust(fit.result$p_conditionB, "BH"),
         lfc = fit.result$logFC_conditionB,
         method = "NEBULA+nobatch")
}

# DEVIL stub: plug your function here. Expect a per-cluster DE test with clustered SE by patient.

run_DEVIL_batch_and_sandwich <- function(Y, meta, cluster_id){
  Y = as.matrix(Y)
  
  design_matrix = model.matrix(~ batch + condition, data = meta)
  
  fit = devil::fit_devil(input_matrix = Y, design_matrix = design_matrix, overdispersion = T, 
                         init_overdispersion = 10, size_factors = "psinorm", parallel.cores = 1, verbose = TRUE)
  
  contrast = rep(0, ncol(design_matrix))
  contrast[length(contrast)] = 1
  res = devil::test_de(devil.fit = fit, contrast = contrast, clusters = meta$patient, parallel.cores = 1)
  
  tibble(gene = paste0("Gene", 1:nrow(res)),
         pval = res$pval,
         p_adj = p.adjust(res$pval, "BH"),
         lfc = res$lfc,
         method = "devil+batch+se")
}




run_DEVIL_batch_and_sandwich <- function(Y, meta, cluster_id){
  
  Y = as.matrix(Y)
  
  # With or withoud batch
  design_matrix = model.matrix(~ batch + condition, data = meta)
  design_matrix = model.matrix(~ condition, data = meta)
  
  # "normed_sum", "psinorm", "edgeR", NULL
  fit = devil::fit_devil(input_matrix = Y, design_matrix = design_matrix, overdispersion = T, 
                         init_overdispersion = 10, size_factors = "psinorm", parallel.cores = 1, verbose = TRUE)
  
  contrast = rep(0, ncol(design_matrix))
  contrast[length(contrast)] = 1
  
  # With clusters 
  # NULL, meta$patient, meta$batch
  res = devil::test_de(devil.fit = fit, contrast = contrast, clusters = meta$patient, parallel.cores = 1)
  
  tibble(gene = paste0("Gene", 1:nrow(res)),
         pval = res$pval,
         p_adj = p.adjust(res$pval, "BH"),
         lfc = res$lfc,
         method = "devil+batch+se")
}


run_devil_grid <- function(
    Y,
    meta,
    include_batch      = c(TRUE, FALSE),
    size_factor_methods = c("normed_sum", "psinorm", "edgeR", NA),  # NA -> NULL
    cluster_by          = c("none", "patient", "batch"),
    init_overdispersion = NULL,
    overdispersion      = "new",
    parallel.cores      = 1,
    verbose_fit         = FALSE
){
  stopifnot(is.matrix(Y) || inherits(Y, "dgCMatrix"))
  Y <- as.matrix(Y)
  
  # Basic checks / coerce
  stopifnot(all(c("condition","batch","patient") %in% colnames(meta)))
  meta <- within(meta, {
    condition <- factor(condition)
    batch     <- factor(batch)
    patient   <- factor(patient)
  })
  if (length(levels(meta$condition)) != 2) {
    stop("meta$condition must have exactly 2 levels.")
  }
  
  # helper: build design safely (drop aliased columns if any)
  make_design <- function(use_batch) {
    if (isTRUE(use_batch)) {
      DM <- model.matrix(~ batch + condition, data = meta)
    } else {
      DM <- model.matrix(~ condition, data = meta)
    }
    # Drop columns with near-zero variance / alias (rare, but safe)
    qrD <- qr(DM)
    DM[, qrD$pivot[seq_len(qrD$rank)], drop = FALSE]
  }
  
  # helper: clusters vector
  get_clusters <- function(what){
    switch(what,
           "none"    = NULL,
           "patient" = meta$patient,
           "batch"   = meta$batch,
           stop("unknown cluster_by: ", what)
    )
  }
  
  # grid of settings
  grid <- expand.grid(
    include_batch      = include_batch,
    size_factors       = size_factor_methods,
    cluster_by         = cluster_by,
    stringsAsFactors   = FALSE
  )
  
  # iterate
  out_list <- lapply(seq_len(nrow(grid)), function(i){
    gb   <- grid[i, ]
    DM   <- make_design(gb$include_batch)
    # contrast on the *last* column (assumes it’s the condition_B coefficient)
    contrast <- rep(0, ncol(DM)); contrast[length(contrast)] <- 1
    
    # map size_factors argument
    sf_arg <- if (is.na(gb$size_factors)) NULL else gb$size_factors
    clus   <- get_clusters(gb$cluster_by)
    
    # run fit + test with safety net
    res_i <- tryCatch({
      fit <- devil::fit_devil(
        input_matrix       = Y,
        design_matrix      = DM,
        overdispersion     = overdispersion,
        init_overdispersion= init_overdispersion,
        size_factors       = sf_arg,
        parallel.cores     = parallel.cores,
        verbose            = verbose_fit
      )
      
      tst <- devil::test_de(
        devil.fit      = fit,
        contrast       = contrast,
        clusters       = clus,
        parallel.cores = parallel.cores
      )
      
      genes <- rownames(Y)
      if (is.null(genes)) genes <- paste0("Gene", seq_len(nrow(tst)))
      
      tibble::tibble(
        gene   = genes,
        pval   = tst$pval,
        p_adj  = p.adjust(tst$pval, method = "BH"),
        lfc    = tst$lfc,
        include_batch    = gb$include_batch,
        size_factors     = if (is.null(sf_arg)) "none" else as.character(sf_arg),
        cluster_by       = gb$cluster_by,
        method           = sprintf("devil%s+%s+%s",
                                   ifelse(gb$include_batch,"+batch",""),
                                   if (is.null(sf_arg)) "sf:none" else paste0("sf:", sf_arg),
                                   paste0("se:", gb$cluster_by))
      )
    }, error = function(e){
      # On error, return a one-row tibble tagging the failure; keeps pipeline alive
      tibble::tibble(
        gene   = NA_character_,
        pval   = NA_real_,
        p_adj  = NA_real_,
        lfc    = NA_real_,
        include_batch    = gb$include_batch,
        size_factors     = if (is.null(sf_arg)) "none" else as.character(sf_arg),
        cluster_by       = gb$cluster_by,
        method           = sprintf("ERROR: %s", conditionMessage(e))
      )
    })
    
    res_i
  })
  
  dplyr::bind_rows(out_list)
}








run_DEVIL_nobatch_sandwichpatients <- function(Y, meta, cluster_id){
  
  Y = as.matrix(Y)
  design_matrix = model.matrix(~ condition, meta)
  fit = devil::fit_devil(input_matrix = Y, design_matrix = design_matrix, overdispersion = T, 
                         init_overdispersion = 10, size_factors = "psinorm", parallel.cores = 1, verbose = TRUE)
  
  contrast = rep(0, ncol(design_matrix))
  contrast[length(contrast)] = 1
  res = devil::test_de(devil.fit = fit, contrast = contrast, clusters = meta$patient, parallel.cores = 1)
  
  tibble(gene = paste0("Gene", 1:nrow(res)),
         pval = res$pval,
         p_adj = p.adjust(res$pval, "BH"),
         lfc = res$lfc,
         method = "devil+se_patients")
}

run_DEVIL_nobatch_sandwichbatches <- function(Y, meta, cluster_id){
  
  Y = as.matrix(Y)
  design_matrix = model.matrix(~ condition, data = meta)
  fit = devil::fit_devil(input_matrix = Y, design_matrix = design_matrix, overdispersion = T, 
                         init_overdispersion = 10, size_factors = "psinorm", parallel.cores = 1, verbose = TRUE)
  
  contrast = rep(0, ncol(design_matrix))
  contrast[length(contrast)] = 1
  res = devil::test_de(devil.fit = fit, contrast = contrast, clusters = meta$batch, parallel.cores = 1)
  
  tibble(gene = paste0("Gene", 1:nrow(res)),
         pval = res$pval,
         p_adj = p.adjust(res$pval, "BH"),
         lfc = res$lfc,
         method = "devil+se_batches")
}
