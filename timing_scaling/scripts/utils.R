
prep_MacaqueBrain_data <- function(N_CELL_TYPES, min_count_per_cell = 100, min_cell_per_gene = 100) {
  data = readRDS("datasets/macaque_brain.rds")
  
  metadata <- data@meta.data
  cell_types = names(sort(table(metadata$cell_type), decreasing = TRUE))[1:N_CELL_TYPES]

  print("Using the following cell types:")
  print(cell_types)

  cell_idx = which(metadata$cell_type %in% cell_types == TRUE)

  metadata <- metadata[cell_idx,]
  metadata$cell_type <- factor(metadata$cell_type, levels = unique(metadata$cell_type))
  #cnt <- as.matrix(data@assays$RNA@counts)[,cell_idx]
  cnt <- data@assays$RNA@counts[,cell_idx]

  cell_idx <- which((colSums(cnt) > min_count_per_cell) == TRUE)
  cnt <- cnt[,cell_idx]

  gene_idx = which((rowSums(cnt) > min_cell_per_gene) == TRUE)
  cnt <- cnt[gene_idx,]

  metadata <- metadata[cell_idx,]
  design_matrix <- model.matrix(~cell_type, metadata)
  donors = metadata$donor_id

  return(list(cnt = cnt, design_matrix = design_matrix, donors = as.factor(donors)))
}

prep_HumanBlood_data <- function(N_CELL_TYPES, min_count_per_cell = 100, min_cell_per_gene = 100) {
  data = readRDS("datasets/HumanBlood/HumanBloodData.rds")

  metadata <- data$meta
  cell_types = names(sort(table(metadata$cell_type), decreasing = TRUE))[1:N_CELL_TYPES]

  print("Using the following cell types:")
  print(cell_types)

  cell_idx = which(metadata$cell_type %in% cell_types == TRUE)

  metadata <- metadata[cell_idx,]
  metadata$cell_type <- factor(metadata$cell_type, levels = unique(metadata$cell_type))
  #cnt <- as.matrix(data@assays$RNA@counts)[,cell_idx]
  cnt <- data$rna[,cell_idx]

  cell_idx <- which((colSums(cnt) >= min_count_per_cell) == TRUE)
  cnt <- cnt[,cell_idx]

  gene_idx = which((rowSums(cnt) >= min_cell_per_gene) == TRUE)
  cnt <- cnt[gene_idx,]

  metadata <- metadata[cell_idx,]
  metadata$Age <- (metadata$Age - mean(metadata$Age)) / sd(metadata$Age)
  design_matrix <- model.matrix(~Age + cell_type, metadata)

  rownames(design_matrix) <- colnames(cnt)

  return(list(cnt = cnt, design_matrix = design_matrix, clusters = metadata$Donor_id))
}

prep_data_small <- function(min_count_per_cell = 100, min_cell_per_gene = 100) {
  data = readRDS("datasets/baronPancreas.rds")
  cnt <- as.matrix(data$counts)
  metadata <- data$metadata
  cell_idx <- (metadata$label %in% c("alpha", "beta"))

  cnt <- cnt[, cell_idx]
  cnt <- cnt[rowSums(cnt > 0) > 100, ]
  metadata <- metadata[cell_idx,]

  design_matrix <- model.matrix(~label, data = metadata)

  return(list(cnt = cnt, design_matrix = design_matrix, donors = factor(metadata$donor)))
}

filter_input <- function(cnt, design_matrix, p_genes, p_cells, n.sub.genes = NULL, n.sub.cells = NULL, clusters = NULL) {
  if (is.null(n.sub.genes) & is.null(n.sub.cells)) {
    n.genes <- dim(cnt)[1]
    n.cells <- dim(cnt)[2]

    n.sub.genes <- as.integer(p_genes * n.genes)
    n.sub.cells <- as.integer(p_cells * n.cells)
  }

  # Find cells with highest expression
  if (n.sub.cells > ncol(cnt)) {
    cell.idxs <- sample(1:ncol(cnt), size = n.sub.cells, replace = TRUE)
  } else {
    cell.idxs <- order(colSums(cnt), decreasing = TRUE)[1:n.sub.cells]
  }
  
  cnt <- cnt[,cell.idxs]
  design_matrix <- design_matrix[cell.idxs,]

  # Find genes with highest expression
  gene.idxs <- order(rowSums(cnt), decreasing = TRUE)[1:n.sub.genes]
  cnt <- cnt[gene.idxs,]

  if (!(is.null(clusters))) {
    clusters = clusters[cell.idxs]
  }

  print(paste0("Lowest RNA counts for a gene = ", min(rowSums(cnt))))
  if (!(is.null(clusters))) {
    list(c=cnt, d = design_matrix, clusters=clusters)
  } else {
    list(c=cnt, d = design_matrix)
  }

}



my_fit_devil <- function(
    input_matrix,
    design_matrix,
    overdispersion = TRUE,
    init_overdispersion = TRUE,
    offset=1e-6,
    size_factors=TRUE,
    verbose=FALSE,
    max_iter=200,
    tolerance=1e-3,
    #eps=1e-6,
    CUDA = TRUE,
    batch_size = 1024L,
    parallel.cores=NULL) {

  max.cores <- parallel::detectCores()
  if (is.null(parallel.cores)) {
    n.cores = max.cores
  } else {
    if (parallel.cores > max.cores) {
      message(paste0("Requested ", parallel.cores, " cores, but only ", max.cores, " available."))
    }
    n.cores = min(max.cores, parallel.cores)
  }

  # Check if CUDA is available
  CUDA_is_available <- FALSE
  if (CUDA) {
    message("Check CUDA availability function need to be implemented")
    CUDA_is_available <- TRUE
  }

  # Add epsilon to input_matrix to avoid non invertible matrices
  # start_time <- Sys.time()
  #
  # input_matrix <- input_matrix + eps
  # input_mat <- devil:::handle_input_matrix(input_matrix, verbose=verbose)
  #
  # end_time <- Sys.time()
  # message("INPUT MATRIX HANDLING:")
  # message(as.numeric(difftime(end_time, start_time, units = "secs")))

  gene_names <- rownames(input_matrix)
  # counts_per_cell <- rowMeans(input_mat)
  # cell_per_genes <- rowSums(input_mat > 0)
  # filter_genes <- (counts_per_cell <= avg_counts) | (cell_per_genes <= min_cells)
  # n_low_genes <- sum(filter_genes)
  # if (n_low_genes > 0) {
  #   message(paste0("Removing ", n_low_genes, " lowly expressed genes."))
  #   input_mat <- matrix(input_mat[!filter_genes, ], ncol = nrow(design_matrix), nrow = sum(!filter_genes))
  #   input_matrix <- input_matrix[!filter_genes, ]
  #   gene_names <- gene_names[!filter_genes]
  # }

  start_time <- Sys.time()
  if (size_factors) {
    if (verbose) { message("Compute size factors") }
    sf <- devil:::calculate_sf(input_matrix, verbose = verbose)
  } else {
    sf <- rep(1, nrow(design_matrix))
  }
  end_time <- Sys.time()
  message("Size factor computing:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  start_time <- Sys.time()
  offset_matrix = devil:::compute_offset_matrix(offset, input_matrix, sf)
  #offset_matrix = devil:::compute_offset_matrix(offset, input_matrix, sf)[1,]
  #offset_matrix = matrix(offset_matrix, nrow = 1)
  end_time <- Sys.time()
  message("Offset matrix computing:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  estimate_dispersion <- function (y, offset_matrix) {
    xim <- 1/mean(DelayedMatrixStats::colMeans2(exp(offset_matrix),
                                                useNames = T))
    bv <- DelayedMatrixStats::rowVars(y, useNames = T)
    bm <- DelayedMatrixStats::rowMeans2(y, useNames = T)
    disp <- (bv - xim * bm)/bm^2
    ifelse(is.na(disp) | disp < 0, 100, disp)
  }

  start_time <- Sys.time()
  if (init_overdispersion) {
    dispersion_init <- c(estimate_dispersion(input_matrix, offset_matrix))
  } else {
    dispersion_init <- rep(100, nrow(input_mat))
  }

  end_time <- Sys.time()
  message("Dispersion init:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  ngenes <- nrow(input_matrix)
  nfeatures <- ncol(design_matrix)



  start_time <- Sys.time()
  if (verbose) { message("Initialize beta estimate") }
  groups <- devil:::get_groups_for_model_matrix(design_matrix)
  end_time <- Sys.time()
  message("get_groups_for_model_matrix:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  init_beta <- function (y, design_matrix, offset_matrix) {
    qrx <- qr(design_matrix)
    Q <- qr.Q(qrx)[seq_len(nrow(design_matrix)), , drop = FALSE]
    R <- qr.R(qrx)
    norm_log_count_mat <- t(log1p((y/exp(offset_matrix[1, ]))))
    t(solve(R, as.matrix(t(Q) %*% norm_log_count_mat)))
  }

  init_beta_groups <- function(y, groups, offset_matrix) {
    #norm_Y <- y / exp(offset_matrix)
    norm_Y <- y / exp(offset_matrix[1,])
    do.call(cbind, lapply(unique(groups), function(gr){
      log(DelayedMatrixStats::rowMeans2(norm_Y, cols = groups == gr, useNames=TRUE))
    }))
  }

  start_time <- Sys.time()
  # if (!is.null(groups)) {
  #   beta_0 <- devil:::init_beta_groups(input_matrix, groups, offset_matrix)
  #   #beta_0 <- devil:::init_beta(input_mat, design_matrix, offset_matrix)
  # } else {
  #   beta_0 <- devil:::init_beta(input_matrix, design_matrix, offset_matrix)
  # }
  beta_0 <- init_beta(input_matrix, design_matrix, offset_matrix)
  #beta_0 <- init_beta_groups(input_matrix, groups, offset_matrix)
  end_time <- Sys.time()
  message("init_beta:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  #beta_0 <- devil:::init_beta(input_mat, design_matrix, offset_matrix)

  if (CUDA & CUDA_is_available) {
    message("Messing with CUDA! Implementation still needed")

    remainder = ngenes %% batch_size
    extra_genes = remainder
    genes_batch = ngenes - extra_genes

    message("Fit beta CUDA")

    start_time <- Sys.time()
    res_beta_fit <- devil:::beta_fit_gpu(
      input_matrix[1:genes_batch,],
      design_matrix,
      beta_0[1:genes_batch,],
      offset_matrix[1:genes_batch,],
      dispersion_init[1:genes_batch],
      max_iter = max_iter,
      eps = tolerance,
      batch_size = batch_size
    )
    res_beta_fit_extra <- devil:::beta_fit_gpu(
      input_matrix[(genes_batch):ngenes,],
      design_matrix,
      beta_0[(genes_batch):ngenes,],
      offset_matrix[(genes_batch):ngenes,],
      dispersion_init[(genes_batch):ngenes],
      max_iter = max_iter,
      eps = tolerance,
      batch_size = batch_size
    )
    end_time <- Sys.time()
    message("BETA GPU RUNTIME:")
    message(as.numeric(difftime(end_time, start_time, units = "secs")))

    beta = res_beta_fit$mu_beta
    beta_extra = res_beta_fit_extra$mu_beta
    beta <- rbind(beta, beta_extra)
    if (is.null(dim(beta))) {
      beta = matrix(beta, ncol = 1)
    }
    
    iterations=c(res_beta_fit$iter, res_beta_fit_extra$iter)

  } else {
    start_time <- Sys.time()
    if (verbose) { message("Fit beta coefficients") }
    tmp <- parallel::mclapply(1:ngenes, function(i) {
      devil:::beta_fit(input_mat[i,], design_matrix, beta_0[i,], offset_matrix[1,], dispersion_init[i], max_iter = max_iter, eps = tolerance)
      #devil:::beta_fit(input_mat[i,], design_matrix, beta_0[i,], offset_matrix[i,], 1, max_iter = max_iter, eps = tolerance)
    }, mc.cores = n.cores)
    beta <- lapply(1:ngenes, function(i) { tmp[[i]]$mu_beta }) %>% do.call("rbind", .)
    rownames(beta) <- gene_names
    end_time <- Sys.time()
    message("BETA CPU RUNTIME:")
    message(as.numeric(difftime(end_time, start_time, units = "secs")))

    iterations <- lapply(1:ngenes, function(i) { tmp[[i]]$iter }) %>% unlist()
  }


  if (overdispersion) {
    if (verbose) { message("Fit overdispersion") }

    theta <- parallel::mclapply(1:ngenes, function(i) {
      devil:::fit_dispersion(beta[i,], design_matrix, input_matrix[i,], offset_matrix[i,], tolerance = tolerance, max_iter = max_iter)
    }, mc.cores = n.cores) %>% unlist()

  } else {
    theta = rep(0, ngenes)
  }

  return(
    list(
      beta=beta,
      overdispersion=theta,
      iterations=iterations,
      size_factors=sf,
      offset_matrix=offset_matrix,
      design_matrix=design_matrix,
      input_matrix=input_matrix,
      input_parameters=list(max_iter=max_iter, tolerance=tolerance, parallel.cores=n.cores)
    )
  )
}


make_example_dds <- function(n_genes, n_replicates, n_conditions){

  dispMeanRel <-  function(x) 4/x + 0.1

  beta <- matrix(rnorm(n_genes*n_conditions, mean = 4, sd = 2), ncol = n_conditions)
  dispersion <- dispMeanRel(2^(beta[, 1]))
  colData <- data.frame(condition = factor(rep(paste0("cond", 1:n_conditions), n_replicates)))
  if (n_conditions == 1) {
    x <- matrix(data = 1, nrow = n_replicates, ncol = 1)
  } else {
    x <- model.matrix.default(~colData$condition)
  }

  mu <- t(2^(x %*% t(beta)))

  countData <- matrix(rnbinom(mu, mu = mu, size = 1/dispersion), ncol = ncol(mu))
  mode(countData) <- "integer"

  #design <- as.formula("~ condition", env = .GlobalEnv)
  object <- list(countData = countData, design = x)
  object
}

fit_devil_gpu <- function(
    input_matrix,
    design_matrix,
    overdispersion = TRUE,
    init_overdispersion = NULL,
    do_cox_reid_adjustment = TRUE,
    offset=1e-6,
    size_factors=TRUE,
    verbose=FALSE,
    max_iter=200,
    tolerance=1e-3,
    CUDA = FALSE,
    batch_size = 1024L,
    parallel.cores=NULL) {

  # Read general info about input matrix and design matrix
  gene_names <- rownames(input_matrix)
  ngenes <- nrow(input_matrix)
  nfeatures <- ncol(design_matrix)

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

  # Check if CUDA is available
  CUDA_is_available <- FALSE
  if (CUDA) {
    message("Check CUDA availability function need to be implemented")
    CUDA_is_available <- TRUE
  }

  # Compute size factors
  if (size_factors) {
    if (verbose) { message("Compute size factors") }
    sf <- devil:::calculate_sf(input_matrix, verbose = verbose)
  } else {
    sf <- rep(1, nrow(design_matrix))
  }

  # Calculate offset vector
  offset_vector = devil:::compute_offset_vector(offset, input_matrix, sf)
  #offset_matrix = devil:::compute_offset_matrix(offset, input_matrix, sf)

  # Initialize overdispersion
  if (is.null(init_overdispersion)) {
    dispersion_init <- c(devil:::estimate_dispersion(input_matrix, offset_vector))
  } else {
    dispersion_init <- rep(init_overdispersion, nrow(input_matrix))
  }

  # if (verbose) { message("Initialize beta estimate") }
  # groups <- devil:::get_groups_for_model_matrix(design_matrix)

  if (verbose) { message("Initialize beta estimate") }
  beta_0 <- devil:::init_beta(input_matrix, design_matrix, offset_vector)

  if (CUDA & CUDA_is_available) {
    message("Messing with CUDA! Implementation still needed")

    remainder = ngenes %% batch_size
    extra_genes = remainder
    genes_batch = ngenes - extra_genes

    message("Fit beta CUDA")
    start_time <- Sys.time()

    res_beta_fit <- devil:::beta_fit_gpu(
      input_matrix[1:genes_batch,],
      design_matrix,
      beta_0[1:genes_batch,],
      offset_vector,
      dispersion_init[1:genes_batch],
      max_iter = max_iter,
      eps = tolerance,
      batch_size = batch_size
    )

    if (remainder > 0) {
      res_beta_fit_extra <- devil:::beta_fit_gpu(
        input_matrix[(genes_batch+1):ngenes,],
        design_matrix,
        beta_0[(genes_batch+1):ngenes,],
        offset_vector,
        dispersion_init[(genes_batch+1):ngenes],
        max_iter = max_iter,
        eps = tolerance,
        batch_size = extra_genes
      )
    }

    end_time <- Sys.time()
    message("BETA GPU RUNTIME:")
    message(as.numeric(difftime(end_time, start_time, units = "secs")))

    beta = res_beta_fit$mu_beta

    if (remainder > 0) {
      beta_extra = res_beta_fit_extra$mu_beta
      beta <- rbind(beta, beta_extra)
      iterations=c(res_beta_fit$iter, res_beta_fit_extra$iter)
    } else {
      iterations=c(res_beta_fit$iter)
    }

    if (is.null(dim(beta))) {
      beta = matrix(beta, ncol = 1)
    }

    bad_idx = which(rowSums(is.na(beta)) > 0)
    tmp = lapply(bad_idx, function(i) {
      beta[i,] <<- devil:::beta_fit(input_matrix[i,], design_matrix, beta_0[i,], offset_vector, dispersion_init[i], max_iter = max_iter, eps = tolerance)$mu_beta
    })

    rm(tmp)

  } else {

    if (verbose) { message("Fit beta coefficients") }

    tmp <- parallel::mclapply(1:ngenes, function(i) {
      devil:::beta_fit(input_matrix[i,], design_matrix, beta_0[i,], offset_vector, dispersion_init[i], max_iter = max_iter, eps = tolerance)
    }, mc.cores = n.cores)

    beta <- lapply(1:ngenes, function(i) { tmp[[i]]$mu_beta }) %>% do.call("rbind", .)
    rownames(beta) <- gene_names
    iterations <- lapply(1:ngenes, function(i) { tmp[[i]]$iter }) %>% unlist()

  }

  s <- Sys.time()
  if (overdispersion) {
    if (verbose) { message("Fit overdispersion") }

    theta <- parallel::mclapply(1:ngenes, function(i) {
      devil:::fit_dispersion(beta[i,], design_matrix, input_matrix[i,], offset_vector, tolerance = tolerance, max_iter = max_iter, do_cox_reid_adjustment = do_cox_reid_adjustment)
    }, mc.cores = n.cores) %>% unlist()

  } else {
    theta = rep(0, ngenes)
  }

  return(list(
    beta=beta,
    overdispersion=theta,
    iterations=iterations,
    size_factors=sf,
    offset_vector=offset_vector,
    design_matrix=design_matrix,
    input_matrix=input_matrix,
    input_parameters=list(max_iter=max_iter, tolerance=tolerance, parallel.cores=n.cores)
  )
  )
}

get_results <- function(results_folder) {
  results_paths <- list.files(results_folder)
  results_paths <- results_paths[results_paths != "fits"]

  results <- lapply(results_paths, function(p) {
    print(p)
    info <- unlist(strsplit(p, "_"))
    res = readRDS(paste0(results_folder, p))


    dplyr::tibble(
      model_name = paste(info[1], info[2]),
      p_gene = info[3],
      p_cells = info[5],
      n_cell_types = info[7],
      time = as.numeric(unlist(res$time), units = "secs"),
      memory = res$mem_alloc
    )
  }) %>% do.call("bind_rows", .) %>%
    dplyr::mutate(p_cells = as.numeric(p_cells))
}

plot_time_and_memory_comparison = function(results) {
  p1 <- results %>%
    dplyr::group_by(p_gene, p_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time), sd =sd(time)) %>%
    ggplot(mapping = aes(x = p_cells, y = y, ymin=y-sd, ymax=y+sd, col = model_name)) +
    geom_pointrange() +
    geom_line() +
    ggh4x::facet_nested(~"Gene percentage"+p_gene, scales = "free_x") +
    theme_bw() +
    labs(x = "Percentage of cells", y = "Time (s)", col = "Model")

  p2 <- results %>%
    ggplot(mapping = aes(x = p_cells, y = memory * 1e-9, col=model_name)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"Gene percentage"+p_gene, scales = "free_x") +
    theme_bw() +
    labs(x = "Cells percentage", y = "Memory (GB)", fill = "Model")

  p3 <- results %>%
    dplyr::group_by(p_gene, p_cells, n_cell_types, model_name) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::group_by(p_gene, p_cells, n_cell_types) %>%
    dplyr::mutate(n = sum(n)) %>%
    #dplyr::ungroup() %>%
    #dplyr::filter(n == max(n)) %>%
    dplyr::group_by(p_gene, p_cells, n_cell_types) %>%
    dplyr::mutate(ratio_time = time / time[model_name == "gpu devil"]) %>%
    dplyr::group_by(p_gene, p_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(ratio_time), sd =sd(ratio_time)) %>%
    ggplot(mapping = aes(x = p_cells, y = y, ymin=y-sd, ymax=y+sd, col = model_name)) +
    geom_pointrange() +
    geom_line() +
    ggh4x::facet_nested(~"Gene percentage"+p_gene, scales = "free_x") +
    theme_bw() +
    labs(x = "Percentage of cells", y = "Time ratio", col = "Model")

  list(time=p1, memory=p2, ratio_time=p3)
}


plot_correlations <- function(fits_folder) {
  fits <- list.files(fits_folder, full.names = T)

  devil.res <- readRDS(fits[grepl("/cpu_devil_", fits)])
  gpu.devil.res <- readRDS(fits[grepl("/gpu_devil_", fits)])
  gpu.devilnotheta.res <- readRDS(fits[grepl("/gpu_devilnotheta_", fits)])
  gpu.devilnoinit.res <- readRDS(fits[grepl("/gpu_noinitdevil", fits)])
  glm.res <- readRDS(fits[grepl("/cpu_glmGam", fits)])

  p1 <- dplyr::bind_rows(
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = devil.res$lfc,
      y_name = "cpu devil"
    ),
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = glm.res$lfc,
      y_name = "cpu glmGamPoi"
    ),
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = gpu.devilnotheta.res$lfc,
      y_name = "gpu no theta"
    ),
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = gpu.devilnoinit.res$lfc,
      y_name = "gpu no init"
    )
  ) %>%
    dplyr::group_by(x_name, y_name) %>%
    dplyr::filter(abs(x) < 10 & abs(y) < 10) %>%
    ggplot(mapping = aes(x=x, y=y)) +
    geom_point() +
    ggpubr::stat_cor() +
    facet_grid(x_name~y_name, switch = "both") +
    theme_bw() +
    labs(x = bquote(LFC[1]), y=bquote(LFC[2])) +
    ggtitle("Log fold change correlation")
  p1

  p2 <- dplyr::bind_rows(
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = devil.res$theta,
      y_name = "cpu devil"
    ),
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = glm.res$theta,
      y_name = "cpu glmGamPoi"
    ),
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = gpu.devilnotheta.res$theta,
      y_name = "gpu no theta"
    ),
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = gpu.devilnoinit.res$theta,
      y_name = "gpu no init"
    )
  ) %>%
    dplyr::group_by(x_name, y_name) %>%
    ggplot(mapping = aes(x=x, y=y)) +
    geom_point() +
    ggpubr::stat_cor() +
    facet_grid(x_name~y_name, switch = "both") +
    theme_bw() +
    labs(x = bquote(theta[1]), y=bquote(theta[2])) +
    ggtitle("Overdisperions correlation")
  p2

  # p3 <- dplyr::bind_rows(
  #   dplyr::tibble(
  #     x = gpu.devil.res$pval,
  #     x_name = "gpu devil",
  #     y = devil.res$pval,
  #     y_name = "cpu devil"
  #   ),
  #   dplyr::tibble(
  #     x = gpu.devil.res$pval,
  #     x_name = "gpu devil",
  #     y = glm.res$pval,
  #     y_name = "cpu glmGamPoi"
  #   ),
  #   dplyr::tibble(
  #     x = gpu.devil.res$pval,
  #     x_name = "gpu devil",
  #     y = gpu.devilnotheta.res$pval,
  #     y_name = "gpu no theta"
  #   ),
  #   dplyr::tibble(
  #     x = gpu.devil.res$pval,
  #     x_name = "gpu devil",
  #     y = gpu.devilnoinit.res$pval,
  #     y_name = "gpu no init"
  #   )
  # ) %>%
  #   dplyr::group_by(x_name, y_name) %>%
  #   ggplot(mapping = aes(x=x, y=y)) +
  #   geom_point() +
  #   ggpubr::stat_cor() +
  #   facet_grid(x_name~y_name, switch = "both") +
  #   theme_bw() +
  #   labs(x = bquote(p[1]), y=bquote(p[2])) +
  #   ggtitle("P-values correlation")

  list(lfc=p1, theta=p2)
}
