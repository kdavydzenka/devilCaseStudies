sim_treatment_clusters <- function(
    ngenes = 1000,
    cells_per_patient = 2000,
    n_batches = 4,
    patients_per_batch = 3,
    # batches (nuisance)
    batch_lfc_sd = 0.4,
    libsize_cv   = 0.4,
    disp_mult_sd = 0.5,
    # treatment assignment mode
    assignment = c("within_patient","by_patient"),
    # pattern of treatment vs batch (only used for by_patient)
    batch_treat_pattern = c("balanced", "confounded"),
    # link cluster labels to treatment? (TRUE => clusters are literally Untreated/Treated)
    link_cluster_to_treatment = FALSE,
    # treatment DE settings (applies within each *responding* cluster)
    de_prop_treat   = 0.15,
    lfc_treat_loc   = log(1.5),   # on log scale (ln FC)
    lfc_treat_sd    = 0.2,
    # optional: baseline cluster markers via Splatter (kept small)
    de_prop_cluster = 0.05,
    # fraction of clusters that truly respond to treatment (e.g. 1.0, 0.5, ...)
    prop_clusters_respond = 1.0,
    seed = 1
){
  set.seed(seed)
  assignment         <- match.arg(assignment)
  batch_treat_pattern <- match.arg(batch_treat_pattern)
  
  n_clusters <- 2L
  n_patients <- n_batches * patients_per_batch
  cells_total <- n_patients * cells_per_patient
  
  # --- baseline via Splatter: 2 clusters + batches ---
  params <- newSplatParams()
  params <- setParam(params, "nGenes", ngenes)
  params <- setParam(params, "batchCells", rep(cells_total / n_batches, n_batches))
  params <- setParam(params, "batch.facLoc", 0)                 # mean 1
  params <- setParam(params, "batch.facScale", libsize_cv)      # lib-size variability
  params <- setParam(params, "group.prob", rep(1 / n_clusters, n_clusters))
  params <- setParam(params, "de.prob", de_prop_cluster)        # faint cluster markers
  params <- setParam(params, "de.downProb", 0.5)
  params <- setParam(params, "de.facLoc", 0.15)
  params <- setParam(params, "de.facScale", 0.15)
  
  sim <- splatSimulate(params, method = "groups",
                       batch.rmEffect = FALSE, verbose = FALSE)
  
  # --- batch effects (mean + dispersion) ---
  B <- n_batches
  G <- ngenes
  
  batch_effects <- matrix(
    rnorm(G * B, mean = 0, sd = batch_lfc_sd),
    nrow = G, ncol = B
  )
  disp_mult <- exp(rnorm(B, mean = 0, sd = disp_mult_sd))
  
  batch <- as.character(SummarizedExperiment::colData(sim)$Batch)
  batch_index <- as.integer(factor(batch, levels = sort(unique(batch))))
  stopifnot(length(unique(batch_index)) == B)
  
  mu_counts <- SummarizedExperiment::assays(sim)$counts
  libsz <- Matrix::colSums(mu_counts); libsz[libsz == 0] <- 1
  log_mu <- log1p(mu_counts / rep(libsz, each = nrow(mu_counts)))
  
  for (b in seq_len(B)) {
    cols_b <- which(batch_index == b)
    log_mu[, cols_b] <- log_mu[, cols_b] + batch_effects[, b]
  }
  
  # --- patient ids nested in batches ---
  # Patients are nested inside batches: patients_per_batch per batch.
  patient <- paste0("P", rep(seq_len(n_patients), each = cells_per_patient))
  batch_of_patient <- rep(seq_len(n_batches), each = patients_per_batch)
  batch_by_cell <- batch_of_patient[as.integer(sub("P", "", patient))]
  
  # check nesting consistency: Splatter batches vs our patient -> batch mapping
  stopifnot(all(batch_index == batch_by_cell))
  
  # --- clusters from Splatter (2 levels) ---
  cluster <- as.character(SummarizedExperiment::colData(sim)$Group)
  
  # --- assign treatment according to mode & batch_treat_pattern ---
  condition <- rep(NA_character_, length(cluster))
  
  if (assignment == "within_patient") {
    # each patient has both A and B cells
    for (p in unique(patient)) {
      idx <- which(patient == p)
      nB  <- length(idx) %/% 2
      condition[idx] <- c(rep("B", nB), rep("A", length(idx) - nB))
      condition[idx] <- sample(condition[idx])  # mix within patient
    }
  } else { # by_patient
    pts <- unique(patient)
    
    if (batch_treat_pattern == "balanced") {
      # random half of patients treated
      pts_B <- sample(pts, length(pts) %/% 2)
    } else { # confounded: later batches mostly treated
      batch_ids <- seq_len(n_batches)
      treated_batches <- batch_ids > (n_batches / 2)
      # get a batch for each patient (via first occurrence)
      patient_index <- as.integer(sub("P", "", pts))
      patient_batches <- batch_of_patient[patient_index]
      pts_B <- pts[patient_batches %in% batch_ids[treated_batches]]
      if (length(pts_B) == 0L || length(pts_B) == length(pts)) {
        # fallback if something degenerate happens
        pts_B <- sample(pts, length(pts) %/% 2)
      }
    }
    
    condition <- ifelse(patient %in% pts_B, "B", "A")
  }
  
  condition <- factor(condition, levels = c("A","B"))
  
  # Optionally make clusters == treatment labels (usually keep FALSE for your sims)
  if (link_cluster_to_treatment) {
    cluster <- ifelse(condition == "B", "Treated", "Untreated")
  }
  
  cluster_levels <- levels(factor(cluster))
  
  # --- choose which clusters truly respond to treatment ---
  n_resp <- max(1L, round(length(cluster_levels) * prop_clusters_respond))
  responding_clusters <- sample(cluster_levels, size = n_resp, replace = FALSE)
  
  # --- add TRUE treatment DE within *responding* clusters only ---
  G_de <- round(de_prop_treat * G)
  de_genes <- if (G_de > 0L) sample.int(G, G_de, replace = FALSE) else integer(0)
  
  lfc_g <- if (G_de > 0L) {
    rnorm(G_de, mean = lfc_treat_loc, sd = lfc_treat_sd) *
      sample(c(-1, 1), G_de, replace = TRUE)
  } else {
    numeric(0)
  }
  
  # truth_treat: gene × cluster table
  truth_treat <- tidyr::expand_grid(
    gene = rownames(mu_counts),
    cluster = cluster_levels
  ) %>%
    dplyr::mutate(
      is_de_treatment = FALSE,
      lfc_treatment   = 0
    )
  
  for (cl in cluster_levels) {
    cells_cl_B <- which(cluster == cl & condition == "B")
    
    if (cl %in% responding_clusters && length(cells_cl_B) > 0 && length(de_genes) > 0) {
      # apply log-FC shift only in responding clusters, treated cells
      log_mu[de_genes, cells_cl_B] <- sweep(
        log_mu[de_genes, cells_cl_B, drop = FALSE],
        1, lfc_g, `+`
      )
      
      # update truth table for this cluster
      rows_cl <- which(truth_treat$cluster == cl)
      gene_rows <- rows_cl[de_genes]
      truth_treat$is_de_treatment[gene_rows] <- TRUE
      truth_treat$lfc_treatment[gene_rows]   <- lfc_g / log(2)
    }
  }
  
  # --- finalize mean & simulate counts with batch-specific dispersion ---
  mu_adj <- expm1(log_mu) * rep(libsz, each = nrow(mu_counts))
  mu_adj[mu_adj < 0] <- 0
  
  counts <- matrix(0, nrow = nrow(mu_counts), ncol = ncol(mu_counts))
  alpha_base <- 1 / 0.5  # base dispersion (can be tuned)
  
  for (b in seq_len(B)) {
    cols_b <- which(batch_index == b)
    m_b <- mu_adj[, cols_b, drop = FALSE]
    alpha_b <- alpha_base / disp_mult[b]
    counts[, cols_b] <- matrix(
      rnbinom(length(m_b), size = alpha_b, mu = as.vector(m_b)),
      nrow = nrow(m_b)
    )
  }
  
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  colnames(counts) <- colnames(mu_counts)
  rownames(counts) <- rownames(mu_counts)
  
  # --- truth: optional cluster markers from Splatter ---
  truth_cluster <- SummarizedExperiment::rowData(sim) %>%
    as.data.frame() %>%
    dplyr::mutate(gene = rownames(sim)) %>%
    dplyr::select(gene, dplyr::starts_with("DEFac")) %>%
    tidyr::pivot_longer(-gene, names_to = "decol", values_to = "fac") %>%
    dplyr::mutate(
      cluster = stringr::str_remove(decol, "^DEFac"),
      cluster = dplyr::if_else(
        stringr::str_detect(cluster, "^[0-9]+$"),
        levels(factor(SummarizedExperiment::colData(sim)$Group))[as.integer(cluster)],
        cluster
      ),
      is_de_cluster = !is.na(fac) & fac != 1,
      lfc_cluster   = dplyr::if_else(is.na(fac), NA_real_, log2(fac))
    ) %>%
    dplyr::select(gene, cluster, is_de_cluster, lfc_cluster)
  
  # --- meta data ---
  meta <- tibble::tibble(
    cell      = colnames(counts),
    batch     = paste0("B", batch_index),
    patient   = patient,
    cluster   = factor(cluster, levels = cluster_levels),
    condition = condition
  )
  
  list(
    counts          = counts,
    meta            = meta,
    truth_treatment = truth_treat,   # gene × cluster truth for treatment effect
    truth_cluster   = truth_cluster  # (optional) Splatter cluster markers
  )
}
