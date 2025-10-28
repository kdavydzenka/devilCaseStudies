#DATASET_NAMES <- c("MuscleRNA", "MuscleATAC")

read_data <- function(dataset_name, data_path) {
  if (dataset_name == "MuscleRNA") {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    tissue <- "muscle"
  } else if (dataset_name == "MuscleATAC") {
    atac <- readRDS(data_path)
    counts <- atac@assays@data@listData[["PeakMatrix"]]
    granges <- atac@rowRanges
    granges %>% mutate(ranges = paste(granges$seqname,granges$start,granges$end, sep = ":"))
    rownames(peak_counts) <- granges$ranges
    metadata <- atac@colData %>% as.data.frame()
    tissue <- "muscle"
  }  else {
    stop("Dataset name not recognized")
  }

  return(list(counts=counts, metadata=metadata, tissue=tissue))
}

prepare_rna_input <- function(input_data) {
  metadata <- input_data$metadata
  metadata <- metadata[ (metadata$tech %in% c("snRNA") & metadata$Annotation %in% c("Type I", "Type II")) , ]
  metadata <- metadata %>% dplyr::rename(cell_type = Annotation)
  metadata$cell_type <- as.factor(metadata$cell_type)
  metadata <- metadata %>%
    mutate(age_cluster = case_when(
      age_pop == "old_pop"  ~ '1',
      age_pop == "young_pop" ~ '0'
    ))
  metadata$age_cluster <- as.factor(metadata$age_cluster)
  counts <- input_data$counts
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)
  
  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  feat_mad_filter <- total_features > 5 * mad(total_features)
  
  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .1
  cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter #|  mit_prop_filter
  
  counts <- counts[, !cell_outliers_filter]
  metadata <- metadata[!cell_outliers_filter, ]
  
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  tissue = "muscle"
  return(list(counts=counts, metadata=metadata, tissue=tissue))
}


subsample_balanced_cells <- function(input_data, seed = 123) {
  set.seed(seed)

  metadata <- input_data$metadata
  counts <- input_data$counts

  # Check that required columns exist
  if (!all(c("age_pop", "cell_type") %in% colnames(metadata))) {
    stop("metadata must contain 'age_pop' and 'cell_type' columns")
  }

  # Determine smallest group size across (age_cluster, cell_type)
  group_sizes <- metadata %>%
    dplyr::count(age_cluster, cell_type)

  min_cells <- min(group_sizes$n)
  message("Balancing to ", min_cells, " cells per (age_cluster × cell_type) group")

  # Subsample each group to the same number of cells
  balanced_metadata <- metadata %>%
    dplyr::group_by(age_cluster, cell_type) %>%
    dplyr::sample_n(min_cells) %>%
    dplyr::ungroup()
  
  rownames(balanced_metadata) <- balanced_metadata$cell_index

  # Subset count matrix
  balanced_counts <- counts[, rownames(balanced_metadata), drop = FALSE]

  list(
    counts = balanced_counts,
    metadata = balanced_metadata,
    tissue = input_data$tissue
  )
}


fit_de <- function(input_data, 
                   method = "devil", 
                   design_type = "age_only"
                   ) {

  if (!(method %in% c('devil', "glmGamPoi", 'nebula'))) {
    stop('method not recognized')
  }

  metadata <- input_data$metadata
  
  counts <- as.matrix(input_data$counts)
  counts <- round(counts)

  ## Design Matrix
  design_matrix <- switch(
    design_type,
    "age_only" = model.matrix(~ age_cluster, metadata),
    "interaction" = model.matrix(~ age_cluster * cell_type, metadata),
    stop("design_type must be one of: 'age_only', 'interaction'")
  )

  ## Fit coefficients
  if (method == 'devil') {
    
    fit <- devil::fit_devil(
      input_matrix = counts,
      design_matrix = design_matrix,
      overdispersion = TRUE,
      offset = 1e-6,
      size_factors = "normed_sum",
      parallel.cores = 1,
      verbose = TRUE, 
      init_overdispersion = NULL, 
      max_iter = 500, 
      tolerance = 1e-3, 
      batch_size = 1
    )
    

  } else if (method == "glmGamPoi") {

    fit <- glmGamPoi::glm_gp(counts, design_matrix, size_factors = TRUE, verbose = TRUE)
    
  } else if (method == 'nebula') {

    metadata$patient <- as.numeric(as.factor(metadata$sample))
    sf <- devil:::calculate_sf(counts)
    #data_g = group_cell(count = counts, id = metadata$patient, pred = design_matrix, offset = sf)
    #print(str(data_g))
    #fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, offset = data_g$offset)
    
    fit <- nebula::nebula(counts, id = metadata$patient, pred = design_matrix, offset = sf, covariance = TRUE)
  
  }

  fit
}


de_test <- function(input_data,
                    fit_res, 
                    method = "devil", 
                    design_test = "age_type1") {
 
  # Check method validity
  if (!(method %in% c('devil', 'glmGamPoi', 'nebula'))) {
    stop('Method not recognized. Choose one of: devil, glmGamPoi, nebula.')
  }

  metadata <- input_data$metadata
  res_de <- NULL  # initialize result

  # Define contrasts based on design
  if (design_test == "age_type1") {
    contrast <- c(0, 1, 0, 0)
  } else if (design_test == "age_type2") {
    contrast <- c(0, 1, 0, 1)
  } else if (design_test == "interaction") {
    contrast <- c(0, 0, 0, 1)
  } else {
    stop("design_test not recognized")
  }

  # DE test by method
  if (method == 'devil') {

    clusters <- as.numeric(as.factor(metadata$sample))

    res_de <- devil::test_de(
      fit_res,
      contrast = contrast,
      clusters = clusters,
      max_lfc = 10
    )

  } else if (method == 'glmGamPoi') {
    
    res_de <- glmGamPoi::test_de(
      fit_res,
      contrast = contrast
    )

    # Keep only key columns if available
    if (all(c("name", "pval", "adj_pval", "lfc") %in% colnames(res_de))) {
      res_de <- dplyr::select(res_de, name, pval, adj_pval, lfc)
    }

  } else if (method == 'nebula') {
    
    # Nebula-like structure: fit_res$summary and fit_res$covariance
    n_genes <- nrow(fit_res$summary)
    effects <- numeric(n_genes)
    p_values <- numeric(n_genes)
    
    n_variables = sum(grepl("logFC", colnames(nebula_fit$summary)))
    
    if (length(contrast) != n_variables) {
      stop("Passed contrast with wrong number of elements.")
    }
    
    for (gene_i in seq_len(n_genes)) {
      # Build covariance matrix
      cov <- matrix(NA, n_variables, n_variables)
      cov[lower.tri(cov, diag = TRUE)] <- as.numeric(nebula_fit$covariance[gene_i, ])
      cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
      
      # Compute the contrast effect
      eff <- sum(contrast * nebula_fit$summary[gene_i, 1:n_variables])
      
      # Compute p-value
      p <- pchisq(eff^2 / (t(contrast) %*% cov %*% contrast), df = 1, lower.tail = FALSE)
      
      effects[gene_i] <- eff
      p_values[gene_i] <- p
    }
    
    # Add results to summary
    fit_res$summary$lfc <- effects
    fit_res$summary$pval <- p_values
    fit_res$summary$adj_pval <- p.adjust(p_values, method = "BH")  
    
    res_de <- dplyr::select(fit_res$summary, gene = 1, lfc, pval, adj_pval)
  }

  return(res_de)
}



#grange_annot <- function(input_data, data_path) {
  #counts <- input_data$counts
  #metadata <- input_data$metadata
  #grange <- input_data$grange
  #txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  #grange_annot <- ChIPseeker::annotatePeak(
    #peak = grange,
    #tssRegion = c(-3000, 3000),
    #TxDb = txdb,
    #level = "gene",
    #assignGenomicAnnotation = TRUE,
    #genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                 # "Downstream", "Intergenic"),
    #annoDb = "org.Hs.eg.db",
    #addFlankGeneInfo = FALSE,
    #flankDistance = 5000,
    #sameStrand = FALSE,
    #ignoreOverlap = FALSE,
    #ignoreUpstream = FALSE,
    #ignoreDownstream = FALSE,
    #overlap = "TSS",
    #verbose = TRUE,
    #columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
  #) %>% as.data.frame()
  #tissue = "muscle"
  #return(list(counts=counts, metadata=metadata, grange=grange_annot, tissue=tissue))
#}

#prepare_atac_input <- function(input_data) {
  #metadata <- input_data$metadata
  #grange <- input_data$grange_annot
  #metadata <- metadata[ (metadata$group %in% c("young", "old") & metadata$cell_type %in% c("Type I", "Type II")),]
  #metadata <- metadata %>%
    #mutate(age_cluster = case_when(
      #group == "old"  ~ '1',
      #group == "young" ~ '0'
    #))
  #peak_counts <- input_data$counts
  #peak_counts <- peak_counts[ ,colnames(peak_counts) %in% rownames(metadata) ]
  #grange_annot <- grange_annot[ (grange_annot$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")), ]
  #grange_annot <- grange_annot %>% mutate(ranges = paste(grange_annot$seqnames, grange_annot$start, grange_annot$end))
  #peak_counts <- peak_counts[ rownames(peak_counts) %in% grange_annot$ranges , ]
  #total_counts <- colSums(peak_counts)
  #total_features <- colSums(peak_counts > 0)

  #mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  #feat1000_filter <- total_features < 1000
  #feat_mad_filter <- total_features > 5 * mad(total_features)

  #cell_outliers_filter <- mad5_filter | feat1000_filter | feat_mad_filter

  #peak_counts <- peak_counts[, !cell_outliers_filter]
  #metadata <- metadata[!cell_outliers_filter, ]

  #non_expressed_genes <- rowMeans(peak_counts) <= 0.01
  #peak_counts <- peak_counts[!non_expressed_genes, ]
  #peak_counts <- peak_counts[,colnames(peak_counts) %in% rownames(metadata)]
  #tissue = "muscle"
  #return(list(counts=counts, metadata=metadata, grange=grange_annot, tissue=tissue))
#}


#perform_analysis_atac <- function(input_data, method = "devil") {
  #if (!(method %in% c('devil', "glmGamPoi", 'nebula'))) {stop('method not recognized')}

  #if (method == 'devil') {
    #metadata <- input_data$metadata
    #peak_counts <- as.matrix(input_data$counts)
    #design_matrix <- model.matrix(~age_cluster, metadata)
    #fit <- devil::fit_devil(peak_counts, design_matrix, verbose = F, size_factors = T)
    #clusters <- as.numeric(as.factor(metadata$patient))
    #res <- devil::test_de(fit, contrast = c(0,1), clusters = clusters, max_lfc = Inf)

  #} else if (method == "glmGamPoi") {
    #metadata <- input_data$metadata
    #peak_counts <- as.matrix(input_data$counts)
    #design_matrix <- model.matrix(~age_cluster, metadata)
    #fit <- glmGamPoi::glm_gp(peak_counts, design_matrix, size_factors = F, verbose = T)
    #res <- glmGamPoi::test_de(fit, contrast = c(0,1))
    #res <- res %>% select(name, pval, adj_pval, lfc)

  #} else if (method == 'nebula') {
    #metadata <- input_data$metadata
    #peak_counts <- as.matrix(input_data$counts)
    #design_matrix <- model.matrix(~age_cluster, metadata)
    #sf <- devil:::calculate_sf(counts)
    #metadata$patient_id <- as.factor(metadata$patient_id)
    #data_g = group_cell(count=peak_counts,id=clusters,pred=design_matrix)
    #fit <- nebula::nebula(data_g$count,id = data_g$id, pred = data_g$pred, ncore = 1)
    #res <- dplyr::tibble(
      #name = fit$summary$gene,
      #pval = fit$summary$p_groupTRUE,
      #adj_pval = p.adjust(fit$summary$p_groupTRUE, "BH"),
      #lfc=fit$summary$logFC_groupTRUE)
  #}
  #res
#}
