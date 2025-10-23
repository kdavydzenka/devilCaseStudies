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


perform_analysis_rna <- function(input_data, 
                                 method = "devil", 
                                 design_type = "age_only",
                                 cell_type_of_interest = NULL) {

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
    "age_plus_celltype" = model.matrix(~ age_cluster + cell_type, metadata),
    "interaction" = model.matrix(~ age_cluster * cell_type, metadata),
    stop("design_type must be one of: 'age_only', 'age_plus_celltype', 'interaction'")
  )

  ## Make contrast helper 
  make_contrast <- function(design, from, to = NULL) {
    c <- rep(0, ncol(design)); names(c) <- colnames(design)
    if (!is.null(to)) {
      stopifnot(from %in% names(c), to %in% names(c))
      c[from] <- 1; c[to] <- -1
    } else {
      stopifnot(from %in% names(c))
      c[from] <- 1
    }
    as.numeric(c)
  }

  age_cols <- grep("^age_cluster1", colnames(design_matrix), value = TRUE)[1]

  if (length(age_cols) == 0)
    stop("No age_cluster columns found — check metadata$age_cluster levels.")
  
  ## Fit and test 
  if (method == 'devil') {
    fit <- devil::fit_devil(
      input_matrix = counts,
      design_matrix = design_matrix,
      overdispersion = TRUE,
      offset = 1e-6,
      size_factors = "normed_sum",
      parallel.cores = 1,
      verbose = TRUE
    )

    clusters <- as.numeric(as.factor(metadata$sample))

    # iDefine contrasts depending on design_type
    if (design_type == "age_only" || design_type == "age_plus_celltype") {
      contrast <- make_contrast(design_matrix, from = age_cols[1])
    } else if (design_type == "interaction") {
      if (is.null(cell_type_of_interest)) {
        stop("Please specify cell_type_of_interest for interaction design")
      }
      # Find correct column name for the interaction term
      interaction_term <- paste0(age_cols, ":cell_type", cell_type_of_interest)
      if (!(interaction_term %in% colnames(design_matrix))) {
        stop(paste("Interaction term", interaction_term, "not found in design matrix"))
      }
      # Contrast for age effect within that cell type
      contrast <- make_contrast(design_matrix, from = "age_cluster1", to = interaction_term)
    }

    res <- devil::test_de(
      fit,
      contrast = contrast,
      clusters = clusters,
      max_lfc = 10
    )

  } else if (method == "glmGamPoi") {

    fit <- glmGamPoi::glm_gp(counts, design_matrix, size_factors = TRUE, verbose = TRUE)

    if (design_type == "age_only" || design_type == "age_plus_celltype") {
      contrast <- make_contrast(design_matrix, from = age_cols[1])
    } else if (design_type == "interaction") {
      if (is.null(cell_type_of_interest)) stop("Please specify cell_type_of_interest for interaction design")
      interaction_term <- paste0(age_cols, ":cell_type", cell_type_of_interest)
      contrast <- make_contrast(design_matrix, from = "age_cluster1", to = interaction_term)
    }

    res <- glmGamPoi::test_de(fit, contrast = contrast)
    res <- res %>% dplyr::select(name, pval, adj_pval, lfc)

  } else if (method == 'nebula') {

    metadata$patient <- as.numeric(as.factor(metadata$sample))
    sf <- devil:::calculate_sf(counts)
    #data_g = group_cell(count = counts, id = metadata$patient, pred = design_matrix, offset = sf)
    #print(str(data_g))
    #fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, offset = data_g$offset)
    fit <- nebula::nebula(counts, id = metadata$patient, pred = design_matrix, offset = sf)
    res <- fit$summary
    #res <- dplyr::tibble(
      #name = fit$summary$gene,
      #pval = fit$summary$p_age_cluster1,
      #adj_pval = p.adjust(fit$summary$p_age_cluster1, "BH"),
      #lfc = fit$summary$logFC_age_cluster1
    #)
  }

  res
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
