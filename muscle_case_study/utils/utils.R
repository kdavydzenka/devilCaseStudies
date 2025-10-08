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


perform_analysis_rna <- function(input_data, method = "devil") {
  if (!(method %in% c('devil', "glmGamPoi", 'nebula'))) {stop('method not recognized')}
  
  if (method == 'devil') {
    metadata <- input_data$metadata
    counts <- as.matrix(input_data$counts)
    design_matrix <- model.matrix(~age_cluster, metadata)
    fit <- devil::fit_devil(counts, design_matrix, verbose = T, size_factors = T)
    clusters <- as.numeric(as.factor(metadata$sample))
    res <- devil::test_de(fit, contrast = c(0,1), clusters = clusters, max_lfc = Inf)
    
  } else if (method == "glmGamPoi") {
    metadata <- input_data$metadata
    counts <- as.matrix(input_data$counts)
    design_matrix <- model.matrix(~age_cluster, metadata)
    fit <- glmGamPoi::glm_gp(counts, design_matrix, size_factors = T, verbose = T)
    res <- glmGamPoi::test_de(fit, contrast = c(0,1))
    res <- res %>% select(name, pval, adj_pval, lfc)
    
  } else if (method == 'nebula') {
    metadata <- input_data$metadata
    counts <- as.matrix(input_data$counts)
    design_matrix <- model.matrix(~age_cluster, metadata)
    metadata$patient <- as.numeric(as.factor(metadata$sample))
    sf <- devil:::calculate_sf(counts)
    #data_g = group_cell(count=counts,id=metadata$orig.ident,pred=design_matrix)
    fit <- nebula::nebula(counts, id = metadata$patient, pred = design_matrix, offset = sf)
    res <- dplyr::tibble(
      name = fit$summary$gene,
      pval = fit$summary$p_age_cluster1,
      adj_pval = p.adjust(fit$summary$p_age_cluster1, "BH"),
      lfc = fit$summary$logFC_age_cluster1
    )
  }
  res
}

grange_annot <- function(input_data, data_path) {
  counts <- input_data$counts
  metadata <- input_data$metadata
  grange <- input_data$grange
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  grange_annot <- ChIPseeker::annotatePeak(
    peak = grange,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    level = "gene",
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                  "Downstream", "Intergenic"),
    annoDb = "org.Hs.eg.db",
    addFlankGeneInfo = FALSE,
    flankDistance = 5000,
    sameStrand = FALSE,
    ignoreOverlap = FALSE,
    ignoreUpstream = FALSE,
    ignoreDownstream = FALSE,
    overlap = "TSS",
    verbose = TRUE,
    columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
  ) %>% as.data.frame()
  tissue = "muscle"
  return(list(counts=counts, metadata=metadata, grange=grange_annot, tissue=tissue))
}

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

input_data 
