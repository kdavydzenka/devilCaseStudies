
DATASET_NAMES <- c("BaronPancreasData", 'ZhaoImmuneLiverDataBlood', 'DarmanisBrainData', 'ZhaoImmuneLiverDataLiver',
                   "BigBloodData", "BigLiverData", "pbmc", "liver")

SEED <- 12345

read_data <- function(dataset_name, data_path=NULL) {
  if (dataset_name == "BaronPancreasData") {
    data <- scRNAseq::BaronPancreasData(which = "human")
    metadata <- data@colData
    metadata$cell_type <- metadata$label
    tissue <- "pancreas"
    counts <- data@assays@data[[1]]
  } else if (dataset_name == "DarmanisBrainData") {
    data <- scRNAseq::DarmanisBrainData()
    metadata <- data@colData
    metadata$cell_type <- metadata$cell.type
    metadata$donor <- metadata$experiment_sample_name
    tissue <- "brain"
    counts <- data@assays@data[[1]]
  } else if (dataset_name == 'ZhaoImmuneLiverDataBlood') {
    data <- scRNAseq::ZhaoImmuneLiverData()

    filter_cell <- grepl('blood', data$sample)
    metadata <- dplyr::tibble(
      cell_type = data$fine[filter_cell],
      donor = data$sample[filter_cell]
    )
    counts <- data@assays@data[[1]]
    counts <- counts[,filter_cell]
    tissue <- "blood"
  } else if (dataset_name == 'ZhaoImmuneLiverDataLiver') {
    data <- scRNAseq::ZhaoImmuneLiverData()

    filter_cell <- grepl('liver', data$sample)
    metadata <- dplyr::tibble(
      cell_type = data$fine[filter_cell],
      donor = data$sample[filter_cell]
    )
    counts <- data@assays@data[[1]]
    counts <- counts[,filter_cell]
    tissue <- "liver"
  } else if (dataset_name == 'BigBloodData') {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    metadata <- dplyr::rename(metadata,c("cell_type" = "celltype.l2", "cell_type2" = "cell_type", "donor" = "donor_id"))
    tissue = "blood"
  } else if (dataset_name == 'BigLiverData') {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    meta_features <- seurat_data@assays[["RNA"]]@meta.features
    rownames(counts) <- meta_features$feature_name
    metadata <- dplyr::rename(metadata,c("donor" = "donor_id"))
    tissue = "liver"
  } else if (dataset_name == "pbmc") {
    seurat_data <- readRDS("datasets/pbmc.rds")
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    metadata <- dplyr::rename(metadata %>% dplyr::select(cell_type, donor_id),c("cell_type" = "cell_type", "donor" = "donor_id"))
    tissue = "blood"
  } else if (dataset_name == "liver") {
    seurat_data <- readRDS("datasets/liver.rds")
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    metadata <- dplyr::rename(metadata %>% dplyr::select(donor_id, BroadCellType), c("cell_type" = "BroadCellType", "donor" = "donor_id"))
    tissue = "liver"
  } else {
    stop("Dataset name not recognized")
  }

  return(list(counts=counts, metadata=metadata, tissue=tissue))
}

prep_seurat_object <- function(input_data, NPC=50, cluster_res=1) {
  counts <- input_data$counts
  metadata <- input_data$metadata
  metadata$cell_type %>% unique()
  metadata$donor %>% unique()

  # Filter cells
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)

  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  feat_mad_filter <- total_features > 5 * mad(total_features)

  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .1

  # ribosomal_genes <- grepl("^RPS", rownames(counts)) | grepl("^RPL", rownames(counts))
  # ribosomal_prop <- colSums(counts[ribosomal_genes, ]) / colSums(counts)
  # rib_prop_filter <- ribosomal_prop < .1

  #cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter | rib_prop_filter | mit_prop_filter
  cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter |  mit_prop_filter

  counts <- counts[, !cell_outliers_filter]
  metadata <- metadata[!cell_outliers_filter, ]
  rownames(metadata) <- colnames(counts)

  # Filter genes
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]

  seurat_obj = SeuratObject::CreateSeuratObject(counts = counts, meta.data = as.tibble(metadata))
  seurat_obj <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = NPC) %>%
    FindNeighbors(dims = 1:NPC) %>%
    FindClusters(resolution = cluster_res) %>%
    RunUMAP(dims=1:NPC)

  new_clusters <- as.numeric(seurat_obj$seurat_clusters)
  new_clusters <- factor(new_clusters, levels = sort(unique(new_clusters)))
  seurat_obj$seurat_clusters <- new_clusters

  seurat_obj
}

perform_analysis <- function(seurat_obj, method = "devil") {
  if (!(method %in% c('devil', "glmGamPoi", 'nebula'))) {stop('method not recognized')}

  counts <- as.matrix(seurat_obj@assays$RNA$counts)

  whole_res <- dplyr::tibble()
  for (c in unique(seurat_obj$seurat_clusters)) {

    print(c)

    idx_cluster <- which(seurat_obj$seurat_clusters == c)
    idx_others <- which(!(seurat_obj$seurat_clusters == c))

    if (length(idx_others) > length(idx_cluster)) {
      set.seed(SEED)
      #idx_others <- sample(idx_others, as.integer(length(idx_others) * .05), replace = F)
      idx_others <- sample(idx_others, length(idx_cluster), replace = F)
    }

    cell_idx <- c(idx_cluster, idx_others)
    clusters <- as.numeric(as.factor(seurat_obj$donor))
    design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_obj$seurat_clusters == c))

    # First filter
    dm <- design_matrix[cell_idx,]
    cc <- counts[,cell_idx]
    clusters <- clusters[cell_idx]

    # Second filter
    cell_idx <- which((colSums(cc) > 20) == TRUE)
    dm <- dm[cell_idx,]
    cc <- cc[,cell_idx]
    clusters <- clusters[cell_idx]

    # Third filter
    gene_idx = which((rowSums(cc) > 20) == TRUE)
    cc <- cc[gene_idx,]

    rownames(dm) <- colnames(cc)

    if (method == 'devil') {
      s <- Sys.time()
      fit <- devil::fit_devil(cc, dm, size_factors = "normed_sum", overdispersion = "old", 
                              init_overdispersion = 100, offset = 1e-6, verbose = TRUE, tolerance = 1e-3, 
                              max_iter = 200, parallel.cores = 1)
      e <- Sys.time()
      res <- devil::test_de(fit, contrast = c(0,1), clusters = clusters, max_lfc = 50) %>% dplyr::mutate(cluster = c)
      #res <- devil::test_de(fit, contrast = c(0,1), clusters = 1:length(idxs), max_lfc = Inf) %>% dplyr::mutate(cluster = c)
    } else if (method == "glmGamPoi") {
      s <- Sys.time()
      fit <- glmGamPoi::glm_gp(cc, dm, size_factors = "normed_sum", verbose = T)
      e <- Sys.time()
      #fit <- glmGamPoi::glm_gp(cc, dm, size_factors = FALSE, verbose = T)
      res <- glmGamPoi::test_de(fit, contrast = c(0,1))
      res <- res %>% dplyr::as_tibble() %>% dplyr::select(name, pval, adj_pval, lfc) %>% dplyr::mutate(cluster = c)
    } else if (method == "nebula") {
      s <- Sys.time()
      sf <- devil:::calculate_sf(Y = cc, method = "normed_sum")
      data_g = nebula::group_cell(count=cc,id=clusters,pred=dm)
      fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, ncore = 1, mincp = 0, cpc = 0, offset = sf)
      e <- Sys.time()
      #fit <- nebula::nebula(data_g$count, id = data_g$id, pred = data_g$pred, ncore = 1, mincp = 0, cpc = 0)
      res <- dplyr::tibble(
        name = fit$summary$gene,
        pval = fit$summary$p_groupTRUE,
        adj_pval = p.adjust(fit$summary$p_groupTRUE, "BH"),
        lfc=fit$summary$logFC_groupTRUE
      ) %>% dplyr::mutate(cluster = c)
    } else {
      stop("method not recognized")
    }

    res <- res %>% dplyr::mutate(delta_time = e - s)
    whole_res <- dplyr::bind_rows(whole_res, res)
  }
  whole_res
}

prepScMayoInput <- function(de_res, count_matrix, seurat_clusters, n_markers=50, lfc_cut=1, pval_cut=1e-10) {
  if (!(all(colnames(de_res) == c("name", "pval", "adj_pval", "lfc", "cluster")))) {stop('input de_res have wrong column names')}

  cluster_values <- de_res$cluster %>% unique()
  remove_genes <- grepl("^ENS", de_res$name)
  de_res <- de_res[!remove_genes, ]

  scMayoInput <- dplyr::tibble()
  for (c in cluster_values) {
    top_res <- de_res %>%
      dplyr::filter(cluster == c, lfc > lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::slice(1:n_markers)

    if (nrow(top_res) > 0) {
      names <- rownames(count_matrix)
      if (all(grepl("ENSG", names))) {
        ensembl_names <- names
        symbol_names <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=ensembl_names,column='SYMBOL', keytype="ENSEMBL", multiVals="first")
      } else {
        symbol_names <- names
        ensembl_names <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=symbol_names,column="ENSEMBL", keytype='SYMBOL', multiVals="first")
      }

      gene_percentage <- dplyr::tibble()
      for (gene_name in top_res$name) {
        if (length(which(gene_name == symbol_names)) > 0) {
          f <- gene_name == symbol_names
          f[is.na(f)] <- F
        } else {
          f <- gene_name == ensembl_names
          f[is.na(f)] <- F
        }

        g_counts <- count_matrix[f, seurat_clusters == c]
        non_g_counts <- count_matrix[f, !(seurat_clusters == c)]

        pct.1 <- sum(g_counts > 0) / length(g_counts)
        pct.2 <- sum(non_g_counts > 0) / length(non_g_counts)

        gene_percentage <- dplyr::bind_rows(gene_percentage, dplyr::tibble(name = gene_name, pct.1=pct.1, pct.2=pct.2))
      }

      top_res <- top_res %>% dplyr::left_join(gene_percentage, by="name")
      scMayoInput <- dplyr::bind_rows(scMayoInput, top_res)
    }
  }

  scMayoInput <- scMayoInput %>%
    dplyr::rename(gene=name, p_val_adj=adj_pval, p_val=pval, pct.1=pct.1, pct.2=pct.2, cluster=cluster, avg_log2FC=lfc)

  #rownames(scMayoInput) <- scMayoInput$gene
  scMayoInput
}

plot_scMayoOutput <- function(scMayoObj) {
  pred_res <- lapply(1:nrow(scMayoObj$res), function(i) {
    r <- scMayoObj$res[i,]
    s <- strsplit(r$ES.norm, split = ";") %>% unlist() %>% as.numeric()
    c <- strsplit(r$celltype, split = ";") %>% unlist()
    c <- lapply(c, function(x) {stringr::str_replace_all(x, "cell", "")}) %>% unlist()
    c <- lapply(c, function(x) {stringr::str_replace_all(x, " ", "")}) %>% unlist()
    s
    dplyr::tibble(cell_type_pred = c, score = s, cluster = r$cluster)
  }) %>% do.call('bind_rows', .)

  ground_truth <- computeGroundTruth(seurat_obj)

  ground_truth %>%
    dplyr::left_join(pred_res, by="cluster") %>%
    dplyr::arrange(cluster) %>%
    dplyr::mutate(ground_truth = as.factor(paste0(true_cell_type, ' (', cluster, ')'))) %>%
    #ggplot(mapping = aes(x=cell_type_pred, y=ground_truth, fill=score)) +
    ggplot(mapping = aes(y=cell_type_pred, x=ground_truth, fill=score)) +
    geom_tile() +
    scale_fill_continuous(type = "viridis") +
    theme_bw() +
    #labs(y = "Annotated cell type", x = "Assigned cell type", fill="Score") +
    labs(y = "Prediction", x = "Ground truth", fill="Score")
}

computeGroundTruth <- function(seurat_obj) {
  lapply(unique(seurat_obj$seurat_clusters), function(c) {
    table_c <- seurat_obj@meta.data %>%
      dplyr::filter(seurat_clusters == c) %>%
      dplyr::pull(cell_type) %>%
      table()
    table_c <- table_c / sum(table_c)
    max_idx <- which.max(table_c)
    dplyr::tibble(cluster = c, true_cell_type = names(table_c)[max_idx], pct = table_c[max_idx])
  }) %>% do.call('bind_rows', .)
}


# cell_type_names_to_scMayo_names <- function(ct, tissue) {
#   # scMayo_names <- colnames(scMayoMap::scMayoMapDatabase)[grepl(tissue, colnames(scMayoMap::scMayoMapDatabase))]
#   # scMayo_names <- str_replace_all(scMayo_names, tissue, "")
#   # scMayo_names <- str_replace_all(scMayo_names, ":", "")
#   # scMayo_names <- str_replace_all(scMayo_names, " cell", "")
#
#   if (tissue == "pancreas") {
#     conversion <- list(
#       "acinar" = "Acinar",
#       'beta' = 'Beta',
#       'delta'='Delta',
#       'activated_stellate'="Stellate",
#       'ductal'="Ductal",
#       'alpha'="Alpha",
#       "epsilon"="Epsilon",
#       "gamma"="Gamma",
#       "endothelial"="Endothelial",
#       'quiescent_stellate'="Stellate",
#       "macrophage"="Macrophage",
#       "schwann"="Schwann",
#       "mast"="Mast",
#       "t_cell"="T"
#     )
#   } else if (tissue == "blood") {
#     conversion = list(
#       'CD14 Mono'="CD14 Monocyte",
#       'CD4 TCM'='CD4 Central Memory T',
#       'CD8 Naive'='CD8 Naive T',
#       'NK'='Natural killer',
#       'CD8 TEM'='CD8 Effector Memory T',
#       'CD16 Mono'="CD16 Monoctye",
#       'B intermediate'="Intermediate B",
#       'CD4 Naive'="CD4 Naive T",
#       'CD4 CTL'="CD4 Cytotoxic T",
#       'B naive'="Naive B",
#       'MAIT'="Mucosal-associated invariant T",
#       'gdT'='Gamma delta T',
#       'gamma-delta T cell'='Gamma delta T',
#       'CD8 TCM'="CD8 Central Memory T",
#       'dnT'="Double-negative T",
#       'B memory'="Memory B",
#       'Doublet'="",
#       'pDC'='Plasmacytoid dendritic',
#       'CD8 Proliferating'='CD8 Proliferating T',
#       'Treg'="Regulatory T",
#       'Plasmablast'="Plasmablast",
#       'CD4 TEM'="CD4 Effector Memory T",
#       'cDC2'="Dendritic",
#       'NK Proliferating'='Natural killer',
#       'ASDC'="Dendritic",
#       'HSPC'='Hematopoietic stem',
#       'Platelet'='Platelet',
#       'NK_CD56bright'='CD56-bright natural killer',
#       'CD4 Proliferating'='CD4 Proliferating T',
#       'Eryth'="Erythroid",
#       'cDC1'='Dendritic',
#       'ILC'='Innate lymphoid',
#       "natural killer cell" = "Natural killer",
#       'central memory CD4-positive, alpha-beta T cell'="CD4 Central Memory T",
#       'effector memory CD8-positive, alpha-beta T cell' = 'CD8 Effector Memory T',
#       'CD4-positive, alpha-beta cytotoxic T cell'="CD4 Central Memory T",
#       "classical monocyte"="Monocyte",
#       "naive B cell"="Naive B",
#       "platelet"="Platelet",
#       "non-classical monocyte"='CD16 Monocyte',
#       "memory B cell"="Memory B",
#       "erythrocyte"="Erythroid precursor",
#       "plasmacytoid dendritic cell, human"="Plasmacytoid dendritic",
#       "naive thymus-derived CD4-positive, alpha-beta T cell"='CD4 Naive T',
#       "mucosal invariant T cell"='Mucosal-associated invariant T',
#       "lymphocyte"="Lymphoid",
#       "regulatory T cell"="Regulatory T",
#       'transitional stage B cell'="Intermediate B",
#       "naive thymus-derived CD8-positive, alpha-beta T cell"='CD8 Naive T',
#       "CD16-negative, CD56-bright natural killer cell, human"="CD56-bright natural killer",
#       "plasmablast"="Plasmablast",
#       "neutrophil"="Neutrophil",
#       "conventional dendritic cell"="Dendritic",
#       "hematopoietic stem cell"="Hematopoietic stem",
#       "immature neutrophil"="Neutrophil"
#     )
#   } else if (tissue == "liver") {
#     conversion <- list(
#       'Erythroid'='Erythroid',
#       'monocyte'='Monocyte',
#       'conventional dendritic cell'='Dendritic',
#       'plasmacytoid dendritic cell'='Dendritic',
#       'macrophage'='Macrophage',
#       'plasma cell'='Plasma',
#       'T cell'='T',
#       'B-cell'='B',
#       'B cell'='B',
#       'Kupffer' = "Kupffer",
#       'natural killer cell'='Natural killer',
#       'cholangiocyte'='Cholangiocyte',
#       'Cholangiocyte'='Cholangiocyte',
#       'T/NK-cell' = "T",
#       'Hepatocyte'='Hepatocyte',
#       'hepatocyte'='Hepatocyte',
#       'basophil'='Basophil',
#       'neutrophil'='Neutrophil',
#       'endothelial cell'='Endothelial',
#       'Endothelial'='Endothelial',
#       'b-cell' = "B",
#       'Stellate' = "Stellate",
#       'fibroblast'='Myofibroblast'
#     )
#   } else {
#     stop("tissue name not recognised")
#   }
#
#   conv <- dplyr::tibble(start = names(conversion), end=conversion %>% unlist() %>% unname())
#   new_cell_names <- as.character(ct)
#
#   for (c in unique(ct)) {
#     print(c)
#     new_name <- conv %>% dplyr::filter(start == c) %>% pull(end)
#     new_cell_names[which(new_cell_names == c)] <- new_name
#   }
#   new_cell_names
# }

my_large_palette <- c(
  "steelblue4",
  "#D98880",
  "goldenrod",
  "indianred3",
  "mediumpurple",
  "brown",
  "plum",
  "tan",
  "darkseagreen",
  'cyan3',
  "lightsteelblue",
  "peru",
  "olivedrab",
  "palevioletred",
  "firebrick"
)
my_large_palette <- c(my_large_palette, my_large_palette, my_large_palette)

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "nebula" =  "steelblue",
  "NEBULA" =  "steelblue",
  "devil" = "#099668"
)
