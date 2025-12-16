# Define GO clusters

GO_CLUSTERS <- list(
  `Immune Response` = c(
    "defense response",
    "immune system process",
    "myeloid leukocyte migration",
    "leukocyte migration",
    "cell surface toll-like receptor signaling pathway",
    "T cell mediated immunity",
    "positive regulation of cytokine production",
    "cellular response to cytokine stimulus",
    "lymphocyte activation",
    "innate immune response",
    "cell killing"
  ),
  
  `Gene Expression Regulation` = c(
    "positive regulation of gene expression",
    "positive regulation of transcription by RNA polymerase II",
    "ncRNA processing"
  ),
  
  `Cellular Adhesion and Communication` = c(
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "cell-cell adhesion via plasma-membrane adhesion molecules",
    "regulation of substrate adhesion-dependent cell spreading",
    "cell-cell signaling",
    "leukocyte cell-cell adhesion",
    "regulation of cell communication"
  ),
  
  `Cell Division, Development & Organelle Biogenesis ` = c(
   "epithelial cell proliferation",
   "positive regulation of cell differentiation",
   "cellular component assembly involved in morphogenesis",
   "non-membrane-bounded organelle assembly", 
   "organelle disassembly",
   "regulation of cell cycle process"
  ),
  
  `Muscle Function, Development & Cytoskeleton Dynamics` = c(
    "muscle contraction",
    "actin filament-based movement",
    "actin filament-based process",
    "actin-mediated cell contraction",
    "myofibril assembly",
    "muscle system process",
    "skeletal muscle cell differentiation",
    "striated muscle cell development",
    "muscle cell development"
  ),
  
  `Cellular Transport` = c(
    "vesicle-mediated transport",
    "regulation of protein secretion",
    "endocytosis"
  ),
  
  `Cell Signaling` = c(
    "regulation of signaling"
  ),
  
  `Metabolic Processes` = c(
    "aromatic compound catabolic process",
    "cellular nitrogen compound catabolic process",
    "proteolysis",
    "negative regulation of cellular metabolic process",
    "nucleic acid metabolic process",
    "protein metabolic process",
    "negative regulation of metabolic process",
    "RNA metabolic process",
    "nucleobase-containing compound catabolic process"
  ),
  
  `Stress Response` = c(
    "response to gamma radiation",
    "cellular response to hypoxia",
    "cellular response to chemical stimulus"
  )
)

#REACTOME_CLUSTERS <- list(
  #`Immune Response` = c(
    #"Toll-like Receptor Cascades",
    #"Innate Immune System"
  #),
  
  #`Gene Expression Regulation` = c(
    #"Transcriptional regulation by RUNX2",
    #"NGF-stimulated transcription",
    #"Estrogen-dependent gene expression",
    #"Generic Transcription Pathway",
    #"RNA Polymerase II Transcription",
    #"tRNA processing"
  #),

  #`Muscle and Movement Processes` = c(
    #"Striated Muscle Contraction",
    #"Muscle contraction",
    #"Smooth Muscle Contraction"
  #),
  
  #`Cell Death & Cellular Maintenance` = c(
    #"Programmed Cell Death",
    #"DNA Repair"
  #),
  
  #`Signaling and Regulation` = c(
    #"Signaling by MET",
    #"Interleukin-1 family signaling",
    #"Negative regulation of the PI3K/AKT network",
    #"PIP3 activates AKT signaling",
    #"Negative regulation of MAPK pathway",
    #"PTEN Regulation",
    #"ER-Phagosome pathway",
    #"RAS processing",
    #"Cytokine Signaling in Immune system",
    #"Signaling by Interleukins",
    #"Signal Transduction"
  #),
  
  #`Metabolic processes` = c(
    #"Metabolism of RNA",
    #"Metabolism of lipids",
    #"Metabolism of steroids",
    #"Iron uptake and transport"
  #),
  
  #`Stress Response` = c(
    #"Senescence-Associated Secretory Phenotype (SASP)",
    #"Oxidative Stress Induced Senescence"
  #)
#)

#remove_redundant_terms <- function(data, enrichment_col = "enrichmentScore", core_col = "core_enrichment", desc_col = "Description", threshold = 0.5) {
  #filtered_data <- data %>%
    #mutate(genes = strsplit(!!sym(core_col), "/"))
  #n_terms <- nrow(filtered_data)
  #overlap_matrix <- matrix(0, nrow = n_terms, ncol = n_terms,
                           #dimnames = list(filtered_data[[desc_col]], filtered_data[[desc_col]]))
  #for (i in 1:n_terms) {
    #for (j in i:n_terms) {
      #shared_genes <- length(intersect(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      #total_genes <- length(union(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      #jaccard_index <- shared_genes / total_genes

      # Fill overlap matrix with Jaccard index
      #overlap_matrix[i, j] <- jaccard_index
      #overlap_matrix[j, i] <- jaccard_index
    #}
  #}
  #redundant_terms <- c()
  #for (i in 1:(n_terms - 1)) {
    #for (j in (i + 1):n_terms) {
      #if (overlap_matrix[i, j] > threshold) {
        #redundant_terms <- c(redundant_terms, filtered_data[[desc_col]][j])
      #}
    #}
  #}
  #non_redundant_data <- filtered_data %>%
    #filter(!(!!sym(desc_col) %in% redundant_terms))

  #redundant_data <- filtered_data %>%
    #filter(!!sym(desc_col) %in% redundant_terms)

  #return(list(
    #non_redundant_data = non_redundant_data,
    #redundant_data = redundant_data
  #))
#}


get_simplified_GOterms = function(GO_devil, GO_glm, GO_nebula, by=.05) {
  sdevil = lapply(seq(0, 1, by =by), function(c) {
    gseGO_devil_s = clusterProfiler::simplify(GO_devil, cutoff=c)
    dplyr::tibble(model = "devil", n_simplified = nrow(gseGO_devil_s@result), n_total = nrow(GO_devil@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  sglm = lapply(seq(0, 1, by =by), function(c) {
    gseGO_glm_s = clusterProfiler::simplify(GO_glm, cutoff=c)
    dplyr::tibble(model = "glmGamPoi", n_simplified = nrow(gseGO_glm_s@result), n_total = nrow(GO_glm@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  snebula = lapply(seq(0, 1, by = .05), function(c) {
    gseGO_nebula_s = clusterProfiler::simplify(GO_nebula, cutoff=c)
    dplyr::tibble(model = "NEBULA", n_simplified = nrow(gseGO_nebula_s@result), n_total = nrow(GO_nebula@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  dplyr::bind_rows(sdevil, sglm, snebula)
}


plot_dotplot_GO = function(devil_res, glm_res, nebula_res, top_K = 15) {

  devil_res <- devil_res %>% dplyr::mutate(method = "devil", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))
  glm_res <- glm_res %>% dplyr::mutate(method = "glmGamPoi", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))
  nebula_res <- nebula_res %>% dplyr::mutate(method = "nebula", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))
  
  replacement_map <- c(
    "innate immune response" = "immune system process",
    "actin filament-based process" = "actin-mediated cell contraction",
    "muscle cell development" = "striated muscle cell development",
    "striated striated muscle cell development" = "striated muscle cell development",
    "homophilic cell adhesion via plasma membrane adhesion molecules" = "homophilic cell adhesion"
  )
  
  combined_data <- bind_rows(devil_res, glm_res, nebula_res) %>%
    dplyr::filter(!(Description %in% c("organic cyclic compound catabolic process", "defense response to other organism",
                                       "biological process involved in interspecies interaction between organisms",
                                       "regulation of response to external stimulus", "actin filament-based movement",
                                       "biological process involved in interspecies interaction between organisms",
                                       "negative regulation of nitrogen compound metabolic process",
                                       "positive regulation of response to stimulus", "positive regulation of multicellular organismal process",
                                       "macromolecule modification", "anatomical structure morphogenesis", "response to organic substance",
                                       "positive regulation of gene expression", "regulation of response to external stimulus", 
                                       "positive regulation of cellular biosynthetic process", "response to chemical",
                                       "cell population proliferation", "system development", "regulation of biological process", 
                                       "multicellular organismal process", "cellular component organization or biogenesis",
                                       "cell-cell adhesion via plasma-membrane adhesion molecules",
                                       "response to gamma radiation", "cognition", "inner ear development",
                                       "T cell activation"
                                       ))) 

  combined_data$GO_cluster = lapply(combined_data$Description, function(go_term) {
    for (cluster_name in names(GO_CLUSTERS)) {
      if (go_term %in% GO_CLUSTERS[[cluster_name]]) {
        return(cluster_name)
      }
    }
    print(go_term)
    return("GO term not found in any cluster")
  }) %>% unlist()

  combined_data$GO_cluster <- str_wrap(combined_data$GO_cluster, width = 20)
  
  # Select top 10 terms per method and DE_type based on enrichmentScore
  filtered_data <- combined_data %>%
    group_by(method) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_K) %>%
    ungroup()
  
  filtered_data <- filtered_data %>% 
    dplyr::mutate(Description = str_replace_all(Description, replacement_map))

  # Enrichment plot

  plot_GO = filtered_data %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = filtered_data %>%
          arrange(factor(method, levels = c("devil", "glmGamPoi", "nebula")), enrichmentScore) %>%
          pull(Description) %>%
          unique()
      ),
      DE_type = factor(DE_type, levels = c("Up-regulated", "Down-regulated")) # Ensure DE_type is ordered
    ) %>%
    ggplot(aes(x = method, y = Description, size = setSize, color = p.adjust)) +
    geom_point() +
    facet_grid(GO_cluster~DE_type, space = "free", scales = "free") +
    scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p-value") +
    theme_bw() +
    labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Count") +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  plot_GO
}


#enrichmentGO <- function(rna_deg_data) {
  #rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] = min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0])
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * rna_deg_data$lfc
  #rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  #genes <- rna_deg_data$RankMetric
  #names(genes) <- rna_deg_data$geneID

  #gseGO <- clusterProfiler::gseGO(
    #genes,
    #ont = "BP",  
    #OrgDb = org.Hs.eg.db,
    #minGSSize = 10,
    #maxGSSize = 350,
    #keyType = "SYMBOL",
    #pvalueCutoff = 0.05,
    #pAdjustMethod = "BH",
    #verbose = TRUE,
    #eps = 0,
    #nPermSimple = 10000
    #by ="fgsea",
    #nPerm = 10000,
    #seed = 123
  #)
  #return(gseGO)
  #return(gseGO@result %>% as.data.frame())
#}

runGO <- function(df,
                  padj_col,
                  lfc_col,
                  method = c("GSE", "ENRICH"),
                  pval_cut = 0.05,
                  lfc_cut  = 0,
                  ont = "BP") {
  
  method   <- match.arg(method)
  padj_sym <- rlang::sym(padj_col)
  lfc_sym  <- rlang::sym(lfc_col)
  
  # --- 1. Clean input & standardize column names ---
  rna_deg_data <- df %>%
    dplyr::filter(
      !is.na(!!padj_sym),
      !is.na(!!lfc_sym),
      !!padj_sym > 0
    ) %>%
    dplyr::rename(
      adj_pval = !!padj_sym,
      lfc      = !!lfc_sym,
      geneID   = name
    )
  
  # handle zeros in adj_pval
  min_p <- min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0], na.rm = TRUE)
  rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] <- min_p
  
  # --- 2. Branch by method ---
  if (method == "GSE") {
    # rank metric for GSEA
    rna_deg_data <- rna_deg_data %>%
      dplyr::mutate(RankMetric = -log10(adj_pval) * sign(lfc)) %>%
      #dplyr::mutate(RankMetric = lfc) %>%
      #dplyr::mutate(RankMetric = -log(pval) * sign(lfc)) %>%
      dplyr::arrange(dplyr::desc(RankMetric))
    
    genes <- rna_deg_data$RankMetric
    names(genes) <- rna_deg_data$geneID
    
    res <- clusterProfiler::gseGO(
      geneList      = genes,
      ont           = ont,
      OrgDb         = org.Hs.eg.db,
      minGSSize     = 10,
      maxGSSize     = 350,
      keyType       = "SYMBOL",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = TRUE, 
      eps           = 0,
      nPermSimple   = 50000
    )
    
  } else if (method == "ENRICH") {
    # significant genes for over-representation
    up_sig_genes <- rna_deg_data %>%
      dplyr::filter(adj_pval < pval_cut, lfc >= lfc_cut) %>%
      dplyr::pull(geneID)
    
    down_sig_genes <- rna_deg_data %>%
      dplyr::filter(adj_pval < pval_cut, lfc <= -lfc_cut) %>%
      dplyr::pull(geneID)
    
    universe_genes <- rna_deg_data$geneID
    
    up_res <- clusterProfiler::enrichGO(
      gene          = up_sig_genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 10,
      maxGSSize     = 350
    )
    
    down_res <- clusterProfiler::enrichGO(
      gene          = down_sig_genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 10,
      maxGSSize     = 350
    )
    
    res = list(up_res = up_res, down_res = down_res, 
               up_res_df = up_res@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Up-regulated"),
               down_res_df = down_res@result %>% dplyr::filter(p.adjust <= pval_cut) %>% dplyr::mutate(regulation = "Down-regulated"))
  }
  
  return(res)
}

enrichmentGO <- function(df, padj_col, lfc_col) {
  
  # drop NAs and zeros safely
  rna_deg_data <- df %>%
    filter(!is.na(.data[[padj_col]]),
           !is.na(.data[[lfc_col]]),
           .data[[padj_col]] > 0) %>%
    dplyr::rename(adj_pval = !!padj_col,
                  lfc = !!lfc_col,
                  geneID = name)
  
  # handle zeros
  min_p <- min(rna_deg_data$adj_pval[rna_deg_data$adj_pval != 0], na.rm = TRUE)
  rna_deg_data$adj_pval[rna_deg_data$adj_pval == 0] <- min_p
  
  # rank metric
  rna_deg_data <- rna_deg_data %>%
    mutate(RankMetric = -log10(adj_pval) * sign(lfc)) %>%
    arrange(desc(RankMetric))
  
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID
  
  # run enrichment
  gseGO <- clusterProfiler::gseGO(
    geneList = genes,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    maxGSSize = 350,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 0,
    nPermSimple = 10000
  )
  
  return(gseGO)
}


#enrichmentReactomePA <- function(rna_deg_data) {
  #rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  #rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  #genes <- rna_deg_data$RankMetric
  #names(genes) <- rna_deg_data$geneID
  
  #entrez_ids <- mapIds(org.Hs.eg.db, keys = names(genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  #entrez_ids <- as.data.frame(entrez_ids)
  #names(genes) <- entrez_ids$entrez_ids
  
  #gseReactome <- ReactomePA::gsePathway(genes, 
                                        #pvalueCutoff = 0.2,
                                        #pAdjustMethod = "BH", 
                                        #verbose = FALSE)
  
  #return(gseReactome@result %>% as.data.frame())
#}


# Function to compute the area under the curve (AUC) using the trapezoidal rule

auc <- function(x, y) {
  # Check if the lengths of x and y are the same
  if (length(x) != length(y)) {
    stop("Vectors x and y must have the same length")
  }

  # Sort the points based on x values (in case they are not ordered)
  ord <- order(x)
  x_ordered <- x[ord]
  y_ordered <- y[ord]

  # Compute the area using the trapezoidal rule
  n <- length(x_ordered)
  area <- sum((y_ordered[-1] + y_ordered[-n]) * diff(x_ordered)) / 2

  return(area)
}

get_tissue_specific_res = function(gse_result, method) {
  # Extract gene sets (significant ones)
  if (nrow(gse_result@result) == 0) return(dplyr::tibble())
  significant_terms <- subset(gse_result@result, qvalue < 0.05)
  if (nrow(gse_result@result) == 0) return(dplyr::tibble())
  
  # Get the genes for each GO term (pathway)
  # Map ENTREZ IDs to SYMBOLs if needed
  all_gene_sets <- lapply(significant_terms$ID, function(go_id) {
    genes <- DOSE::geneInCategory(gse_result)[[go_id]]
    genes
  })
  names(all_gene_sets) <- significant_terms$Description
  
  # Run TissueEnrich on each gene set
  results <- lapply(seq_along(all_gene_sets), function(i) {
    gene_set <- all_gene_sets[[i]]
    gs <- GeneSet(geneIds=gene_set,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
    if (length(gene_set) >= 5) { # Minimum size
      TissueEnrich::teEnrichment(gs)
    } else {
      NULL
    }
  })
  names(results) <- names(all_gene_sets)
  
  # Extract enrichment score (p-value or fold change) for the tissue of interest
  res_df = lapply(1:length(results), function(i) {
    path_name = names(results)[i]
    te_result = results[[i]]
    if (is.null(te_result)) {
      print(i)
      return(NULL)
    }
    seEnrichmentOutput<-te_result[[1]]
    enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
    enrichmentOutput$Tissue<-row.names(enrichmentOutput)
    enrichmentOutput %>% dplyr::mutate(path_name = path_name)
  }) %>% do.call("bind_rows", .)
  rownames(res_df) = NULL
  res_df %>% dplyr::mutate(method = method)
}


bayes_de_prob <- function(df, lfc_cut = 1.0, alpha = 0.05, k = 8) {
  # avoid zeros in adj_pval
  if (any(df$adj_pval == 0, na.rm = TRUE)) {
    min_nonzero <- min(df$adj_pval[df$adj_pval > 0], na.rm = TRUE)
    df$adj_pval[df$adj_pval == 0] <- min_nonzero
  }
  
  df %>%
    dplyr::mutate(
      # p-value based DE "probability"
      p_de_p = 1 / (1 + adj_pval / alpha),
      
      # LFC-based DE "probability": 0.5 at |lfc| = lfc_cut
      p_de_lfc = 1 / (1 + exp(-k * (abs(lfc) - lfc_cut))),
      
      # combined DE probability (you can change comb rule if you like)
      p_de = p_de_p * p_de_lfc
    )
}
