
custom_labels <- function(x) {
  sapply(x, function(xi) {
    if (xi == floor(xi)) {
      as.character(xi)       # No decimals for integers
    } else {
      format(round(xi, 2), nsmall = 2)  # 2 decimals otherwise
    }
  })
}

method_colors = c(
  #"glmGamPoi (Pb)" = "#A22E29",
  #"edgeR" = "#7D629E",
  "limma" = "#B96461",
  "limma (Pb)" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "glmGamPoi" = "#EAB578",
  "Nebula" =  'steelblue', #"#B0C4DE",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "Devil (base)" = "#099668",
  "Devil (mixed)" = "#099668",
  "Devil" = "#099668",
  "devil" = "#099668"
)

method_levels = c("devil", "NEBULA", "glmGamPoi", "limma")


all_null_plots <- function(author, is.pb, algos = c("glmGamPoi (Pb)", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)"),
                                                    ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL, only_tibble=FALSE) {
  if (is.pb) {
    head_foler = "nullpower/null_subject"
  } else {
    head_foler = "nullpower/null_cell"
  }

  fl <- list.files(head_foler, full.names = TRUE)[grepl(author, list.files(head_foler, full.names = TRUE))]

  if (is.null(n_samples_vec)) {n_samples_vec <- str_extract_all(fl, "(?<=n\\.)\\d+(?=\\.ngenes)") %>% unlist() %>% unique()}
  n_samples <- 4
  plots <- lapply(n_samples_vec, function(n_samples) {
    print(n_samples)
    fl_samples <- fl[grepl(paste0("n.", n_samples), fl)]

    if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.prob)") %>% unlist() %>% unique() }
    if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    iter.values <- str_extract_all(fl_samples, "(?<=iter.)\\d+(?=\\.cs)") %>% unlist() %>% unique()
    if (is.null(pde.values)) {
      pde.values <- lapply(fl_samples, function(l) {
        unlist(strsplit(unlist(strsplit(l, "probde."))[2], ".iter"))[1]
      }) %>% unlist() %>% unique()
    }


    dd <- lapply(ct.indexes, function(ct.index) {
      dd_inner <- lapply(pde.values, function(pde) {
        dd_inner_inner <- lapply(iter.values, function(i.iter) {
          n_genes <- as.integer(as.numeric(pde) * 1000)

          file_name <- paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".probde.",pde,".iter.",i.iter,".csv")
          if (file.exists(file_name)) {
            d <- read.delim(file_name, sep = ",")
            colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")

            mask <- colnames(d) %in% algos
            mask[1] <- TRUE
            d <- d[,mask]

            cols <- colnames(d)
            d <- lapply(2:ncol(d), function(c) {
              values = d[,c] %>% sort()
              x = seq(0,1,length = length(values))
              dplyr::tibble(x = x, observed_p_value = values, name = colnames(d)[c])
            }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
            return(d)
          }
        }) %>% do.call("bind_rows", .)
        dd_inner_inner
      }) %>% do.call('bind_rows', .)
    }) %>% do.call("bind_rows", .)

    colnames(dd)

    dd <- dd %>%
      #group_by(name, ct.index, n.genes) %>%
      group_by(name, n.genes) %>%
      dplyr::arrange(observed_p_value) %>%
      dplyr::mutate(x = row_number() / n())
    xy_lines <- dd %>%
      group_by(ct.index, n.genes) %>%
      dplyr::summarise(min = 0, max = 1) %>%
      pivot_longer(!c(ct.index, n.genes))

    if (only_tibble) {
      return(dd %>% dplyr::mutate(n.samples = n_samples))
    }

    p <- dd %>%
      #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
      ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
      #geom_line() +
      geom_point(size = .1) +
      geom_line(data = xy_lines, mapping = aes(x=value, y=value), col="black", linetype="dashed") +
      #geom_abline(slope = 1, intercept = 0, col = "black", linetype="dashed") +
      #ggtitle(paste0("Cell type ", ct.index, " - ", n_samples, " patients"), subtitle = paste0("Patient hierarchy ", is.pb)) +
      theme_bw() +
      labs(x = "Uniform quantiles", y="Observed p-value", col="Algorithm") +
      #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
      scale_color_manual(values = method_colors) +
      #facet_wrap(~paste0(n.genes, " genes")) +
      ggh4x::facet_nested(~paste0(n_samples, " patient")+"N genes"+n.genes) +
      theme(legend.position = "bottom") +
      scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n=3))

    p #+ ggtitle(paste("N patients = ", n_samples))
  })

  plots
}


all_pow_plots <- function(author, is.pb, algos = c("glmGamPoi (Pb)", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)"),
                          ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL, only_tibble=FALSE) {
  if (is.pb) {
    head_foler = "nullpower/pow_subject/"
  } else {
    head_foler = "nullpower/pow_cell/"
  }

  fl <- list.files(head_foler, full.names = TRUE)[grepl(author, list.files(head_foler, full.names = TRUE))]

  if (is.null(n_samples_vec)) {n_samples_vec <- str_extract_all(fl, "(?<=n\\.)\\d+(?=\\.ngenes)") %>% unlist() %>% unique()}

  plots <- lapply(n_samples_vec, function(n_samples) {
    fl_samples <- fl[grepl(n_samples, fl)]

    # if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.fc)") %>% unlist() %>% unique() }
    # if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    #
    # dd <- lapply(ct.indexes, function(ct.index) {
    #   dd_inner <- lapply(genes.values, function(n_genes) {
    #     if (file.exists(paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".fc.0.5.csv"))) {
    #       #d <- read.delim(paste0(head_foler, "/", author ,".n.", n_samples, ".ct.",ct.index,".fc.0.5.csv"), sep = ",")
    #       d <- read.delim(paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".fc.0.5.csv"), sep = ",")
    #       colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")
    #
    #       mask <- colnames(d) %in% algos
    #       mask[1] <- TRUE
    #       d <- d[,mask]
    #
    #       cols <- colnames(d)
    #       d <- lapply(2:ncol(d), function(c) {
    #         values = d[,c] %>% sort(decreasing = TRUE)
    #         values[values <= 1e-300] = 1e-300
    #         x = seq(0,1,length = length(values))
    #         dplyr::tibble(x = x, observed_p_value = -log10(values), name = colnames(d)[c])
    #       }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
    #       return(d)
    #     }
    #   }) %>% do.call("bind_rows", .)
    # }) %>% do.call("bind_rows", .)

    if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.prob)") %>% unlist() %>% unique() }
    if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    iter.values <- str_extract_all(fl_samples, "(?<=iter.)\\d+(?=\\.cs)") %>% unlist() %>% unique()
    if (is.null(pde.values)) {
      pde.values <- lapply(fl_samples, function(l) {
        unlist(strsplit(unlist(strsplit(l, "probde."))[2], ".iter"))[1]
      }) %>% unlist() %>% unique()
    }

    dd <- lapply(ct.indexes, function(ct.index) {
      dd_inner <- lapply(pde.values, function(pde) {
        dd_inner_inner <- lapply(iter.values, function(i.iter) {
          n_genes <- as.integer(as.numeric(pde) * 1000)
          file_name <- paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".probde.",pde,".iter.",i.iter,".csv")
          if (file.exists(file_name)) {
            d <- read.delim(file_name, sep = ",")
            colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")

            mask <- colnames(d) %in% algos
            mask[1] <- TRUE
            d <- d[,mask]

            cols <- colnames(d)
            d <- lapply(2:ncol(d), function(c) {
              values = d[,c] %>% sort(decreasing = TRUE)
              values[values <= 1e-300] = 1e-300
              x = seq(0,1,length = length(values))
              dplyr::tibble(x = x, observed_p_value = -log10(values), name = colnames(d)[c])
            }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
            return(d)
          }
        }) %>% do.call("bind_rows", .)
        dd_inner_inner
      }) %>% do.call('bind_rows', .)
    }) %>% do.call("bind_rows", .)

    colnames(dd)
    dd <- dd %>%
      #group_by(name, ct.index, n.genes) %>%
      group_by(name, n.genes) %>%
      dplyr::arrange(observed_p_value) %>%
      dplyr::mutate(x = row_number() / n())

    if (only_tibble) {
      return(dd %>% dplyr::mutate(n.samples = n_samples))
    }

    p <- dd %>%
      #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
      ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
      #geom_line() +
      geom_point(size = .1) +
      #theme_minimal() +
      theme_bw() +
      scale_y_continuous(trans = "log10") +
      labs(x = "Uniform quantiles", y="Observed -log10(p-value)", col="Algorithm") +
      theme(legend.position = "bottom") +
      #facet_wrap(~ct.index) +
      #ggh4x::facet_nested("Cell type index"+ct.index~"N genes"+n.genes) +
      ggh4x::facet_nested(~paste0(n_samples, " patient")+"N genes"+n.genes) +
      #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n=3))

    p #+ ggtitle(paste("N patients = ", n_samples))
  })

  plots
}

get_result = function(author, is.pb, n_patients = 4, ngenes = 50, cell_index = 1, iter = 2, stop_on_error = T) {
  # Read non-DE res
  if (is.pb) {
    head_folder = "nullpower/null_subject"
  } else {
    head_folder = "nullpower/null_cell"
  }
  
  paths = list.files(head_folder)
  f = grepl(author, paths) & 
    grepl(paste0("n.", n_patients), paths, fixed = T) & 
    grepl(paste0("ngenes.", ngenes, "."), paths, fixed = T) & 
    grepl(paste0("ct.", cell_index), paths) & 
    grepl(paste0("iter.", iter), paths)
  
  if (!sum(f) == 1) {
    if (stop_on_error) {
      stop("found either zero or too many (>2) paths")  
    } else {
      return(dplyr::tibble())
    }
  }
  
  fp = paths[f]
  d <- read.delim(file.path(head_folder, fp), sep = ",")
  d_non_de = d %>% dplyr::mutate(status = 0) %>% 
    tidyr::pivot_longer(!c(X, status), values_to = "p_val")
  
  # Read DE res
  if (is.pb) {
    head_folder = "nullpower/pow_subject/"
  } else {
    head_folder = "nullpower/pow_cell/"
  }
  
  d <- read.delim(file.path(head_folder, fp), sep = ",")
  d_de = d %>% dplyr::mutate(status = 1) %>% 
    tidyr::pivot_longer(!c(X, status), values_to = "p_val")
  
  # Merge results and compute adj.pvalues
  dplyr::bind_rows(d_non_de, d_de) %>%
    dplyr::group_by(name) %>% 
    dplyr::mutate(adj_pval = p.adjust(p_val, method = "BH")) %>% 
    dplyr::group_by(name) %>% 
    dplyr::mutate(X = row_number()) %>% 
    dplyr::mutate(author = author, is.pb = is.pb, n_patients = n_patients, ngenes = ngenes, cell_index = cell_index, iter = iter) %>% 
    dplyr::ungroup()
}

calculate_roc_data <- function(df, status_col, name_col, score_col) {
  unique_names <- unique(df[[name_col]])
  roc_data_list <- list()
  auc_data <- data.frame()
  
  for (name in unique_names) {
    # Filter data for this name
    name_data <- df[df[[name_col]] == name, ]
    
    # Check if we have both classes
    status_values <- unique(name_data[[status_col]])
    if (length(status_values) < 2 || !all(status_values %in% c(0, 1))) {
      cat("Warning:", name, "doesn't have both classes (0 and 1). Skipping.\n")
      next
    }
    
    
    
    # Create ROC object
    roc_obj <- roc(name_data[[status_col]], 
                   1 - name_data[[score_col]], # Higher score = more likely positive
                   direction = "<", 
                   quiet = TRUE)
    
    # Extract coordinates
    roc_coords <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      name = name
    )
    
    roc_data_list[[name]] <- roc_coords
    
    # Store AUC information
    auc_value <- round(auc(roc_obj), 3)
    auc_data <- rbind(auc_data, data.frame(name = name, auc = auc_value))
  }
  
  # Combine all ROC data
  roc_data <- do.call(rbind, roc_data_list)
  
  return(list(roc_data = roc_data, auc_data = auc_data))
}

compute_tpr_fdr <- function(df, threshold) {
  # Predictions: significant if adj_pval < threshold
  predictions <- df$adj_pval < threshold
  ground_truth <- df$status == 1
  
  # Confusion matrix components
  tp <- sum(predictions & ground_truth)        # True Positives
  fp <- sum(predictions & !ground_truth)       # False Positives  
  fn <- sum(!predictions & ground_truth)       # False Negatives
  tn <- sum(!predictions & !ground_truth)      # True Negatives
  
  # Calculate TPR (Sensitivity/Recall) and FDR
  tpr <- ifelse(tp + fn > 0, tp / (tp + fn), 0)  # TPR = TP / (TP + FN)
  fdr <- ifelse(tp + fp > 0, fp / (tp + fp), 0)  # FDR = FP / (TP + FP)
  
  return(list(tpr = tpr, fdr = fdr, tp = tp, fp = fp, fn = fn, tn = tn))
}

plot_pr_with_alpha <- function(df, tool_names, alpha = 0.05, show_alpha_point = TRUE) {
  pr_data_list <- list()
  auc_data <- data.frame()
  
  for (tool in tool_names) {
    d <- df %>% dplyr::filter(name == tool)
    
    # ---- PR curve from ranking ----
    pred_scores <- 1 - d$adj_pval
    o <- order(pred_scores, decreasing = TRUE)
    y <- d$status[o]                   # 1 = positive, 0 = negative
    tp <- cumsum(y)
    fp <- cumsum(1 - y)
    precision <- tp / (tp + fp)
    precision[is.nan(precision)] <- 0
    recall <- tp / sum(y)
    
    pr_coords <- data.frame(
      Precision = precision,
      Recall = recall,
      tool = tool
    )
    pr_data_list[[tool]] <- pr_coords
    
    # Trapezoidal AUC-PR on the stepwise curve
    auc_pr <- sum(diff(recall) * precision[-1], na.rm = TRUE)
    auc_data <- rbind(auc_data, data.frame(tool = tool, auc_pr = round(auc_pr, 3)))
  }
  
  pr_data <- do.call(rbind, pr_data_list)
  pr_data <- merge(pr_data, auc_data, by = "tool")
  pr_data$label <- paste0(pr_data$tool, " (AUC = ", pr_data$auc_pr, ")")
  
  # ---- Alpha (e.g., 0.05) point per tool ----
  pr_points <- NULL
  if (show_alpha_point) {
    pr_points <- lapply(tool_names, function(tool) {
      d <- df %>% dplyr::filter(name == tool)
      y <- d$status == 1
      pred <- d$adj_pval <= alpha
      TP <- sum(pred & y)
      FP <- sum(pred & !y)
      P  <- sum(y)
      prec <- ifelse(TP + FP > 0, TP / (TP + FP), NA_real_)
      rec  <- ifelse(P > 0, TP / P, NA_real_)
      data.frame(tool = tool, Precision = prec, Recall = rec)
    }) %>% dplyr::bind_rows() %>% dplyr::left_join(auc_data, by = "tool") %>%
      dplyr::mutate(label = paste0(tool, " (AUC = ", auc_pr, ")"),
                    alpha = alpha)
  }
  
  # ---- Plot ----
  p <- ggplot(pr_data, aes(x = Recall, y = Precision, color = label)) +
    geom_line(size = 1.2) +
    scale_color_brewer(type = "qual", palette = "Set2") +
    labs(title = "Precision-Recall Curves",
         x = "Recall (Sensitivity)", y = "Precision (PPV)",
         color = "Tool") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  if (show_alpha_point && !is.null(pr_points)) {
    p <- p +
      geom_point(data = pr_points, aes(x = Recall, y = Precision, color = label),
                 size = 3, stroke = 1.1) +
      ggrepel::geom_text_repel(
        data = pr_points,
        aes(x = Recall, y = Precision,
            label = paste0("α=", format(pr_points$alpha[1])), color = label),
        size = 3, show.legend = FALSE, max.overlaps = 20
      )
  }
  
  p
}


create_other_curves <- function(df, tool_names) {
  pr_data_list <- list()
  auc_data <- data.frame()
  
  tool = tool_names[1]
  r = lapply(tool_names, function(tool) {
    d = df %>% dplyr::filter(name == tool)
    
    # Calculate precision-recall
    pred_scores <- 1 - d$adj_pval
    
    # Sort by prediction score (descending)
    sorted_indices <- order(pred_scores, decreasing = TRUE)
    true_labels <- d$status[sorted_indices]
    xs = cumsum(true_labels == 0)
    dplyr::tibble(x = 1:length(xs), y = xs, tool = tool)
  }) %>% do.call("bind_rows", .)
  
  ggplot(r, mapping = aes(x = x, y = y, colour = tool)) +
    geom_line() +
    scale_y_continuous(transform = "log10") +
    theme_bw()
  
  p <- ggplot(pr_data, aes(x = Recall, y = Precision, color = label)) +
    geom_line(size = 1.2) +
    scale_color_brewer(type = "qual", palette = "Set2") +
    labs(title = "Precision-Recall Curves",
         x = "Recall (Sensitivity)", y = "Precision (PPV)",
         color = "Tool") +
    theme_bw() +
    theme(legend.position = "bottom")
  p
  return(p)
}



plot_tpr_fdr = function(author, is.pb, n_patients, ngenes, cell_index, iter) {
  df = get_result(author, is.pb, n_patients, ngenes, cell_index, iter) %>% 
    dplyr::mutate(X = paste0(X, author, is.pb, n_patients, ngenes, cell_index, iter))       
  
  df = df %>% na.omit()

  # Prep cobra object
  pval = df %>%
    dplyr::select(X, name, p_val) %>%
    tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
    dplyr::arrange(X) %>%
    dplyr::select(!X) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))

  truth = df %>%
    dplyr::select(X, status) %>%
    dplyr::distinct() %>%
    dplyr::arrange(X) %>%
    dplyr::select(!X) %>%
    as.data.frame()

  rownames(pval) = paste0("Gene", 1:nrow(pval))
  rownames(truth) = paste0("Gene", 1:nrow(pval))
  truth$feature = rownames(truth)

  library(iCOBRA)
  cobradata = iCOBRA::COBRAData(pval, truth = truth)
  cobradata <- calculate_adjp(cobradata, method = "BH")
  cobraperf <- calculate_performance(cobradata, binary_truth = "status",
                                     cont_truth = "logFC", splv = "none",
                                     maxsplit = 4)

  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2")
  plot_fdrtprcurve(cobraplot)
  
  # Check if required columns exist
  status_col <- "status"
  name_col <- "name"
  score_col <- "adj_pval"
  
  required_cols <- c(status_col, name_col, score_col)
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("Missing columns:", paste(setdiff(required_cols, colnames(df)), collapse = ", ")))
  }
  
  # Calculate ROC data
  roc_results <- calculate_roc_data(df, status_col, name_col, score_col)
  roc_data <- roc_results$roc_data
  auc_data <- roc_results$auc_data
  
  # Create labels with AUC values
  auc_data$label <- paste0(auc_data$name, " (AUC = ", auc_data$auc, ")")
  
  # Merge labels back to roc_data
  roc_data <- merge(roc_data, auc_data[, c("name", "label")], by = "name")
  
  # Get number of unique names for color palette
  n_names <- length(unique(roc_data$name))
  
  # Prepare colors
  if (n_names <= 11) {
    colors <- RColorBrewer::brewer.pal(max(3, n_names), "Spectral")
  } else {
    colors <- rainbow(n_names)
  }
  
  # Create the ggplot
  p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = label)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed", alpha = 0.7) +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = "ROC Curves by Name",
      x = "False Positive Rate (FPR)",
      y = "True Positive Rate (TPR)",
      color = "Method"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
      panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", size = 0.25),
      axis.title = element_text(face = "bold"),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 1.5)))
  
  p

  dplyr::left_join(
    df %>%
      dplyr::mutate(TP = adj_pval <= .05 & status == 1) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(TP = sum(TP)),
    df %>%
      dplyr::mutate(FP = adj_pval <= .05 & status == 0) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(FP = sum(FP))
  )
  # 
  # 
  # 
  # #df = df %>% na.omit()
  # 
  # # Prep cobra object
  # pval = df %>% 
  #   dplyr::select(X, name, p_val) %>% 
  #   tidyr::pivot_wider(values_from = p_val, names_from = name) %>% 
  #   dplyr::arrange(X) %>% 
  #   dplyr::select(!X) %>% 
  #   as.data.frame()
  # rownames(pval) = paste0("Gene", 1:nrow(pval))
  # 
  # truth = df %>% 
  #   dplyr::select(X, status) %>% 
  #   dplyr::distinct() %>% 
  #   dplyr::arrange(X) %>% 
  #   dplyr::select(!X) %>% 
  #   as.data.frame()
  # 
  # rownames(pval) = paste0("Gene", 1:nrow(pval))
  # rownames(truth) = paste0("Gene", 1:nrow(pval)) 
  # truth$feature = rownames(truth)
  # 
  # library(iCOBRA)
  # cobradata = iCOBRA::COBRAData(pval, truth = truth)
  # cobradata <- calculate_adjp(cobradata, method = "BH")
  # cobraperf <- calculate_performance(cobradata, binary_truth = "status", 
  #                                    cont_truth = "logFC", splv = "none",
  #                                    maxsplit = 4)
  # 
  # cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", facetted = TRUE)
  # plot_fdrtprcurve(cobraplot)
}

map_name = c(
  "edgeR..Pb." = "edgeR (Pb)",
  "edgeR..cell." = "edgeR (cell)",
  "MAST..cell." = "MAST (cell)",
  "Seurat..cell." = "Seurat (cell)",
  "limma..Pb." = "limma (Pb)",
  "limma..cell." = "limma (cell)",
  "glmGamPoi..Pb." = "glmGamPoi (Pb)",
  "glmGamPoi..cell." = "glmGamPoi (cell)",
  "NEBULA" = "Nebula",
  "devil..base." = "Devil (base)",
  "devil..mixed." = "Devil (mixed)",
  "devil..sf.base." = "Devil (SF/base)",
  "devil..sf.mixed." = "Devil (SF/mixed)"
)

author = "hsc"
plot_pvals = function(author, is.pb, iters = c(1:5), ngenes = c(12, 25, 50), ct.indexes = 1:6, n.patients = 8) {
  df = lapply(iters, function(i) {
    lapply(ngenes, function(ng) {
      lapply(ct.indexes, function(ct) {
        lapply(n.patients, function(np) {
          get_result(author, is.pb = is.pb, n_patients = np, ngenes = ng, cell_index = ct, iter = i, stop_on_error = F)
        }) %>% do.call("bind_rows", .)
      }) %>% do.call("bind_rows", .)
    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .)
  
  df$name = map_name[df$name]
  
  df %>% 
    dplyr::filter(status == 0) %>% 
    # dplyr::filter(name %in% method_patientwise) %>% 
    ggplot(mapping = aes(x = p_val)) +
    geom_histogram(binwidth = .01) +
    facet_wrap(~name, scales = "free_y")
  
  df %>% 
    dplyr::filter(status == 0) %>% 
    group_by(name, iter, n_patients, ngenes, cell_index) %>%
    arrange(p_val, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1)
    ) %>% 
    ggplot(mapping = aes(x = -log10(expected), y =-log10(p_val), col = name)) +
    # ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
    geom_smooth() +
    #geom_point() +
    theme_bw() +
    #ggsci::scale_color_nejm() +
    scale_x_continuous(limits = c(0, 5)) +
    scale_y_continuous(limits = c(0, 5)) +
    geom_abline(slope = 1, intercept = 0) +
    coord_fixed()
  
  df %>% 
    dplyr::filter(status == 0) %>% 
    group_by(name, iter, n_patients, ngenes, cell_index) %>%
    arrange(p_val, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1)
    ) %>% 
    ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
    # ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
    #geom_smooth() +
    geom_point() +
    theme_bw() +
    facet_wrap(~name)
  
  df %>% 
    dplyr::filter(status == 0) %>% 
    group_by(name, iter, n_patients, ngenes, cell_index) %>%
    arrange(p_val, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1)
    ) %>% 
    ggplot(mapping = aes(x = -log10(expected), y =-log10(p_val), col = name)) +
    # ggplot(mapping = aes(x = expected, y = p_val, col = name)) +
    #geom_smooth() +
    geom_point() +
    theme_bw() +
    #ggsci::scale_color_nejm() +
    scale_x_continuous(limits = c(0, 5)) +
    scale_y_continuous(limits = c(0, 5)) +
    geom_abline(slope = 1, intercept = 0) +
    coord_fixed()
  
  
  # TPR vs FDR curve
  df = df %>% na.omit()
  df = df %>% 
    dplyr::filter(n_patients == 20, iter == 1, cell_index == 1, ngenes == 50) %>% 
    dplyr::filter(name %in% method_patientwise)
  
  # Prep cobra object
  pval = df %>%
    dplyr::select(X, name, p_val) %>%
    tidyr::pivot_wider(values_from = p_val, names_from = name) %>%
    dplyr::arrange(X) %>%
    dplyr::select(!X) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))
  
  padj = df %>%
    dplyr::select(X, name, adj_pval) %>%
    tidyr::pivot_wider(values_from = adj_pval, names_from = name) %>%
    dplyr::arrange(X) %>%
    dplyr::select(!X) %>%
    as.data.frame()
  rownames(pval) = paste0("Gene", 1:nrow(pval))
  
  truth = df %>%
    dplyr::select(X, status) %>%
    dplyr::distinct() %>%
    dplyr::arrange(X) %>%
    dplyr::select(!X) %>%
    as.data.frame()
  
  df %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(TP = sum(adj_pval <= .05 & status == 1))
  
  df %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(FP = sum(adj_pval <= .05 & status == 0))
  
  rownames(pval) = paste0("Gene", 1:nrow(pval))
  rownames(padj) = paste0("Gene", 1:nrow(pval))
  rownames(truth) = paste0("Gene", 1:nrow(pval))
  truth$feature = rownames(truth)
  
  library(iCOBRA)
  cobradata = iCOBRA::COBRAData(pval, truth = truth, padj = padj)
  #cobradata <- calculate_adjp(cobradata, method = "BH")
  cobraperf <- calculate_performance(cobradata, binary_truth = "status")
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2")
  plot_fdrtprcurve(cobraplot)
  plot_fpr(cobraplot)
  plot_tpr(cobraplot)
  
  
}


plot_pvalues = function(author, method_cellwise, method_patientwise) {
  d1 <- all_null_plots(author, FALSE, algos = method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble=TRUE)[[1]] %>%
    dplyr::mutate(ytype = "p-value", xtype="Cell-wise")
  d2 <- all_null_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "p-value", xtype="Patient-wise")
  d3 <- all_pow_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "-log10 p-value", xtype="Cell-wise")
  d4 <- all_pow_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "-log10 p-value", xtype="Patient-wise")
  
  pA = dplyr::bind_rows(d1, d2) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>% 
    filter(!is.na(observed_p_value)) %>%
    group_by(name) %>%
    arrange(observed_p_value, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      expected = rank / (n() + 1)
    ) %>% 
    ggplot(aes(x = expected, y = observed_p_value, color = name)) +
    geom_point(size = 0.1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(
      x = "Expected uniform quantiles",
      y = "Observed p-values"
    ) +
    theme_bw() +
    facet_wrap(~xtype)

  # dplyr::bind_rows(d1, d2) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  #   dplyr::group_by(xtype, ytype) %>%
  #   dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  #   ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
  #   ggridges::geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = F, alpha = .7) +
  #   scale_color_manual(values = sort(method_colors)) +
  #   scale_fill_manual(values = sort(method_colors)) +
  #   facet_grid(~xtype, scales = "free") +
  #   theme_bw() +
  #   labs(x = "p-value", y="", col="Algorithm") +
  #   scale_color_manual(values = method_colors) +
  #   #facet_wrap(~paste0(n.genes, " genes")) +
  #   theme(legend.position = "bottom") +
  #   scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
  #   theme(legend.position = "none") +
  #   theme()

  # dplyr::bind_rows(d3, d4) %>%
  #   dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  #   dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  #   dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  #   dplyr::group_by(xtype, ytype) %>%
  #   dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  #   ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
  #   ggridges::geom_density_ridges(stat = "binline", bins = 100, scale = 0.95, draw_baseline = F, alpha = .7) +
  #   scale_color_manual(values = sort(method_colors)) +
  #   scale_fill_manual(values = sort(method_colors)) +
  #   facet_grid(~xtype, scales = "free") +
  #   theme_bw() +
  #   labs(x = bquote(-log[10]~ "(p-value)"), y="", col="Algorithm") +
  #   scale_color_manual(values = method_colors) +
  #   #facet_wrap(~paste0(n.genes, " genes")) +
  #   theme(legend.position = "bottom") +
  #   scale_x_continuous(transform = "log10") +
  #   #scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
  #   scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
  #   theme(legend.position = "none") +
  #   theme()


  # pB <- dplyr::bind_rows(d1, d2) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  #   dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  #   dplyr::group_by(xtype, ytype) %>%
  #   dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  #   ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
  #   ggridges::geom_density_ridges(alpha = .7, scale = 1) +
  #   #ggridges::geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = F) +
  #   scale_color_manual(values = sort(method_colors)) +
  #   scale_fill_manual(values = sort(method_colors)) +
  #   facet_grid(~xtype, scales = "free") +
  #   theme_bw() +
  #   labs(x = "p-value", y="", col="Algorithm") +
  #   scale_color_manual(values = method_colors) +
  #   #facet_wrap(~paste0(n.genes, " genes")) +
  #   theme(legend.position = "bottom") +
  #   scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
  #   scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
  #   theme(legend.position = "none") +
  #   theme()

  pC <- dplyr::bind_rows(d3, d4) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = bquote(-log[10]~ "(p-value)"), y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(transform = "log10") +
    #scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme()

  list(null_pvalue = pB, de_pvalue = pC)
}

plot_MCCs = function(author, method_cellwise, method_patientwise) {
  res <- readRDS("nullpower/final_res/results.rds")
  a <- author

  pfalse <- res %>%
    na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::filter(author == a) %>%
    dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = patients, Ngenes = paste0(ngenes, " genes")) %>%
    #dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
    dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
    dplyr::group_by(Npatients, ngenes) %>%
    ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of DE genes", y = "MCC", col="Model", linetype = "N patients", shape = "N patients") +
    ggh4x::facet_nested(~"Number of samples"+Npatients, scales = "free_y") +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "bottom")

  ptrue <- res %>%
    na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::filter(author == a) %>%
    dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = patients, Ngenes = paste0(ngenes, " genes")) %>%
    #dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
    dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
    dplyr::group_by(Npatients, ngenes) %>%
    ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of DE genes", y = "MCC", col="Model", linetype = "N patients", shape = "N patients") +
    ggh4x::facet_nested(~"Number of samples"+Npatients, scales = "free_y") +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "bottom")

  list(cellwise = pfalse, patientwise = ptrue)
}

plot_ks = function(author, method_cellwise, method_patientwise) {
  res <- readRDS("nullpower/final_res/results.rds")
  a = author

  {
    dx <- .15

    r_ks <- res %>%
      na.omit() %>%
      dplyr::filter(author == a) %>%
      dplyr::filter(is.pb == FALSE, name %in% method_cellwise)

    r_ks_tot <- lapply(unique(r_ks$name), function(n) {
      MCCs <- r_ks %>%
        dplyr::filter(name==n) %>%
        dplyr::pull(MCC)
      ECDFs <- lapply(MCCs, function(mcc) {
        sum(MCCs <= mcc) / length(MCCs)
      }) %>% unlist()

      dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
    }) %>% do.call("rbind", .)

    ks_pvals <- lapply(method_cellwise[method_cellwise!="Devil (base)"], function(m) {
      print(m)
      pval = ks.test(filter(r_ks_tot, name=="Devil (base)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
      if (pval == 0) {
        pval <- "< 2e-16"
      } else {
        pval <- paste0("= ", round(pval, 2))
      }
      dplyr::tibble(m=m, pval=pval)
    }) %>% do.call("bind_rows", .)

    ks_false <- r_ks_tot %>%
      dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
      dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
      dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
      dplyr::mutate(name = factor(name, levels = method_levels)) %>%
      ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
      geom_line(linewidth = .8, position = position_dodge(width = .02)) +
      scale_color_manual(values = sort(method_colors)) +
      #ggtitle(paste0("Cell-wise")) +
      labs(x = "MCC", y = "Empirical CDF", col="Model") +
      theme_bw()

    L = ggplot_build(ks_false)$layout$panel_params[[1]]
    Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
    Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

    ks_false <- ks_false +
      annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
      annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
      annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
      theme(legend.position = "left")
  }

  # kolmogorov smirnof plots - TRUE ####
  {
    r_ks <- res %>%
      na.omit() %>%
      dplyr::filter(author == a) %>%
      dplyr::filter(is.pb == TRUE, name %in% method_patientwise)

    r_ks_tot <- lapply(unique(r_ks$name), function(n) {
      MCCs <- r_ks %>%
        dplyr::filter(name==n) %>%
        dplyr::pull(MCC)
      ECDFs <- lapply(MCCs, function(mcc) {
        sum(MCCs <= mcc) / length(MCCs)
      }) %>% unlist()

      dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
    }) %>% do.call("rbind", .)

    ks_pvals <- lapply(method_patientwise[method_patientwise!="Devil (mixed)"], function(m) {
      print(m)
      pval = ks.test(filter(r_ks_tot, name=="Devil (mixed)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
      true_pval <- pval
      if (pval <= 1e-6) {
        pval <- "< 1e-6"
      } else {
        pval <- paste0("= ", round(pval, 3))
      }
      dplyr::tibble(m=m, pval=pval, true_pval=true_pval)
    }) %>% do.call("bind_rows", .)

    ks_true <- r_ks_tot %>%
      dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
      dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
      dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
      dplyr::mutate(name = factor(name, levels = method_levels)) %>%
      ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
      geom_line(linewidth = .8, position = position_dodge(width = .02)) +
      scale_color_manual(values = method_colors) +
      #ggtitle(paste0("Patient-wise")) +
      #ggtitle(paste0("Cell-wise")) +
      labs(x = "MCC", y = "Empirical CDF", col="Model") +
      theme_bw()

    L = ggplot_build(ks_true)$layout$panel_params[[1]]
    Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
    Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

    ks_true <- ks_true +
      annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
      annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
      annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
      theme(legend.position = "left")
  }

  list(cellwise = ks_false, patientwise = ks_true)
}

plot_timing = function(author, competitor = "NEBULA") {
  a = author
  
  lfs = list.files("nullpower/timing_results/")[grepl(paste0(a, "_"), list.files("nullpower/timing_results/"))]
  df_time = lapply(lfs, function(lf) {
    readRDS(paste0("nullpower/timing_results/", a,".rds")) %>%
      dplyr::filter(algo %in% c("Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
      dplyr::rename(name = algo) %>%
      dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
      dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
      dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
      dplyr::rename(algo = name) %>%
      dplyr::group_by(author, is.pb, n.sample, n.gene, int.ct, iter) %>%
      dplyr::mutate(time_fold = timings[algo == competitor] / timings) %>%
      dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
      dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k")))
  }) %>% do.call("bind_rows", .)
  
  
  
  timing_plot <- df_time %>% 
    ggplot(mapping = aes(x=cell_order, y=time_fold, col=algo)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of cells", y=paste0("Speedup (vs. ",competitor,")"), col = "Model") +
    theme_bw() +
    geom_hline(yintercept = 1, color = "darkslategray", linetype = 'dashed') +
    coord_flip() +
    scale_y_continuous(transform = "log10")
  timing_plot
}


plot_all_models = function() {
  method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)")
  method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)")
  method_colors = c(
    "glmGamPoi (fixed)" = "#A22E29",
    "edgeR" = "#7D629E",
    "edgeR (Pb)" = "#7D629E",
    "limma" = "#B96461",
    "limma (Pb)" = "#B96461",
    "glmGamPoi (cell)" = "#EAB578",
    "glmGamPoi" = "#EAB578",
    "Nebula" =  'steelblue', #"#B0C4DE",
    "NEBULA" =  'steelblue', #"#B0C4DE",
    "Devil (base)" = "#099668",
    "Devil (mixed)" = "#099668",
    "Devil" = "#099668",
    "devil" = "#099668"
  )

  res = readRDS("nullpower/final_res/results.rds") %>%
    dplyr::filter((is.pb == TRUE & name %in% method_patientwise) | (is.pb == FALSE & name %in% method_cellwise))

  all_timing = lapply(list.files("nullpower/timing_results/"), function(p) {
    readRDS(file.path("nullpower/timing_results/", p)) %>% dplyr::mutate(dataset = p)
  }) %>% do.call("bind_rows", .) %>%
    dplyr::filter((is.pb == TRUE & algo %in% method_patientwise) | (is.pb == FALSE & algo %in% method_cellwise))

  all_MCC_boxplots = res %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    ggplot(mapping = aes(x=name, y=MCC, col=name)) +
    geom_boxplot() +
    ggh4x::facet_nested("Dataset"+author~is.pb) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(color = "Model") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )

  all_time_boxplots = all_timing %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = algo) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    ggplot(mapping = aes(x=name, y=timings, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(color = "Model") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_y_continuous(transform = "log10") +
    coord_flip() +
    facet_wrap(~author) +
    labs(y = "Runtime (seconds)")

  failure_rate_plot = all_timing %>%
    dplyr::filter(is.pb) %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = algo) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    dplyr::mutate(is_bad = is.na(timings)) %>%
    dplyr::group_by(name, is.pb, author) %>%
    dplyr::summarise(`Failure rate` = sum(is_bad) / n()) %>%
    ggplot(mapping = aes(x=author, y=`Failure rate`, fill=name, col=name)) +
    geom_col(position = "dodge") +
    ylim(c(0,1)) +
    ggh4x::facet_nested(~is.pb) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    labs(color = "", fill="", x = "Dataset name") +
    theme_bw()

  # failure_rate_plot = all_timing %>%
  #   dplyr::filter(is.pb) %>%
  #   dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
  #   dplyr::mutate(name = algo) %>%
  #   dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
  #   dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
  #   dplyr::mutate(is_bad = is.na(timings)) %>%
  #   dplyr::group_by(name, is.pb, author) %>%
  #   dplyr::summarise(`Failure rate` = sum(is_bad) / n()) %>%
  #   ggplot(mapping = aes(x=name, y=`Failure rate`, fill=name, col=name)) +
  #   geom_col(col="black") +
  #   ylim(c(0,1)) +
  #   ggh4x::facet_nested("Dataset"+author~is.pb) +
  #   scale_fill_manual(values = method_colors) +
  #   labs(color = "") +
  #   theme_bw() +
  #   coord_flip() +
  #   theme(
  #     legend.position = "none",
  #     axis.title.y = element_blank()
  #   )

  list(MCC = all_MCC_boxplots, timing = all_time_boxplots, failure_rate = failure_rate_plot)
}

plot_MCCs_boxplots = function(a = NULL) {
  res = readRDS("nullpower/final_res/results.rds")
  if (!is.null(a)) {
    res = res %>% dplyr::filter(author == a)
  }
  
  res$name %>% unique()
  
  df_cellwise = res %>%
    na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_patientwise = res %>%
    na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_all = dplyr::bind_rows(df_patientwise, df_cellwise)
  df_all %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    #dplyr::summarise(MCC = mean(MCC)) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = name, y=MCC, col=name)) +
    geom_boxplot() +
    geom_point() +
    coord_flip() +
    ggh4x::facet_nested(ct.index~is.pb+paste0(patients, " patients")) +
    # scale_y_continuous(labels = custom_labels) +
    ggsci::scale_color_nejm() +
    #scale_color_manual(values = method_colors) +
    theme_bw() +
    labs(x = "", col = "")
}


author = "hsc"
is.pb = F
n_patients = 8
ngenes = 12
cell_index = 1
iter = 1
df = get_result(author, is.pb, n_patients, ngenes, cell_index, iter)


