
require(ggupset)

method_colors = c(
  "glmGamPoi - cpu" = "#EAB578",
  "nebula" =  "steelblue",
  "NEBULA" =  "steelblue",
  "devil - cpu" = "#099668",
  "devil - h100" = "#384B41",
  "devil - a100" = "#9BB0A5"
)

get_time_results <- function(results_folder) {
  results_paths <- list.files(results_folder)
  results_paths <- results_paths[!results_paths %in% c("fits", "memory", "large", "no_overdispersion")]

  device_path = "a100"

  results <- lapply(results_paths, function(device_path) {
    lapply(list.files(file.path(results_folder, device_path)), function(p){
      if (!grepl(".rds", p)) return(dplyr::tibble())
      info <- unlist(strsplit(p, "_"))
      res = readRDS(file.path(results_folder, device_path, p))
      dplyr::tibble(
        model_name = paste(info[1], info[2]),
        n_genes = info[3],
        n_cells = info[5],
        n_cell_types = info[7],
        time = as.numeric(unlist(res$time), units = "secs")#,
        #memory = res$mem_alloc
      ) %>% dplyr::mutate(model_name = ifelse(grepl("gpu" ,model_name), device_path, model_name))
    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .) %>%
    dplyr::mutate(n_cells = as.numeric(n_cells), n_genes = as.numeric(n_genes))

  results %>%
    dplyr::mutate(model_name = dplyr::recode(model_name,
                                             "cpu devil" = "devil - cpu",
                                             "cpu glmGamPoi" = "glmGamPoi - cpu",
                                             "a100" = "devil - a100",
                                             "h100" = "devil - h100")) %>%
    dplyr::mutate(Measure = "observed")
}

# get_time_results <- function(results_folder) {
#   # top-level folders: e.g. "a100", "h100", "cpu", "fits", ...
#   results_paths <- list.files(results_folder)
#   results_paths <- results_paths[!results_paths %in% c("fits", "memory", "large", "no_overdispersion")]
#
#   results <- lapply(results_paths, function(device) {
#     device_dir <- file.path(results_folder, device)
#     lvl2 <- list.files(device_dir)
#
#     # Case 1: CPU-like layout -> .rds files directly inside device_dir
#     if (any(grepl("\\.rds$", lvl2))) {
#       files <- lvl2[grepl("\\.rds$", lvl2)]
#
#       lapply(files, function(p) {
#         info <- strsplit(p, "_")[[1]]
#         res  <- readRDS(file.path(device_dir, p))
#
#         dplyr::tibble(
#           device        = device,
#           n_gpus        = NA_integer_,
#           model_name    = paste(info[1], info[2]),   # "cpu devil", "cpu glmGamPoi"
#           n_genes       = info[3],
#           n_cells       = info[5],
#           n_cell_types  = info[7],
#           time          = as.numeric(res$time[[1]], units = "secs")
#         )
#       }) %>% dplyr::bind_rows()
#
#     } else {
#       # Case 2: GPU-like layout -> subdirs "1_gpus", "8_gpus", each with .rds files
#       gpu_folders <- lvl2
#
#       lapply(gpu_folders, function(gpu_folder) {
#         gpu_dir <- file.path(device_dir, gpu_folder)
#         files   <- list.files(gpu_dir, pattern = "\\.rds$")
#         n_gpu_val <- as.integer(strsplit(gpu_folder, "_")[[1]][1])  # "1_gpus" -> 1
#
#         lapply(files, function(p) {
#           info <- strsplit(p, "_")[[1]]
#           res  <- readRDS(file.path(gpu_dir, p))
#
#           dplyr::tibble(
#             device        = device,                    # "a100" or "h100"
#             n_gpus        = n_gpu_val,                 # 1 or 8
#             model_name    = paste(info[1], info[2]),   # "gpu devil"
#             n_genes       = info[3],
#             n_cells       = info[5],
#             n_cell_types  = info[7],
#             time          = as.numeric(res$time[[1]], units = "secs")
#           )
#         }) %>% dplyr::bind_rows()
#       }) %>% dplyr::bind_rows()
#     }
#   }) %>%
#     dplyr::bind_rows() %>%
#     dplyr::mutate(
#       n_cells = as.numeric(n_cells),
#       n_genes = as.numeric(n_genes)
#     )
#
#   results %>%
#     # For GPU runs, replace "gpu devil" with the device name ("a100"/"h100")
#     dplyr::mutate(
#       model_name = ifelse(grepl("^gpu", model_name), device, model_name)
#     ) %>%
#     dplyr::mutate(
#       model_name = dplyr::recode(
#         model_name,
#         "cpu devil"     = "devil - cpu",
#         "cpu glmGamPoi" = "glmGamPoi - cpu",
#         "a100"          = "devil - a100",
#         "h100"          = "devil - h100"
#       ),
#       Measure = "observed"
#     )
# }


predict_time_results = function(results) {
  unique_nc = unique(results$n_cells)
  unique_ng = unique(results$n_genes)
  unique_m = unique(results$model_name)

  lapply(unique_m, function(m) {
    lapply(unique_nc, function(nc) {
      d = results %>% dplyr::filter(n_cells==nc, model_name==m)
      if (length(unique(d$n_genes)) < 3) {
        n_genes_to_predict = unique_ng[!unique_ng %in% d$n_genes]

        ratios = results %>%
          dplyr::filter(n_genes == n_genes_to_predict) %>%
          #dplyr::filter(n_cells == nc) %>%
          dplyr::group_by(model_name, n_genes, n_cells) %>%
          dplyr::summarise(time = mean(time)) %>%
          dplyr::filter(n_cells != nc) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(ratio = time[model_name == m] / time)

        average_ratios = ratios %>% dplyr::filter(model_name != m) %>%
          dplyr::group_by(model_name) %>%
          dplyr::summarise(average_ratio = median(ratio))
        average_ratios_vec = c(average_ratios$average_ratio)
        names(average_ratios_vec) = average_ratios$model_name

        time_predicted = results %>%
          dplyr::filter(n_genes == n_genes_to_predict, n_cells == nc) %>%
          dplyr::select(model_name, time) %>%
          dplyr::mutate(predicted_time = time * average_ratios_vec[model_name]) %>%
          dplyr::pull(predicted_time) %>%
          mean()

        # ratios = results %>% dplyr::filter(n_genes == n_genes_to_predict, model_name==m | model_name == m_compare) %>%
        #   dplyr::group_by(model_name, n_genes, n_cells) %>%
        #   dplyr::summarise(time = mean(time)) %>%
        #   dplyr::filter(n_cells != nc) %>%
        #   dplyr::ungroup() %>%
        #   dplyr::mutate(ratio = time / time[model_name == m_compare]) %>%
        #   dplyr::filter(model_name == m) %>%
        #   dplyr::pull(ratio)
        # ratio = mean(ratios)
        #
        # time_predicted = ratio * (results %>% dplyr::filter(n_genes == n_genes_to_predict, model_name == m_compare, n_cells == nc) %>%
        #   dplyr::group_by(model_name, n_genes, n_cells) %>%
        #   dplyr::summarise(time = mean(time)) %>%
        #   dplyr::pull(time))

        # lin_mod = lm(time~n_genes * n_cells, data = results %>% dplyr::filter(model_name == m))
        # n_genes_to_predict = unique_ng[!unique_ng %in% d$n_genes]
        # time_predicted = predict(lin_mod, data.frame(n_genes = n_genes_to_predict, n_cells = nc))
        #
        # lin_mod = lm(time~n_genes, data = d)
        # n_genes_to_predict = unique_ng[!unique_ng %in% d$n_genes]
        # time_predicted = predict(lin_mod, data.frame(n_genes = n_genes_to_predict, n_cells = nc))

        dplyr::tibble(model_name = m, n_genes = n_genes_to_predict, n_cells = nc, n_cell_types = unique(d$n_cell_types), time = time_predicted, Measure = "predicted")
      }

    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .)
}


get_memory_results <- function(results_folder) {
  results_paths <- list.files(results_folder)
  results_paths <- results_paths[!results_paths %in% c("fits", "memory", "large", "no_overdispersion")]

  device_path = "h100"

  results <- lapply(results_paths, function(device_path) {
    lapply(list.files(file.path(results_folder, device_path)), function(p){
      if (!grepl(".rds", p)) return(NULL)
      info <- unlist(strsplit(p, "_"))
      res = readRDS(file.path(results_folder, device_path, p))
      dplyr::tibble(
        model_name = paste(info[1], info[2]),
        n_genes = info[3],
        n_cells = info[5],
        n_cell_types = info[7],
        #time = as.numeric(unlist(res$time), units = "secs")#,
        memory = res$mem_alloc / 1e9
      ) %>% dplyr::mutate(model_name = ifelse(grepl("gpu" ,model_name), device_path, model_name))
    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .) %>%
    dplyr::mutate(n_cells = as.numeric(n_cells), n_genes = as.numeric(n_genes))

  results %>%
    dplyr::mutate(model_name = dplyr::recode(model_name,
                                             "cpu devil" = "devil - cpu",
                                             "cpu glmGamPoi" = "glmGamPoi - cpu",
                                             "a100" = "devil - a100",
                                             "h100" = "devil - h100")) %>%
    dplyr::mutate(Measure = "observed")
}

predict_memory_results = function(results) {
  unique_nc = unique(results$n_cells)
  unique_ng = unique(results$n_genes)
  unique_m = unique(results$model_name)

  m = unique_m[2]
  nc = unique_nc[3]
  lapply(unique_m, function(m) {
    lapply(unique_nc, function(nc) {
      d = results %>% dplyr::filter(n_cells==nc, model_name==m)
      if (length(unique(d$n_genes)) < 3) {
        n_genes_to_predict = unique_ng[!unique_ng %in% d$n_genes]

        ratios = results %>%
          dplyr::filter(n_genes == n_genes_to_predict) %>%
          #dplyr::filter(n_cells == nc) %>%
          dplyr::group_by(model_name, n_genes, n_cells) %>%
          dplyr::summarise(memory = mean(memory)) %>%
          dplyr::filter(n_cells != nc) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(ratio = memory[model_name == m] / memory)

        average_ratios = ratios %>% dplyr::filter(model_name != m) %>%
          dplyr::group_by(model_name) %>%
          dplyr::summarise(average_ratio = median(ratio))
        average_ratios_vec = c(average_ratios$average_ratio)
        names(average_ratios_vec) = average_ratios$model_name

        memory_predicted = results %>%
          dplyr::filter(n_genes == n_genes_to_predict, n_cells == nc) %>%
          dplyr::select(model_name, memory) %>%
          dplyr::mutate(predicted_memory = memory * average_ratios_vec[model_name]) %>%
          dplyr::pull(predicted_memory) %>%
          mean()

        dplyr::tibble(model_name = m, n_genes = n_genes_to_predict, n_cells = nc, n_cell_types = unique(d$n_cell_types), memory = memory_predicted,Measure = "predicted")
      }

    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .)
}

time_comparison = function(results, ratio = FALSE) {
  if (!isFALSE(ratio)) {
    df = results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name,Measure) %>%
      dplyr::summarise(y = mean(time)) %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types) %>%
      dplyr::mutate(y = y[model_name == ratio] / y) %>%
      dplyr::group_by(n_genes, n_cells) %>%
      dplyr::mutate(is_observed =Measure[model_name == ratio] == "observed") %>%
      dplyr::mutate(Measure = if_else(is_observed,Measure, "predicted"))

    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name)) +
      geom_bar(position = "dodge", stat = "identity", col="black") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T)
  } else {
    df = results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::mutate(y = mean(time), sd =sd(time)) %>%
      dplyr::select(n_genes, n_cells, n_cell_types, model_name, y, sd,Measure) %>%
      dplyr::distinct()

    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name, ymin=y-sd, ymax=y+sd)) +
      geom_bar(position = position_dodge(), stat = "identity", col="black") +
      geom_errorbar(col="black", width = .2, position = position_dodge(.99), linetype = "solid") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T, independent = "y")
  }

  y_name = ifelse(isFALSE(ratio), "Runtime (seconds)", paste0("Speedup (vs. ", ratio, ")"))

  p +
    theme_bw() +
    labs(
      x = "Number of Genes",
      y = y_name,
      color = "Model - Device",
      fill = "Model - Device"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    guides(
      fill = guide_legend(order = 1),
      col = guide_legend(order = 1),
      linetype = guide_legend(order = 2, override.aes = list(fill = rep(NA, length(unique(df$Measure)))))
    )
}

memory_comparison = function(results, ratio=FALSE) {

  if (!isFALSE(ratio)) {
    df = results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name,Measure) %>%
      dplyr::summarise(y = mean(memory)) %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types) %>%
      dplyr::mutate(y = y / y[model_name == ratio]) %>%
      dplyr::group_by(n_genes, n_cells) %>%
      dplyr::mutate(is_observed =Measure[model_name == ratio] == "observed") %>%
      dplyr::mutate(Measure = if_else(is_observed,Measure, "predicted"))

    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name)) +
      geom_bar(position = "dodge", stat = "identity", col="black") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T)
  } else {
    df = results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::mutate(y = mean(memory), sd =sd(memory)) %>%
      dplyr::select(n_genes, n_cells, n_cell_types, model_name, y, sd,Measure) %>%
      dplyr::distinct()

    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name, ymin=y-sd, ymax=y+sd)) +
      geom_bar(position = position_dodge(), stat = "identity", col="black") +
      geom_errorbar(col="black", width = .2, position = position_dodge(.99), linetype = "solid") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T, independent = "y")
  }

  y_name = ifelse(isFALSE(ratio), "Memory (GB)", paste0("Relative memory (vs. ", ratio, ")"))

  p +
    theme_bw() +
    labs(
      x = "Number of Genes",
      y = y_name,
      color = "Model - Device",
      fill = "Model - Device"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    guides(
      fill = guide_legend(order = 1),
      col = guide_legend(order = 1)
    )
}

time_per_gene_comparison = function(results) {

  results %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time / n_genes)) %>%
    ggplot(mapping = aes(x = as.factor(n_genes), y = y, fill = model_name)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
    theme_bw() +
    labs(x = "N genes", y = "Time per gene (s)", col = "Model", fill="Model") +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors)
}

plot_time_and_memory_comparison = function(results, n_extrapolation = 3, ncols=2) {
  d <- results %>%
    dplyr::group_by(model_name, n_genes, n_cells) %>%
    dplyr::summarise(memory = mean(memory), time = mean(time)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(size = n_genes * n_cells)

  # Extraolate unkwonw sizes
  sizes = unique(d$size)
  models <- unique(d$model_name)

  d <- lapply(models, function(m) {
    dd <- d %>% dplyr::filter(model_name == m) %>% dplyr::arrange(-size)

    x = c(log10(dd$size[1:n_extrapolation]))
    y_time = c(log10(dd$time[1:n_extrapolation]))
    y_memory = c(log10(dd$memory[1:n_extrapolation]))

    time_lm <- lm(y_time ~ x)
    memory_lm <- lm(y_memory ~ x)

    ddd = dd %>%
      ungroup() %>%
      dplyr::filter(size %in% sizes) %>%
      dplyr::select(model_name, time, memory, size) %>%
      dplyr::mutate(is_extrapolated = "FALSE")

    sizes_to_extrapolate <- sizes[!(sizes %in% ddd$size)]
    if (length(sizes_to_extrapolate) > 0) {
      s <- sizes_to_extrapolate[1]

      dd_extr <- lapply(sizes_to_extrapolate, function(s) {
        time_pred = 10^predict(time_lm, newdata = data.frame(x = c(log10(s))))
        mem_pred = 10^predict(memory_lm, newdata = data.frame(x = c(log10(s))))
        dplyr::tibble(model_name = m, time = time_pred, memory = mem_pred, size = s, is_extrapolated = "TRUE")
      }) %>% do.call("bind_rows", .)

      ddd <- dplyr::bind_rows(ddd, dd_extr)
    }

    if (any(ddd$is_extrapolated == "TRUE")) {
      fake_row = ddd %>% dplyr::filter(is_extrapolated == "FALSE") %>% dplyr::filter(size == max(size))
      fake_row$is_extrapolated <- "TRUE"
      ddd <- dplyr::bind_rows(ddd, fake_row)
    }

    ddd
  }) %>% do.call("bind_rows", .)

  # Compute Ratios
  d <- d %>%
    dplyr::group_by(size) %>%
    dplyr::mutate(time_ratio = time / time[model_name == "devil (GPU)"]) %>%
    dplyr::mutate(memory_ratio = memory / memory[model_name == "devil (GPU)"]) %>%
    dplyr::mutate(is_extrapolated = ifelse(is_extrapolated == "FALSE", "Observed", "Extrapolated")) %>%
    dplyr::mutate(is_extrapolated = factor(is_extrapolated, levels = c("Observed", "Extrapolated")))


  p = d %>%
    `colnames<-`(c("model_name", "Time (s)", "Memory(GB)", "size", "is_extrapolated", "Speed ratio", "Memory consumption ratio")) %>%
    tidyr::pivot_longer(!c(model_name, is_extrapolated, size)) %>%
    dplyr::mutate(name = factor(name, levels = c("Time (s)", "Memory(GB)", "Speed ratio", "Memory consumption ratio"))) %>%
    ggplot(mapping = aes(x = size, y= value, col = model_name)) +
    geom_point(aes(shape = is_extrapolated), size = 3) +
    geom_line() +
    scale_x_continuous(transform = "log10") +
    scale_y_continuous(transform = "log10") +
    labs(x = "Dataset size", y = "", col = "Model", shape="Measurement") +
    scale_size_continuous(guide = "none") +
    theme_bw() +
    scale_color_manual(values = method_colors) +
    facet_wrap(name~., scales = "free", ncol = ncols, strip.position = "left")

  return(p)

  # p_time = d %>%
  #   ggplot(mapping = aes(x = size, y= time, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Time (s)", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_time_ratio = d %>%
  #   ggplot(mapping = aes(x = size, y= time_ratio, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Time ratio", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_mem <- d %>%
  #   ggplot(mapping = aes(x = size, y= memory * 1e-9, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Memory (GB)", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_mem_ratio <- d %>%
  #   ggplot(mapping = aes(x=size, y= memory_ratio, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Memory ratio", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # list(
  #   p_time=p_time,
  #   p_time_ratio=p_time_ratio,
  #   p_mem=p_mem,
  #   p_mem_ratio=p_mem_ratio
  # )
}

plot_correlations_single <- function(fits_folder) {
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(tibble)
  require(patchwork)

  fits <- list.files(fits_folder, full.names = TRUE)

  devil.res     <- readRDS(fits[grepl("cpu_devil_",     fits)])
  gpu.devil.res <- readRDS(fits[grepl("gpu_devil_",     fits)])
  glm.res       <- readRDS(fits[grepl("cpu_glmGam",     fits)])

  # Helper to reshape comparisons
  make_df <- function(x1, x2, val1, val2, metric) {
    tibble(
      x = x1[[metric]],
      y = x2[[metric]],
      metric = metric,
      comp = paste0(val1, " vs ", val2)
    )
  }

  # Build master dataframe
  df_all <- bind_rows(
    make_df(gpu.devil.res, devil.res, "devil – A100", "devil – CPU", "lfc"),
    make_df(devil.res, glm.res, "devil – CPU", "glmGamPoi – CPU", "lfc"),
    make_df(gpu.devil.res, devil.res, "devil – A100", "devil – CPU", "theta"),
    make_df(devil.res, glm.res, "devil – CPU", "glmGamPoi – CPU", "theta")
  ) %>%
    mutate(
      metric = factor(metric, levels = c("lfc", "theta"), labels = c("LFC", "Theta")),
      comp = factor(comp)
    )

  # Single combined plot
  p <- df_all %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(alpha = 0.35, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
    ggpubr::stat_cor(
      method = "pearson",
      size = 3,
      label.x.npc = "left",
      label.y.npc = "top"
    ) +
    facet_grid(metric ~ comp, scales = "free", switch = "both") +
    #coord_equal() +
    theme_bw() +
    labs(x = "Estimate (model 1)", y = "Estimate (model 2)")

  p
}

plot_correlations_by_model <- function(fits_folder) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  require(tibble)

  fits <- list.files(fits_folder, full.names = TRUE)

  devil.res     <- readRDS(fits[grepl("cpu_devil_", fits)])
  gpu.devil.res <- readRDS(fits[grepl("gpu_devil_", fits)])
  glm.res       <- readRDS(fits[grepl("cpu_glmGam", fits)])

  # Build long data frame
  df_lfc <- tibble(
    devil_cpu = devil.res$lfc,
    `devil - a100`   = gpu.devil.res$lfc,
    `glmGamPoi - cpu` = glm.res$lfc
  ) %>%
    pivot_longer(
      cols = c(`devil - a100`, `glmGamPoi - cpu`),
      names_to = "model",
      values_to = "estimate"
    ) %>%
    mutate(metric = "LFC")

  df_theta <- tibble(
    devil_cpu = devil.res$theta,
    `devil - a100`   = gpu.devil.res$theta,
    `glmGamPoi - cpu` = glm.res$theta
  ) %>%
    pivot_longer(
      cols = c(`devil - a100`, `glmGamPoi - cpu`),
      names_to = "model",
      values_to = "estimate"
    ) %>%
    mutate(metric = "Theta")

  df_pvalue <- tibble(
    devil_cpu = devil.res$adj_pval,
    `devil - a100`   = gpu.devil.res$adj_pval,
    `glmGamPoi - cpu` = glm.res$adj_pval
  ) %>%
    pivot_longer(
      cols = c(`devil - a100`, `glmGamPoi - cpu`),
      names_to = "model",
      values_to = "estimate"
    ) %>%
    mutate(metric = "p-adj")

  df_all <- bind_rows(df_lfc, df_theta, df_pvalue) %>%
    mutate(
      metric = factor(metric, levels = c("LFC", "Theta", "p-adj")),
      model  = factor(model, levels = c("devil - a100", "glmGamPoi - cpu"))
    ) %>%
    filter(
      between(devil_cpu,  -10, 10),
      between(estimate, -10, 10)
    )

  # Plot
  p <- df_all %>%
    ggplot(aes(x = devil_cpu, y = estimate, colour = model)) +
    geom_point(size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
    ggpubr::stat_cor(
      aes(colour = model),
      method = "pearson"
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free") +
    scale_colour_manual(
      name = "Model",
      values = method_colors
    ) +
    labs(
      x = expression(paste("devil-CPU estimate")),
      y = "Other model estimate"
    ) +
    theme_bw()

  p
}




# plot_correlations <- function(fits_folder) {
#   fits <- list.files(fits_folder, full.names = T)
#
#   devil.res <- readRDS(fits[grepl("/cpu_devil_", fits)])
#   gpu.devil.res <- readRDS(fits[grepl("/gpu_devil_", fits)])
#   glm.res <- readRDS(fits[grepl("/cpu_glmGam", fits)])
#
#
#   p1 <- dplyr::bind_rows(
#     dplyr::tibble(
#       x = gpu.devil.res$lfc,
#       x_name = "devil - a100",
#       y = devil.res$lfc,
#       y_name = "devil - cpu"
#     ),
#     dplyr::tibble(
#       x = devil.res$lfc,
#       x_name = "devil",
#       y = glm.res$lfc,
#       y_name = "glmGamPoi"
#     )
#   ) %>%
#     dplyr::group_by(x_name, y_name) %>%
#     dplyr::filter(abs(x) < 10 & abs(y) < 10) %>%
#     ggplot(mapping = aes(x=x, y=y)) +
#     geom_point() +
#     ggpubr::stat_cor() +
#     facet_grid(x_name~y_name, switch = "both") +
#     theme_bw() +
#     labs(x = bquote(LFC[1]), y=bquote(LFC[2]))
#
#   p2 <- dplyr::bind_rows(
#     dplyr::tibble(
#       x = gpu.devil.res$theta,
#       x_name = "devil - gpu",
#       y = devil.res$theta,
#       y_name = "devil - cpu"
#     ),
#     dplyr::tibble(
#       x = devil.res$theta,
#       x_name = "devil",
#       y = glm.res$theta,
#       y_name = "glmGamPoi"
#     )
#   ) %>%
#     dplyr::group_by(x_name, y_name) %>%
#     ggplot(mapping = aes(x=x, y=y)) +
#     geom_point() +
#     ggpubr::stat_cor() +
#     facet_grid(x_name~y_name, switch = "both") +
#     theme_bw() +
#     labs(x = bquote(theta[1]), y=bquote(theta[2]))
#
#   list(lfc=p1, theta=p2)
# }

plot_correlations <- function(fits_folder) {
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(tibble)

  fits <- list.files(fits_folder, full.names = TRUE)

  devil.res     <- readRDS(fits[grepl("cpu_devil_",     fits)])
  gpu.devil.res <- readRDS(fits[grepl("gpu_devil_",     fits)])
  glm.res       <- readRDS(fits[grepl("cpu_glmGam",     fits)])

  ## small helper to avoid repetition
  nice_corr_plot <- function(df, x_lab, y_lab, clip_range = NULL) {
    if (!is.null(clip_range)) {
      df <- df %>%
        filter(
          x >= clip_range[1], x <= clip_range[2],
          y >= clip_range[1], y <= clip_range[2]
        )
    }

    df %>%
      mutate(
        x_name = factor(x_name),
        y_name = factor(y_name)
      ) %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(alpha = 0.35, size = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
      ggpubr::stat_cor(
        method = "pearson",
        label.x.npc = "left",
        label.y.npc = "top",
        size = 3
      ) +
      facet_grid(x_name ~ y_name, switch = "both") +
      coord_equal() +
      labs(x = x_lab, y = y_lab) +
      theme_bw(base_size = 11) +
      theme(
        strip.placement   = "outside",
        strip.background  = element_rect(fill = "grey95", colour = NA),
        panel.grid        = element_blank(),
        axis.title        = element_text(face = "bold"),
        axis.text         = element_text(colour = "black"),
        panel.border      = element_rect(colour = "black"),
        plot.margin       = margin(5.5, 5.5, 5.5, 5.5),
        strip.text.x      = element_text(face = "bold"),
        strip.text.y      = element_text(face = "bold")
      )
  }

  ## LFC plot
  df_lfc <- bind_rows(
    tibble(
      x = gpu.devil.res$lfc,
      x_name = "devil – A100",
      y = devil.res$lfc,
      y_name = "devil – CPU"
    ),
    tibble(
      x = devil.res$lfc,
      x_name = "devil – CPU",
      y = glm.res$lfc,
      y_name = "glmGamPoi – CPU"
    )
  )

  p1 <- nice_corr_plot(
    df_lfc,
    x_lab = bquote(LFC[fit[1]]),
    y_lab = bquote(LFC[fit[2]]),
    clip_range = c(-10, 10)
  )

  ## theta plot
  df_theta <- bind_rows(
    tibble(
      x = gpu.devil.res$theta,
      x_name = "devil – A100",
      y = devil.res$theta,
      y_name = "devil – CPU"
    ),
    tibble(
      x = devil.res$theta,
      x_name = "devil – CPU",
      y = glm.res$theta,
      y_name = "glmGamPoi – CPU"
    )
  )

  p2 <- nice_corr_plot(
    df_theta,
    x_lab = bquote(theta[fit[1]]),
    y_lab = bquote(theta[fit[2]])
    # no clipping here; add clip_range if you want
  )

  list(lfc = p1, theta = p2)
}


# fits_folder <- "results/baronPancreas/fits/"
plot_upset <- function(fits_folder, lfc_cut, pval_cut) {
  fits <- list.files(fits_folder, full.names = T)

  l = list()
  fits = fits[grepl("cpu_devil_", fits) | grepl("gpu_devil_", fits) | grepl("glmGamPoi", fits)]
  f = fits[1]
  res <- lapply(fits, function(f) {
    info = unlist(strsplit(unlist(strsplit(f, "fits//"))[2], "_"))
    name = paste(info[1], info[2])

    de_genes = readRDS(f) %>%
      dplyr::mutate(name = paste0("Gene ", row_number())) %>%
      dplyr::filter(abs(lfc) >= lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::pull(name)
    dplyr::tibble(algorithm = name, de_genes = de_genes)
  }) %>% do.call("bind_rows", .)

  res = res %>%
    dplyr::mutate(algorithm = dplyr::recode(algorithm,
                                             "cpu devil" = "devil - cpu",
                                             "cpu glmGamPoi" = "glmGamPoi - cpu",
                                            "gpu devil" = "devil - a100",
                                             "a100" = "devil - a100",
                                             "h100" = "devil - h100"))
  res %>%
    group_by(de_genes) %>%
    summarize(algorithm = list(algorithm)) %>%
    ggplot(aes(x = algorithm)) +
    geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
    scale_y_continuous(breaks = NULL, lim = c(0, nrow(res) / length(fits)), name = "") +
    scale_x_upset() +
    theme_bw() +
    labs(x = "")
}

plot_venn = function(fits_folder, lfc_cut, pval_cut) {
  fits <- list.files(fits_folder, full.names = T)

  fits = fits[grepl("cpu_devil_", fits) | grepl("gpu_devil_", fits) | grepl("glmGamPoi", fits)]
  res <- lapply(fits, function(f) {
    info = unlist(strsplit(unlist(strsplit(f, "fits//"))[2], "_"))
    name = paste(info[1], info[2])

    de_genes = readRDS(f) %>%
      dplyr::mutate(name = paste0("Gene ", row_number())) %>%
      dplyr::filter(abs(lfc) >= lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::pull(name)
    #dplyr::tibble(algorithm = name, de_genes = de_genes)
    de_genes
  })

  fit_names = lapply(fits, function(f) {
    if (grepl("cpu_devil", f)) return("devil - cpu")
    if (grepl("cpu_glm", f)) return("glmGamPoi - cpu")
    if (grepl("gpu_devil", f)) return("devil - a100")
  }) %>% unlist()

  names(res) = fit_names
  library(ggVennDiagram)
  ggVennDiagram(res, set_color = method_colors[fit_names], label_alpha = 0) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    scale_fill_gradient(low="white",high = "gray80") +
    theme(legend.position = "none")
}

plot_volc = function (
    devil.res,
    lfc_cut = 1,
    pval_cut = 0.05,
    labels = TRUE,
    colors = c("gray", "forestgreen", "steelblue", "indianred"),
    color_alpha = 0.7,
    point_size = 1,
    center = TRUE,
    title = "Volcano plot") {

  if (sum(is.na(devil.res))) {
    message("Warning: some of the reults are unrealiable (i.e. contains NaN)\n Those genes will not be displayed")
    devil.res <- stats::na.omit(devil.res)
  }

  d <- devil.res %>%
    dplyr::mutate(pval_filter = .data$adj_pval <= pval_cut, lfc_filter = abs(.data$lfc) >= lfc_cut) %>%
    dplyr::mutate(class = dplyr::if_else(
      .data$pval_filter & .data$lfc_filter,
      "p-value and lfc",
      dplyr::if_else(
        .data$pval_filter,
        "p-value",
        dplyr::if_else(
          .data$lfc_filter,
          "lfc",
          "non-significant"
          )
        )
      )
      ) %>% dplyr::mutate(class = factor(class, levels = c("non-significant", "lfc", "p-value", "p-value and lfc"))) %>%
    dplyr::mutate(label = dplyr::if_else(class == "p-value and lfc", .data$name, NA))
  d

  if (sum(d$adj_pval == 0) > 0) {
    message(paste0(sum(d$adj_pval == 0), " genes have adjusted p-value equal to 0, will be set to the minimum non-zero value"))
    d$adj_pval[d$adj_pval == 0] <- min(d$adj_pval[d$adj_pval != 0])
  }

  p <- d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$lfc, y = -log10(.data$adj_pval), col = .data$class, label = .data$label)) +
    ggplot2::geom_point(alpha = color_alpha, size = point_size) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "") + ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(pval_cut),
                                                                   linetype = "dashed") + ggplot2::ggtitle(title) + ggplot2::theme(legend.position = "bottom")
  if (center) {
    p <- p + ggplot2::xlim(c(-max(abs(d$lfc)), max(abs(d$lfc))))
  }

  if (labels) {
    p = p + ggrepel::geom_label_repel(max.overlaps = Inf, col = "black")
  }
  p
}

plot_large_test = function(dataset_name, add_competitors = FALSE) {
  results_folder <- file.path("results", dataset_name, "large")
  res_paths = list.files(results_folder)
  p = res_paths[1]
  all_res = read.delim(file.path(results_folder, p), sep = ",") %>%
    dplyr::rename(n_genes=X.genes, n_cells=cells, time=timing_gpu_sec) %>%
    dplyr::mutate(model_name = "devil - h100")

  if (add_competitors) {
    competitors_time = get_time_results(file.path("results/", dataset_name))
    competitors_time = competitors_time %>%
      dplyr::filter(n_genes %in% all_res$n_genes) %>%
      dplyr::group_by(model_name, n_genes, n_cells) %>%
      dplyr::summarise(time = mean(time)) %>%
      dplyr::filter(!grepl("a100", model_name))

    all_res = dplyr::bind_rows(all_res, competitors_time) %>%
      dplyr::group_by(model_name, n_genes, n_cells) %>%
      dplyr::summarise(time = mean(time))

    all_res %>%
      ggplot(mapping = aes(x = n_cells, y=time, col=model_name)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors) +
      #scale_x_continuous(transform = "log10") +
      #scale_y_continuous(transform = "log10") +
      labs(x = "Number of cells", y = "Runtime (seconds)", col = "Model - device")
  } else {
    all_res %>%
      dplyr::mutate(n_cells = factor(n_cells, levels = sort(n_cells))) %>%
      ggplot(mapping = aes(x=n_cells, y=time, fill=model_name)) +
      geom_col(col="black") +
      scale_fill_manual(values = method_colors) +
      ggh4x::facet_nested(~"Number of genes"+n_genes, scales = "free_y", shrink = T, independent = "y") +
      theme_bw() +
      theme_bw() +
      labs(
        x = "Number of Cells",
        y = "Runtime (seconds)",
        color = "Model - Device",
        fill = "Model - Device"
      ) +
      theme(
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors)
  }
}


# Plot overdispersion comparison
plot_overdispersion_comparison = function(dataset_name, ratio) {
  no_disp_path = file.path("results", dataset_name, "no_overdispersion")
  device_paths = list.files(no_disp_path)

  df1 = lapply(device_paths, function(device) {
    lapply(list.files(file.path(no_disp_path, device, dataset_name)), function(p){
      fp = file.path(no_disp_path, device, dataset_name, p)
      if (file.exists(fp) && !dir.exists(fp)) {
        info <- unlist(strsplit(p, "_"))
        res = readRDS(fp)
        dplyr::tibble(
          model_name = paste(info[1], info[2]),
          n_genes = info[3],
          n_cells = info[5],
          n_cell_types = info[7],
          time = as.numeric(unlist(res$time), units = "secs")#,
          #memory = res$mem_alloc
        ) %>% dplyr::mutate(model_name = ifelse(grepl("gpu" ,model_name), device, model_name))
      }
    }) %>% do.call("bind_rows", .) %>%
      dplyr::mutate(model_name = paste0("devil - ", model_name)) %>%
      dplyr::mutate(Overdispersion = "w.o")
  }) %>% do.call("bind_rows", .)
  rm(no_disp_path, device_paths)

  disp_path = file.path("results", dataset_name)
  device_paths = c("a100", "h100")
  df2 = lapply(device_paths, function(device) {
    lapply(list.files(file.path(disp_path, device)), function(p){
      fp = file.path(disp_path, device, p)
      if (file.exists(fp) && !dir.exists(fp)) {
        info <- unlist(strsplit(p, "_"))
        res = readRDS(fp)
        dplyr::tibble(
          model_name = paste(info[1], info[2]),
          n_genes = info[3],
          n_cells = info[5],
          n_cell_types = info[7],
          time = as.numeric(unlist(res$time), units = "secs")#,
          #memory = res$mem_alloc
        ) %>% dplyr::mutate(model_name = ifelse(grepl("gpu" ,model_name), device, model_name))
      }
    }) %>% do.call("bind_rows", .) %>%
      dplyr::mutate(model_name = paste0("devil - ", model_name)) %>%
      dplyr::mutate(Overdispersion = "w")
  }) %>% do.call("bind_rows", .)
  rm(disp_path, device_paths)

  df = dplyr::bind_rows(df1, df2)
  rm(df1, df2)

  df = df %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name, Overdispersion) %>%
    dplyr::mutate(y = mean(time), sd =sd(time)) %>%
    dplyr::select(n_genes, n_cells, n_cell_types, model_name, y, sd, Overdispersion) %>%
    dplyr::distinct() %>%
    dplyr::mutate(n_cells = factor(as.numeric(n_cells), levels = unique(sort(as.numeric(df$n_cells)))))

  if (ratio) {
    df = df %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::mutate(y_diff = y[Overdispersion == "w"] - y) %>%
      dplyr::mutate(y = y[Overdispersion == "w"] / y)
  }


  if (!ratio) {
    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name, ymin=y-sd, ymax=y+sd, linetype = Overdispersion)) +
      geom_bar(position = position_dodge(), stat = "identity", col="black") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T, independent = "y") +
      geom_errorbar(col="black", width = .2, position = position_dodge(.99), show.legend = F)
  } else {
    p = df %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name, ymin=y-sd, ymax=y+sd, linetype = Overdispersion)) +
      geom_bar(position = position_dodge(), stat = "identity", col="black") +
      ggh4x::facet_nested(~"Number of cells"+n_cells, scales = "free_y", shrink = T)
  }

  y_name = ifelse(isFALSE(ratio), "Runtime (seconds)", paste0("Speedup (vs. with Overdispersion)"))

  p +
    theme_bw() +
    labs(
      x = "Number of Genes",
      y = y_name,
      color = "Model - Device",
      fill = "Model - Device"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    guides(
      fill = guide_legend(order = 1),
      col = guide_legend(order = 1),
      linetype = guide_legend(order = 2, override.aes = list(fill = rep(NA, length(unique(df$Overdispersion)))))
    )
}
