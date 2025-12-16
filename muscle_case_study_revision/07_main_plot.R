
# go_semantic_plots_clean.R
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(tidytext)
  library(stringr)
})

source("utils/go_and_ora_classes.R")

# ---------------------------
# Settings
# ---------------------------
# NOTE: keep these as file keys for I/O
models_in_files <- c("devil", "glm", "nebula")
conditions_in_files <- c("age_type1", "age_type2", "interaction")
TOP_K <- 10

# Pretty labels for plots
method_levels <- c("devil", "glmGamPoi", "NEBULA")
method_map <- c(
  "devil"     = "devil",
  "glm"       = "glmGamPoi",
  "glmGamPoi" = "glmGamPoi",
  "nebula"    = "NEBULA",
  "NEBULA"    = "NEBULA"
)

condition_levels <- c("age_type1", "age_type2", "interaction")
condition_labs <- c(
  "age_type1"   = "Aging Type I",
  "age_type2"   = "Aging Type II",
  "interaction" = "Interaction term"
)

method_colors <- c(
  "glmGamPoi" = "#EAB578",
  "NEBULA"    = "steelblue",
  "devil"     = "#099668"
)

pretty_term <- function(x, width = 30) stringr::str_wrap(x, width = width)

# ---------------------------
# 1) GSEA lollipop plots (saved as RDS per method)
# ---------------------------
gsea_df <- map_dfr(conditions_in_files, function(cnd) {
  df <- readRDS(
    file.path("results/MuscleRNA/per_contrast_vector_analysis/full",
              cnd, "gsea_GO_list_df_simp.RDS")
  )
  
  df %>%
    group_by(method) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = TOP_K) %>%
    ungroup() %>%
    mutate(
      Class     = go_program_5[Description],
      condition = cnd
    )
}) %>%
  mutate(
    method = dplyr::recode(method, !!!method_map),
    method = factor(method, levels = method_levels),
    condition = factor(condition, levels = condition_levels, labels = condition_labs)
  )

for (m in levels(gsea_df$method)) {
  plot_df_lolli <- gsea_df %>%
    mutate(
      Class = ifelse(is.na(Class), "Other", Class),
      score = -log10(p.adjust),
      term  = pretty_term(Description)
    ) %>%
    filter(method == m) %>%
    group_by(condition) %>%
    mutate(term_f = reorder_within(term, score, condition)) %>%
    ungroup()
  
  p <- ggplot(plot_df_lolli, aes(x = score, y = term_f)) +
    geom_segment(aes(x = 0, xend = score, yend = term_f),
                 linewidth = 0.5, alpha = 0.8) +
    geom_point(aes(size = abs(enrichmentScore), color = Class), alpha = 0.99) +
    scale_y_reordered() +
    facet_wrap(~ condition, nrow = 1, scales = "free_y") +
    scale_size_continuous(range = c(2.2, 7), name = "|Enrichment|") +
    labs(
      x = expression(-log[10](adjP)),
      y = NULL,
      color = "Program"
    ) +
    theme_bw(base_size = 11) +
    ggsci::scale_color_nejm()
  p
  saveRDS(p, paste0("figures/go_plot_", m, ".rds"))
}

# ---------------------------
# 2) GO semantic embedding plots
# ---------------------------
go_df <- tibble()

for (cnd in conditions_in_files) {
  for (m in models_in_files) {
    r <- readRDS(file.path(
      "results/MuscleRNA/per_contrast_vector_analysis/full",
      cnd, paste0("gseGO_", m, ".RDS")
    ))
    
    go_df <- bind_rows(
      go_df,
      r@result %>%
        mutate(
          method    = m,
          condition = cnd,
          Class     = go_program_5[Description]
        )
    )
  }
}

go_df <- go_df %>%
  mutate(
    Class = ifelse(is.na(Class), "Other", Class),
    method = dplyr::recode(method, !!!method_map),
    method = factor(method, levels = method_levels),
    condition = factor(condition, levels = condition_levels, labels = condition_labs)
  ) %>%
  distinct(ID, Description, method, p.adjust, condition, Class)

# ---- semantic similarity matrix (Wang)
hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
go_ids <- unique(go_df$ID)

S <- mgoSim(go_ids, go_ids, semData = hsGO, measure = "Wang", combine = NULL)
S[is.na(S)] <- 0
D <- as.dist(1 - S)

# ---- 2D embedding (MDS)
xy <- cmdscale(D, k = 2) %>% as.data.frame()
colnames(xy) <- c("dim1", "dim2")
xy$ID <- go_ids

plot_df <- go_df %>%
  left_join(xy, by = "ID") %>%
  mutate(size = -log10(p.adjust))

make_mst_edges <- function(df) {
  n <- nrow(df)
  if (n < 3) return(tibble(x=double(), y=double(), xend=double(), yend=double()))
  
  coords <- as.matrix(df[, c("dim1","dim2")])
  distmat <- as.matrix(dist(coords))
  
  g <- graph_from_adjacency_matrix(distmat, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  mst_g <- mst(g, weights = E(g)$weight)
  el <- as.data.frame(as_edgelist(mst_g))
  colnames(el) <- c("i","j")
  el$i <- as.integer(el$i); el$j <- as.integer(el$j)
  
  tibble(
    x    = df$dim1[el$i],
    y    = df$dim2[el$i],
    xend = df$dim1[el$j],
    yend = df$dim2[el$j]
  )
}

edges_df <- plot_df %>%
  group_by(condition, method) %>%
  group_modify(~ make_mst_edges(.x)) %>%
  ungroup()

# ---- Base semantic plot (overlay by condition)
semantic_plot_base <- ggplot(plot_df, aes(dim1, dim2)) +
  geom_segment(
    data = edges_df,
    aes(x=x, y=y, xend=xend, yend=yend, color=method),
    alpha = 0.25, linewidth = 0.6
  ) +
  geom_point(aes(color = method, size = size), alpha = 0.75) +
  facet_grid(~condition) +
  #facet_wrap(~condition, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    x = "GO semantic dimension 1",
    y = "GO semantic dimension 2",
    color = "Method",
    size  = expression(-log[10](adjP))
  )

# ---------------------------
# 3) Labels: curated list (now keyed by pretty condition labels)
# ---------------------------
terms_to_plot_raw <- list(
  "age_type1" = c(
    "myofibril assembly",
    "sarcomere organization",
    "angiogenesis",
    "muscle cell development",
    "muscle filament sliding",
    "macroautophagy",
    "immune response-activating cell surface receptor signaling pathway",
    "respiratory system development",
    "lung development",
    "mitochondrial respiratory chain complex assembly",
    "ATP synthesis coupled electron transport",
    "mitochondrial ATP synthesis coupled electron transport",
    "regulation of cell division",
    "leukocyte migration",
    "regulation of apoptotic signaling pathway",
    "vesicle organization",
    "positive regulation of erythrocyte differentiation"
  ),
  "age_type2" = c(
    "inner ear development",
    "ear development",
    "ribosome biogenesis",
    "T cell activation",
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "tRNA metabolic process",
    "glycoprotein metabolic process",
    "RNA catabolic process",
    "macroautophagy",
    "regulation of apoptotic signaling pathway",
    "intrinsic apoptotic signaling pathway",
    "signal transduction by p53 class mediator"
  ),
  "interaction" = c(
    "ribosomal large subunit biogenesis",
    "positive regulation of long-term synaptic potentiation"
  )
)

# Map to pretty condition labels
terms_to_plot <- list(
  "Aging Type I"        = terms_to_plot_raw[["age_type1"]],
  "Aging Type II"       = terms_to_plot_raw[["age_type2"]],
  "Interaction term"    = terms_to_plot_raw[["interaction"]]
)

# Label df WITH method+condition so facet_grid(method ~ condition) works correctly
label_df <- imap_dfr(terms_to_plot, function(terms, cond_label) {
  plot_df %>%
    dplyr::filter(condition == cond_label, Description %in% terms) %>%
    dplyr::select(Description, condition, Class, dim1, dim2, method) %>%
    dplyr::distinct()
})

# ---------------------------
# 4) Two final semantic plots
# ---------------------------

# A) Overlay (labels once per condition, across methods) — uses distinct labels w/out method
semantic_plot_overlay <- semantic_plot_base +
  ggrepel::geom_text_repel(
    data = label_df %>% dplyr::select(Description, condition, dim1, dim2) %>% dplyr::distinct(),
    aes(label = Description),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  )

# B) Grid (labels only in the correct method x condition panel)
semantic_plot_grid <- ggplot(plot_df, aes(dim1, dim2)) +
  geom_segment(
    data = edges_df,
    aes(x=x, y=y, xend=xend, yend=yend, color=method),
    alpha = 0.25, linewidth = 0.6
  ) +
  geom_point(aes(color = method, size = size), alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = label_df,   # IMPORTANT: label_df includes method + condition
    aes(label = Description),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  facet_grid(method ~ condition) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    x = "GO semantic dimension 1",
    y = "GO semantic dimension 2",
    color = "Method",
    size  = expression(-log[10](adjP))
  )

# ---------------------------
# Save outputs
# ---------------------------
saveRDS(semantic_plot_overlay, "figures/semantic_plot_overlay.rds")
saveRDS(semantic_plot_grid,    "figures/semantic_plot_grid.rds")

# Print to screen (optional)
print(semantic_plot_overlay)
print(semantic_plot_grid)
