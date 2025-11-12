# Install/load helpers
# install.packages(c("dplyr","readr","stringr","Matrix","pROC"))
# install.packages("decoupler")  # PanglaoDB via OmniPath
library(dplyr)
library(readr)
library(stringr)
library(Matrix)
library(pROC)
library(decoupleR)

# --- 1) CURATED MARKER SOURCES ----------------------------------------------

## 1A) PanglaoDB (programmatic)
pang <- decoupleR::get_resource(name = "PanglaoDB", organism = "human", license = "academic")
# pang columns typically: gene, cell_type, evidence, ...
pang <- pang %>% transmute(source="PanglaoDB", gene=genesymbol, cell_type=cell_type)

## 1B) LM22 (download once from CIBERSORT portal; plain text matrix) 
# https://cibersort.stanford.edu (LM22 signature)  — contains 547 genes across 22 immune types
# Load and turn nonzero entries into markers
lm22_path <- "data/Datasets/LM22.txt"  # <-- set your local path
lm22 <- read_tsv(lm22_path, show_col_types = FALSE)
lm22_long <- lm22 %>%
  dplyr::rename(Gene.symbol = `Gene symbol`) %>% 
  tidyr::pivot_longer(-Gene.symbol, names_to="cell_type", values_to="w") %>%
  filter(w != 0) %>%
  transmute(source="LM22", gene=Gene.symbol, cell_type=cell_type)

## 1C) HPA immune cell export (optional but recommended)
# Visit HPA immune cell section and export "immune cell enriched/enhanced" genes per cell type.
# https://www.proteinatlas.org/humanproteome/single+cell/immune+cell
# Bind them here as a two-column CSV: gene, cell_type (matching your labels or synonyms)
# hpa <- read_csv("HPA_immune_markers.csv") %>% mutate(source="HPA")

# Combine available sources
markers_all <- bind_rows(
  pang,
  lm22_long
  # , hpa
) %>% distinct()

# --- 2) HARMONIZE CELL TYPE NAMES TO YOUR LABELS -----------------------------

# Your target labels
target_labels <- c(
  "NK", "Naive_CD8_T", "Gamma_delta_T", "Naive_B", "Monocytes",
  "Granulocyte", "Memory_CD4_T", "Memory_CD8_T", "Naive_CD4_T", "Memory_B"
)

# Simple synonym map (extend as needed)
syn_map <- tribble(
  ~from,                      ~to,
  "NK cell",                  "NK",
  "NK cells",                 "NK",
  "T cell NK",                "NK",
  "Naive CD4 T cell",         "Naive_CD4_T",
  "Naive CD8 T cell",         "Naive_CD8_T",
  "Memory CD4 T cell",        "Memory_CD4_T",
  "Memory CD8 T cell",        "Memory_CD8_T",
  "B cell naive",             "Naive_B",
  "Naive B cell",             "Naive_B",
  "B cell memory",            "Memory_B",
  "Memory B cell",            "Memory_B",
  "Monocyte",                 "Monocytes",
  "Classical monocyte",       "Monocytes",
  "Non-classical monocyte",   "Monocytes",
  "Granulocyte",              "Granulocyte",
  "Neutrophil",               "Granulocyte",
  "Gamma delta T cell",       "Gamma_delta_T",
  "GdT-cell",                 "Gamma_delta_T"
)

markers_h <- markers_all %>%
  mutate(cell_type = str_replace_all(cell_type, "_", " "),
         cell_type = str_trim(cell_type)) %>%
  left_join(syn_map, by = c("cell_type" = "from")) %>%
  mutate(cell_type = coalesce(to, cell_type)) %>%
  select(-to) %>%
  filter(cell_type %in% target_labels) %>%
  distinct(gene, cell_type, source)

# --- 3) SCORE/RANK MARKERS IN YOUR DATA -------------------------------------

# 3A) helper: compute AUC (one-vs-rest) per gene per type using log1p CPM
libsize <- colSums(cnt_mat)
cpm <- t(t(cnt_mat) / libsize * 1e6)
expr <- log1p(cpm)

compute_auc <- function(g, ct) {
  y <- as.integer(Labels$cell_type == ct)
  s <- as.numeric(expr[g, ])
  if (length(unique(y)) < 2 || sd(s) == 0) return(NA_real_)
  suppressMessages(pROC::roc(y, s, quiet=TRUE)$auc)
}

# 3B) optionally: incorporate your DE results (mean logFC one-vs-rest), else skip
# If you already computed one-vs-rest logFC per gene per ct into a data.frame `de_long`
# with columns: gene, cell_type, logFC  — join it; otherwise set logFC = NA.
de_long <- tibble(gene=character(), cell_type=character(), logFC=numeric())

# 3C) rank per cell type
ranked <- markers_h %>%
  group_by(cell_type, gene) %>%
  summarise(source_count = n_distinct(source), .groups="drop") %>%
  group_by(cell_type) %>%
  # compute AUCs
  mutate(auc = vapply(gene, compute_auc, numeric(1), ct = unique(cell_type))) %>%
  ungroup() %>%
  left_join(de_long, by = c("gene","cell_type")) %>%
  mutate(mean_logFC = logFC) %>%
  group_by(cell_type) %>%
  # score: prefer multiple-source support + higher AUC + higher logFC
  mutate(score = 2*source_count + scales::rescale(coalesce(auc, 0), to=c(0,1)) +
           0.5*scales::rescale(coalesce(mean_logFC, 0), to=c(0,1))) %>%
  arrange(cell_type, desc(score)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 50) %>%
  select(cell_type, gene, rank, auc, mean_logFC, source_count)

# Save & inspect
write_csv(ranked, "marker_top50_by_celltype.csv")
head(ranked, 20)
