setwd("~/GitHub/devilCaseStudies/multiomics_analysis/")
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", "CNAqc")
sapply(pkgs, require, character.only = TRUE)

cell_group_colors = c(
  "old" = "darkorange",
  "young" = "steelblue"
)

# input UMAPs ####
load("results/metadata_rna_umap.Rdata")
load("results/metadata_atac_umap.Rdata")

d_atac <- metadata_atac %>%
  `rownames<-`(1:nrow(metadata_atac)) %>%
  dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  #dplyr::sample_n(5000) %>%
  dplyr::select(UMAP_1, UMAP_2, cell_type, group) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::filter(UMAP_1 >= -1, UMAP_2 >= -6, UMAP_2 <= 7.5)

d_rna <- metadata_rna %>%
  dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  #dplyr::sample_n(5000) %>%
  dplyr::mutate(group = if_else(age_pop == "young_pop", "young", "old")) %>%
  dplyr::select(umap_1, umap_2, cell_type, group) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::filter(umap_1 >= -5, umap_2 >= 0, umap_2 <= 11, umap_1 <= 1)

colnames(d_atac) <- colnames(d_rna)
d_omics <- rbind(d_atac %>% dplyr::mutate(tech = "ATAC"), d_rna %>% dplyr::mutate(tech = "RNA"))
umaps <- ggplot() +
  geom_point(d_omics %>% dplyr::select(!cell_type), mapping = aes(x=umap_1, y=umap_2), size = .2, col="gainsboro", alpha = .3) +
  geom_point(d_omics, mapping = aes(x=umap_1, y=umap_2, col=group), size = .2, alpha = 1) +
  theme_bw() +
  labs(x="UMAP 1", y="UMAP 2", col="Sample age") +
  #facet_wrap(tech~cell_type, scales = "free") +
  ggh4x::facet_nested_wrap(tech~cell_type, scales = 'free') +
  #ggh4x::facet_nested() +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = cell_group_colors) +
  theme(legend.position = 'bottom')
#umaps
rm(d_atac, metadata_atac, d_omics, d_rna, metadata_rna)


# Read results ####
grange_path <- "results/grange_annot_scADA.RDS"
#grange_path <- "multiomics_analysis/results/grange_annot.RDS"

grange <- readRDS(grange_path)
atac_scaDA_path <- "results/MuscleATAC/scADA_res.RDS"
atac_scaDA <- readRDS(atac_scaDA_path)
atac_scaDA <- cbind(grange, atac_scaDA)
atac_scaDA$adj_pval = p.adjust(atac_scaDA$pval, "BH")

atac_scaDA$annotation %>% unique()
atac_scaDA <- atac_scaDA %>%
  dplyr::filter(annotation != 'Promoter (2-3kb)') %>%
  #dplyr::filter(annotation == 'Promoter (<=1kb)') %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::filter(abs(log2fc) == max(abs(log2fc)))
  #dplyr::filter(abs(distanceToTSS) == min(abs(distanceToTSS)))
  # dplyr::mutate(log2fc = mean(log2fc), pval = mean(pval), FDR = mean(FDR)) %>%
  #dplyr::distinct(SYMBOL, log2fc, pval, FDR)
# atac_scaDA <- atac_scaDA %>%
#   dplyr::filter(distanceToTSS <= 100, distanceToTSS >= -1000)
atac_scaDA$log2fc <- (-1 * atac_scaDA$log2fc)
atac_scaDA$geneID = atac_scaDA$SYMBOL

# atac_scaDA <- atac_scaDA %>%
#   group_by(geneID) %>%
#   dplyr::filter(abs(distanceToTSS) == min(abs(distanceToTSS)))

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

# volcano plots ####
lfc_cut <- 1.0
lfc_cut_atac <- .5
pval_cut <- .05
pval_cut_atac <- .05
de_gene_colors <- c("Not significant" = "gainsboro", "Down-regulated" = "steelblue", "Up-regulated"="indianred")

devil_d <- rna_devil %>%
  dplyr::mutate(isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(lfc > 0, "Up-regulated", "Down-regulated"))) %>%
  dplyr::mutate(method = "devil")

glm_d <- rna_glm %>%
  dplyr::mutate(adj_pval = if_else(adj_pval == 0, min(adj_pval[adj_pval!=0]), adj_pval)) %>%
  dplyr::mutate(isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(lfc > 0, "Up-regulated", "Down-regulated"))) %>%
  dplyr::mutate(method = "glmGamPoi")

nebula_d <- rna_nebula %>%
  dplyr::mutate(adj_pval = if_else(adj_pval == 0, min(adj_pval[adj_pval!=0]), adj_pval)) %>%
  dplyr::mutate(isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(lfc > 0, "Up-regulated", "Down-regulated"))) %>%
  dplyr::mutate(method = "NEBULA")

atac_d <- atac_scaDA %>%
  dplyr::mutate(adj_pval = FDR, lfc = log2fc) %>%
  dplyr::mutate(adj_pval = if_else(adj_pval == 0, min(atac_scaDA$FDR[atac_scaDA$FDR != 0]), adj_pval)) %>%
  dplyr::mutate(isDE = (abs(lfc) >= lfc_cut_atac) & (adj_pval <= pval_cut_atac)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(lfc > 0, "Up-regulated", "Down-regulated"))) %>%
  dplyr::ungroup() %>% 
  dplyr::select(geneID, pval, adj_pval, lfc, isDE, DEtype) %>%
  dplyr::mutate(method = "ATAC")

#Remove outliers
row.remove.neb <- c("C21orf91", "AL137246.2")
row.remove.devil <- c("CASP4", "KCTD1")
devil_d <- devil_d[!(devil_d$geneID %in% row.remove.devil), ]
nebula_d <- nebula_d[!(nebula_d$geneID %in% row.remove.neb), ]

gene_sign <- c("FOSL2", "FOSB", "STAT3", "TNNT2", "ID1", "DCLK1", "SAA2-SAA4",
               "SLC2A4", "PPARA", "MYH1", "MYH2", "MYH4")


data <- rbind(devil_d, glm_d, nebula_d, atac_d)
volcanos <-  data %>% 
  ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=DEtype)) +
  geom_point(size=.8) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  facet_wrap(~method, scales = "free", nrow=2) +
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
volcanos

# Add labels for significant genes #
volcanos <- volcanos +
  ggrepel::geom_label_repel(
    data =  data %>% filter(geneID %in% gene_sign),
    aes(label = geneID),
    color = "black",
    size = 1.2,
    box.padding = 0.1,
    point.padding = 0.1,
    segment.color = 'black'
  )
volcanos

ggsave("plot/volcanos.pdf", plot = volcanos, dpi = 300, width = 6, height = 6)

rm(devil_d, glm_d, nebula_d, atac_d)

# Filter differential expressed genes ####
atac_deg <- atac_scaDA %>% dplyr::filter(FDR < pval_cut_atac, abs(log2fc) > lfc_cut_atac)
rna_deg_devil <- rna_devil %>% dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)
rna_deg_glm <- rna_glm %>% dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)
rna_deg_nebula <- rna_nebula %>% dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)

devil <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_devil$geneID)
glm <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_glm$geneID)
nebula <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_nebula$geneID)
all_lists <- list(devil=devil, glmGamPoi=glm, NEBULA=nebula)
all_genes <- unique(c(atac_deg$geneID, rna_deg_devil$geneID, rna_deg_glm$geneID, rna_deg_nebula$geneID))

# Our plot ####
library(biomaRt)
coords <- CNAqc::chr_coordinates_GRCh38
from = coords$from
names(from) = coords$chr
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
tool_deg <- rna_deg_devil

DE_genes <- gene_positions <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = tool_deg$geneID,
  mart = mart
) %>% dplyr::mutate(pos = (start_position + end_position) / 2) %>%
  dplyr::filter(chromosome_name %in% c(1:22, 'X', 'Y'))

DA_genes <- gene_positions <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = atac_deg$geneID,
  mart = mart
) %>% dplyr::mutate(pos = (start_position + end_position) / 2) %>%
  dplyr::filter(chromosome_name %in% c(1:22, 'X', 'Y'))

genes_db <- rbind(DE_genes, DA_genes) %>%
  distinct() %>%
  dplyr::mutate(y1 = dplyr::if_else(hgnc_symbol %in% DE_genes$hgnc_symbol, "DE", "non DE")) %>%
  dplyr::mutate(y2 = dplyr::if_else(hgnc_symbol %in% DA_genes$hgnc_symbol, "DA", "non DA")) %>%
  dplyr::mutate(chr = paste0("chr", chromosome_name)) %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(x = start_position + from[chr])

dd <- dplyr::bind_rows(
  coords %>% dplyr::mutate(class = 'non DE'),
  coords %>% dplyr::mutate(class = 'DA'),
  coords %>% dplyr::mutate(class = 'DE'),
  coords %>% dplyr::mutate(class = 'non DA')
) %>% dplyr::mutate(class = factor(class, levels = c("non DE", "DA", "DE", "non DA")))

ggplot() +
  geom_segment(genes_db, mapping=aes(x=x, xend=x, y=y1, yend=y2, col=y1)) +
  geom_vline(xintercept = coords$from, color = "darkslategray", linetype = "dashed") +
  geom_vline(xintercept = coords$to, color = "darkslategray", linetype = "dashed") +
  geom_line(dd, mapping = aes(x = from, y=class, col=class)) +
  scale_x_continuous(breaks = coords$centromerStart, labels = stringr::str_replace_all(coords$chr, pattern = "chr", "")) +
  theme_bw() +
  theme(legend.position = 'none')

genes_db %>%
  group_by(y1, y2) %>%
  dplyr::summarise(n = n())

# UpSet and Venn plots ####
lt = list(
  ATAC = atac_deg$geneID,
  devil = rna_deg_devil$geneID,
  NEBULA = rna_deg_nebula$geneID,
  glmGamPoi = rna_deg_glm$geneID
)
upset_plot <- UpSetR::upset(UpSetR:::fromList(lt), order.by = "freq", text.scale = 1.8, sets.bar.color = c("#EAB578", "#099668", "purple3", "#E4A6A7"))

pdf("plot/upset.pdf", width = 10, height = 5)
upset_plot %>% print()
dev.off()

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  "steelblue",
  "devil" = "#099668"
)

i <- 1

venn_plots <- lapply(1:length(all_lists), function(i) {
  algo <- names(all_lists)[i]
  ggvenn::ggvenn(all_lists[[i]], c("snATAC", "snRNA"), show_percentage = FALSE, set_name_size = 4, fill_color = unname(c("purple3", method_colors[algo])), fill_alpha = .75) +
    ggplot2::labs(title=algo) +
    theme(plot.title=element_text(hjust=0.5, vjust=0.5))
})

#venn_plots[[1]]

venns <- patchwork::wrap_plots(venn_plots, ncol = 3)
ggsave("plot/venns.pdf", plot = venns, dpi = 300, width = 10, height = 5)

# Correlation plots ####
d_corr_devil <- rna_deg_devil %>%
  dplyr::left_join(atac_deg, by="geneID") %>%
  dplyr::filter(!is.na(log2fc)) %>%
  dplyr::mutate(method = 'devil')

d_corr_glm <- rna_deg_glm %>%
  dplyr::left_join(atac_deg, by="geneID") %>%
  dplyr::filter(!is.na(log2fc)) %>%
  dplyr::mutate(method = 'glmGamPoi')

d_corr_nebula <- rna_deg_nebula %>%
  dplyr::left_join(atac_deg, by="geneID") %>%
  dplyr::filter(!is.na(log2fc)) %>%
  dplyr::mutate(method = 'NEBULA')

corr_plot <- rbind(d_corr_devil, d_corr_glm, d_corr_nebula) %>%
  ggplot2::ggplot(mapping = aes(x = lfc, y = log2fc)) +
  geom_point(shape = 21, fill = 'black', size = 1) +
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE) +
  smplot2::sm_statCorr(fit.params = list(color = "indianred"), separate_by = "\n", corr_method = 'spearman') +
  #smplot2::sm_statCorr(corr_method = "spearman", fit.params = list(color = method_colors["devil"])) +
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_bw() +
  facet_wrap(. ~ method, scales = "free")
corr_plot
dev.off()

ggsave("plot/corr_plot.pdf", dpi = 300, width = 14, height = 5, plot = corr_plot)

corr_plot / venns

#saveRDS(d_corr_devil, file = "results/d_corr_devil.RDS")
#saveRDS(d_corr_glm, file = "results/d_corr_glm.RDS")
#saveRDS(d_corr_nebula, file = "results/d_corr_nebula.RDS")

# Gene set Enrichement analysis ####
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)
deg_list <- list(
  "devil" = d_corr_devil,
  "NEBULA" = d_corr_nebula,
  "glmGamPoi" = d_corr_glm
)

n <- names(deg_list)[1]
res_gse_list <- list()
plots <- lapply(names(deg_list), function(n) {
  print(n)
  overlap_genes <- deg_list[[n]]

  # Convert Gene symbols to EntrezID #
  hs <- org.Hs.eg.db
  my_symbols <- overlap_genes$geneID
  gene_list <- AnnotationDbi::select(hs,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")

  # gene_list <- na.omit(gene_list)
  # gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
  # gene_list <- gene_list[!is.na(gene_list$SYMBOL),]
  # gene_list <- gene_list[gene_list$SYMBOL %in% overlap_genes$geneID,]
  gene_list_rank <- as.vector(overlap_genes$lfc)
  names(gene_list_rank) <- gene_list$SYMBOL
  gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

  gseGO <- clusterProfiler::gseGO(
    gene_list_rank,
    ont = "All",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    eps = 0,
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = "BH",
    pvalueCutoff = .05
  )

  res_gse <- gseGO@result
  res_gse <- res_gse %>%
    mutate(gene_clusters = case_when(NES > 0  ~ 'up-regulated', NES < 0  ~ 'down-regulated'))

  res_gse$log_padjust <- -log10(res_gse$p.adjust)
  res_gse_list[[n]] <<- res_gse
  
})

saveRDS(res_gse_list, file = "results/res_gse_list.RDS")

#GO_pathways <- c(
#"actin filament-based movement",
#"actin-mediated cell contraction",
#"muscle contraction",
#"muscle system process",
#"structural constituent of muscle",
#"response to oxygen-containing compound",
#"response to fibroblast growth factor",
#"ERK1 and ERK2 cascade",
#"response to hydrogen peroxide",
#"striated muscle cell proliferation",
#"cell-cell adhesion",
#"inflammatory response",
#"negative regulation of apoptotic process",
#"negative regulation of cellular metabolic process",
#"cellular oxidant detoxification"
#)

# Select enriched pathways #

devil <- res_gse_list['devil'] %>% as.data.frame()
glm <- res_gse_list['glmGamPoi'] %>% as.data.frame()
nebula <- res_gse_list['NEBULA'] %>% as.data.frame()

GO_pathways <- c("immune system process",
                 "response to oxygen-containing compound",
                 "actin filament binding",
                 "actin cytoskeleton",
                 "actin filament-based process",
                 "regulation of angiogenesis",
                 "molecular function inhibitor activity ",
                 "reactive oxygen species metabolic process",
                 "leukocyte activation",
                 "negative regulation of cellular metabolic process",
                 "negative regulation of programmed cell death",
                 "enzyme inhibitor activity",
                 "regulation of cytokine production",
                 "negative regulation of cellular biosynthetic process",
                 "defence responce")


res_gse_devil <- res_gse_list['devil']$devil %>% 
  dplyr::mutate(method = "devil") %>% 
  filter(Description %in% GO_pathways)

res_gse_glm <- res_gse_list['glmGamPoi']$glmGamPoi %>% 
  dplyr::mutate(method = "glmGamPoi") %>% 
  filter(Description %in% GO_pathways)

res_gse_nebula <- res_gse_list['NEBULA']$NEBULA %>% 
  dplyr::mutate(method = "NEBULA") %>% 
  filter(Description %in% GO_pathways)

res_gse <- rbind(res_gse_devil, res_gse_glm, res_gse_nebula)
res_gse$method <- as.factor(res_gse$method)

# Barplot #
#pbar <- ggplot(res_gse %>% filter(Description %in% GO_pathways), aes(x = log_padjust, y = Description, fill=method)) +
  #geom_bar(stat="identity", position="dodge", width = 0.50) +
  #labs(x = "-Log10(pvalue)", y = "GO Pathways", fill="") +
  #theme_bw() +
  #scale_fill_manual(values = method_colors) +
  #theme(text = element_text(size = 12)) +
  #theme(legend.position = 'bottom')
#pbar

pbar <- ggplot(res_gse, aes(x = log_padjust, y = Description, fill=method)) +
  geom_bar(stat="identity", position="dodge", width = 0.50) +
  labs(x = "-Log10(pvalue)", y = "GO Pathways", fill="") +
  theme_bw() +
  scale_fill_manual(values = method_colors) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = 'bottom')
pbar

# Main figure ####
design <- "
AAABBB
AAABBB
CCCCCC
CCCCCC
LLLLLL
LLLLLL
LLLLLL
"

main_fig <- wrap_plots(
  A=umaps,
  B=volcanos,
  C=corr_plot,
  #F=wrap_elements(upset_plot),
  L=wrap_elements(pbar),
  design = design
) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 12))
ggsave("plot/main_fig.png", dpi = 400, width = 8.3, height = 11.7, plot = main_fig)


