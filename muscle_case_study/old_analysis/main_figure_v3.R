setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis/")
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork")
sapply(pkgs, require, character.only = TRUE)

cell_group_colors = c(
  "old" = "skyblue3",
  "young" = "#FF9999"
)


# input UMAPs ####
load("results/metadata_rna_umap.Rdata")
load("results/metadata_atac_umap.Rdata")

d_atac <- metadata_atac %>%
  `rownames<-`(1:nrow(metadata_atac)) %>%
  dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  dplyr::select(UMAP_1, UMAP_2, cell_type, group) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::filter(UMAP_1 >= -1, UMAP_2 >= -6, UMAP_2 <= 7.5)

d_rna <- metadata_rna %>%
  dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  dplyr::mutate(group = if_else(age_pop == "young_pop", "young", "old")) %>%
  dplyr::select(umap_1, umap_2, cell_type, group) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::filter(umap_1 >= -5, umap_2 >= 0, umap_2 <= 11, umap_1 <= 1)

conditions <- c("Myonuclei TI", "Myonuclei TII")
replacement_values <- c("Myonuclei TI/TII", "Myonuclei TI/TII")

d_atac$cell_type <- as.character(d_atac$cell_type)
d_atac$cell_type <- replace(d_atac$cell_type, d_atac$cell_type %in% conditions, replacement_values)

d_rna$cell_type <- as.character(d_rna$cell_type)
d_rna$cell_type <- replace(d_rna$cell_type, d_rna$cell_type %in% conditions, replacement_values)

colnames(d_atac) <- colnames(d_rna)
d_omics <- rbind(d_atac %>% dplyr::mutate(tech = "ATAC"), d_rna %>% dplyr::mutate(tech = "RNA"))
d_omics$cell_type <- as.factor(d_omics$cell_type)

#umaps <- ggplot() +
  #geom_point(d_omics %>% dplyr::select(!cell_type), mapping = aes(x=umap_1, y=umap_2), size = .2, col="gainsboro", alpha = .3) +
  #geom_point(d_omics, mapping = aes(x=umap_1, y=umap_2, col=group), size = .2, alpha = 1) +
  #theme_bw() +
  #labs(x="UMAP 1", y="UMAP 2", col="Sample age") +
  #facet_wrap(~tech) +
  #ggh4x::facet_nested(tech~cell_type) +
  #ggh4x::facet_nested_wrap(tech~cell_type, scales = 'free')+
  #guides(color = guide_legend(override.aes = list(size=2))) +
  #scale_color_manual(values = cell_group_colors) +
  #theme(legend.position = 'bottom')

umaps <- ggplot(d_omics, aes(x = umap_1, y = umap_2, color = group))+
  geom_point(size=0.2)+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color = "Sample age") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_bw() +
  facet_wrap(~tech, nrow=2, scale = "free") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = cell_group_colors) +
  theme(legend.position = 'bottom')

umaps <- Seurat::LabelClusters(plot = umaps, id = 'cell_type', color="black", 
                               repel = T, position = "median", box = T)
umaps

ggsave("plot/umaps.pdf", plot = umaps, dpi = 300, width = 3, height = 6)

rm(d_atac, metadata_atac, d_omics, d_rna, metadata_rna)

#### Read results ####
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

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

# Volcano plots ####
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
    size = 1.5,
    box.padding = 0.1,
    point.padding = 0.1,
    segment.color = 'black'
  )
volcanos

ggsave("plot/volcanos.pdf", plot = volcanos, dpi = 300, width = 8, height = 7)
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


#### Gene set Enrichement analysis ####

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

deg_list <- list(
  "devil" = d_corr_devil,
  "NEBULA" = d_corr_nebula,
  "glmGamPoi" = d_corr_glm
)

#GO_pathways <- c("regulation of cytokine production",
              #"inflammatory response",
              #"negative regulation of apoptotic process",
              #"immune system process",
              #"regulation of actin filament organization",
              #"actin polymerization or depolymerization",
              #"regulation of supramolecular fiber organization",
              #"reactive oxygen species metabolic process",
              #"response to calcium ion",
              #"regulation of interleukin-6 production",
              #"phagocytosis",
              #"response to calcium ion",
              #"response to oxidative stress",
              #"negative regulation of cellular process",
              #"regulation of angiogenesis")

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


n <- names(deg_list)[1]
res_gse_list <- list()

plots_data <- lapply(names(deg_list), function(n) {
  print(n)
  overlap_genes <- deg_list[[n]]

  # Convert Gene symbols to EntrezID #
  hs <- org.Hs.eg.db
  my_symbols <- overlap_genes$geneID
  gene_list <- AnnotationDbi::select(hs,
                                     keys = my_symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
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
    maxGSSize = 500
  )

  # Select enriched pathways #
  res_gse <- gseGO@result
  res_gse <- res_gse %>%
    mutate(gene_clusters = case_when(NES > 0  ~ 'up-regulated', NES < 0  ~ 'down-regulated'))

  res_gse$log_padjust <- -log10(res_gse$p.adjust)
  res_gse_list[[n]] <<- res_gse

})


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

saveRDS(res_gse , file = "results/res_gse.RDS")


# Barplot #
#res_gse <- readRDS("results/res_gse.RDS")
#res_gse$log_padjust <- -log10(res_gse$p.adjust)

pbar <- ggplot(res_gse, aes(x = log_padjust, y = Description, fill=method)) +
  geom_bar(stat="identity", position="dodge", width = 0.50) +
  labs(x = "-Log10(pvalue)", y = "GO Pathways", fill="") +
  theme_bw() +
  scale_fill_manual(values = method_colors) +
  theme(text = element_text(size = 8)) +
  theme(legend.position = 'bottom')
pbar

ggsave("plot/pbar.pdf", dpi = 300, width = 7, height = 5, plot = pbar)

# Main figure ####
design <- "
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
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
  theme(text = element_text(size = 10))
ggsave("plot/main_fig.png", dpi = 400, width = 8.3, height = 11.7, plot = main_fig)
