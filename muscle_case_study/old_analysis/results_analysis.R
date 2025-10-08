### Results downstream analysis ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")
#setwd("~/GitHub/devilCaseStudies")
#rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork")
sapply(pkgs, require, character.only = TRUE)


### scaDA results analysis ###

grange_scaDA_path <- "multiomics_analysis/results/grange_annot_scADA.RDS"
grange_scaDA <- readRDS(grange_scaDA_path)

grange <- "multiomics_analysis/results/grange_annot.RDS"
grange <- readRDS(grange)

atac_scaDA_path <- "multiomics_analysis/results/MuscleATAC/scADA_res.RDS"
atac_scaDA <- readRDS(atac_scaDA_path)

rna_devil <- "multiomics_analysis/results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "multiomics_analysis/results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "multiomics_analysis/results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

grange_scaDA <- grange_scaDA %>% dplyr::select(SYMBOL) %>%
  dplyr::rename(geneID=SYMBOL)

atac_scaDA <- cbind(grange_scaDA, atac_scaDA)
atac_scaDA$log2fc <- -1 * atac_scaDA$log2fc

grange$annotation %>% unique()

grange <- grange[ (grange$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)" )) , ]
#grange <- grange[ (grange$annotation %in% c("Promoter (1-2kb)")) , ]
#grange <- grange[ (grange$annotation %in% c("Promoter (2-3kb)")) , ]
atac_scaDA <- atac_scaDA[ (atac_scaDA$geneID %in% grange$SYMBOL) , ]

atac_up <- atac_scaDA %>% dplyr::filter(log2fc > 0) %>%
  dplyr::group_by(geneID) %>%
  dplyr::arrange(log2fc)

atac_up <- atac_up[order(atac_up$log2fc, decreasing = TRUE), ]
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]

atac_down <- atac_scaDA %>% dplyr::filter(log2fc < 0) %>%
  dplyr::group_by(geneID) %>%
  dplyr::arrange(log2fc)

atac_down <- atac_down[order(atac_down$log2fc), ]
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
atac_nodup$log2fc_arr <- abs(atac_nodup$log2fc)
atac_nodup <- atac_nodup[order(atac_nodup$log2fc_arr, decreasing = TRUE),] 
atac_nodup <- atac_nodup[!duplicated(atac_nodup$geneID), ]

saveRDS(atac_nodup, file = "multiomics_analysis/results/atac_nodup_scaDA.RDS")

# Gene selection based on LFC & pvalue cutoff #
quantile <- .5
lfc_cut <- 0.5
pval_cut <- .01

atac_deg <- atac_nodup %>%
  #dplyr::filter(FDR < 0.05, abs(log2fc) > stats::quantile(abs(atac_nodup$log2fc), quantile))
  dplyr::filter(FDR < pval_cut, abs(log2fc) > lfc_cut)

rna_deg_devil <- rna_devil %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_devil$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)

rna_deg_glm <- rna_glm %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_glm$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)

rna_deg_nebula <- rna_nebula %>%
  #dplyr::filter(adj_pval < 0.05, abs(lfc) > stats::quantile(abs(rna_nebula$lfc), quantile))
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)


# Vienn diagram #
devil <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_devil$geneID)
glm <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_glm$geneID)
nebula <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_nebula$geneID)
all_lists <- list(devil=devil, glmGamPoi=glm, NEBULA=nebula)

all_genes <- unique(c(atac_deg$geneID, rna_deg_devil$geneID, rna_deg_glm$geneID, rna_deg_nebula$geneID))
tibble_upset <- lapply(all_genes, function(g) {
  in_devil = g %in% rna_deg_devil$geneID
  dplyr::tibble(
    gene = g,
    devil=g %in% rna_deg_devil$geneID,
    NEBULA=g %in% rna_deg_nebula$geneID,
    glmGamPoi=g %in% rna_deg_glm$geneID
  )
}) %>% do.call("bind_rows", .) %>% as.data.frame()

lt = list(
  ATAC = atac_deg$geneID,
  devil = rna_deg_devil$geneID,
  NEBULA = rna_deg_nebula$geneID,
  glmGamPoi = rna_deg_glm$geneID
)


upset_plot <- UpSetR::upset(UpSetR:::fromList(lt), order.by = "freq", text.scale = 1.8, sets.bar.color = c("#EAB578", "#099668", "purple3", "#E4A6A7"))

pdf("multiomics_analysis/plot/upset.pdf", width = 10, height = 5)
upset_plot %>% print()
dev.off()

#method_colors = c(
  #"glmGamPoi" = "#EAB578",
  #"NEBULA" =  "#E4A6A7",
  #"devil" = "#099668"
#)

method_colors = c(
  "glmGamPoi" = "#EAB578",
  "NEBULA" =  "steelblue2",
  "devil" = "#099668"
)

i <- 1

venn_plots <- lapply(1:length(all_lists), function(i) {
  algo <- names(all_lists)[i]
  ggvenn::ggvenn(all_lists[[i]], c("snATAC", "snRNA"), show_percentage = FALSE, set_name_size = 4, fill_color = unname(c("purple3", method_colors[algo])), fill_alpha = .75) +
    ggplot2::labs(title=algo) +
    theme(plot.title=element_text(hjust=0.5, vjust=0.5))
})

venns <- patchwork::wrap_plots(venn_plots, ncol = 3)
ggsave("multiomics_analysis/plot/venns.pdf", plot = venns, dpi = 300, width = 10, height = 5)

### Filter overlapped genes ###
atac_deg_devil <- atac_deg[ atac_deg$geneID %in% rna_deg_devil$geneID,]
rna_deg_devil <- rna_deg_devil[ rna_deg_devil$geneID %in% atac_deg$geneID,]
atac_deg_devil <- atac_deg_devil %>% dplyr::select(geneID, FDR, log2fc)
colnames(atac_deg_devil) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")
rna_deg_devil <- rna_deg_devil %>% dplyr::select(geneID, adj_pval, lfc)
colnames(rna_deg_devil) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")
overlap_devil <- dplyr::full_join(atac_deg_devil, rna_deg_devil, by = "geneID")

saveRDS(overlap_devil, file = "multiomics_analysis/results/overlap_devil.RDS")

atac_deg_glm <- atac_deg[ atac_deg$geneID %in% rna_deg_glm$geneID,]
rna_deg_glm <- rna_deg_glm[ rna_deg_glm$geneID %in% atac_deg_glm$geneID,]
atac_deg_glm<- atac_deg_glm %>% dplyr::select(geneID, FDR, log2fc)
colnames(atac_deg_glm) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")
rna_deg_glm <- rna_deg_glm %>% dplyr::select(geneID, adj_pval, lfc)
colnames(rna_deg_glm) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")
overlap_glm <- dplyr::full_join(atac_deg_glm, rna_deg_glm, by = "geneID")

saveRDS(overlap_glm, file = "multiomics_analysis/results/overlap_glm.RDS")

atac_deg_nebula <- atac_deg[ atac_deg$geneID %in% rna_deg_nebula$geneID,]
rna_deg_nebula <- rna_deg_nebula[ rna_deg_nebula$geneID %in% atac_deg_nebula$geneID,]
atac_deg_nebula <- atac_deg_nebula %>% dplyr::select(geneID, log2fc, FDR)
colnames(atac_deg_nebula) <- c("geneID", "lfc_snATAC", "adj_pval_snATAC")
rna_deg_nebula <- rna_deg_nebula %>% dplyr::select(geneID, lfc, adj_pval)
colnames(rna_deg_nebula) <- c("geneID", "lfc_snRNA", "adj_pval_snRNA")
overlap_nebula <- dplyr::full_join(atac_deg_nebula, rna_deg_nebula, by = "geneID")

saveRDS(overlap_nebula, file = "multiomics_analysis/results/overlap_nebula.RDS")


### Correlation plot ###

corr1 <- ggplot2::ggplot(mapping = aes(x = overlap_devil$lfc_snRNA, y = overlap_devil$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 1) +
  labs(title = "Devil")+
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE) +
  #smplot2::sm_statCorr(corr_method = "spearman", fit.params = list(color = method_colors["devil"])) +
  smplot2::sm_statCorr(fit.params = list(color = method_colors["devil"])) +
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()

corr2 <- ggplot2::ggplot(mapping = aes(x = overlap_glm$lfc_snRNA, y = overlap_glm$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 1) +
  labs(title = "glmGamPoi")+
  xlab("snRNA -log10(FDR)") +
  ylab ("snATAC -log10(FDR)") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  #smplot2::sm_statCorr(corr_method = "spearman", fit.params = list(color = method_colors["glmGamPoi"])) +
  smplot2::sm_statCorr(fit.params = list(color = method_colors["glmGamPoi"])) +
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()


corr3 <- ggplot2::ggplot(mapping = aes(x = overlap_nebula$lfc_snRNA, y = overlap_nebula$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 1) +
  labs(title = "Nebula")+
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  #smplot2::sm_statCorr(corr_method = "spearman", fit.params = list(color = method_colors["NEBULA"])) +
  smplot2::sm_statCorr(fit.params = list(color = method_colors["NEBULA"])) +
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()

corr_plot <- (corr1 | corr2 | corr3) + plot_layout(axis_titles = 'collect')
ggsave("multiomics_analysis/plot/corr_plot.pdf", dpi = 300, width = 16, height = 8)


### Gene set Enrichement analysis ###

#rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

overlap_genes <- overlap_nebula

# Convert Gene symbols to EntrezID #
hs <- org.Hs.eg.db
my_symbols <- overlap_genes$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
gene_list <- gene_list[gene_list$SYMBOL %in% overlap_genes$geneID,]
gene_list_rank <- as.vector(overlap_genes$lfc)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

gseGO <- clusterProfiler::gseGO(gene_list_rank, ont = "All", OrgDb = org.Hs.eg.db,
                                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05)

# Select enriched pathways #
res_gse_devil <- gseGO@result %>% mutate(method = "devil")
res_gse_glm <- gseGO@result %>% mutate(method = "glmGamPoi")
res_gse_nebula <- gseGO@result %>% mutate(method = "NEBULA")

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
              #"regulation of angiogenesis",
              #"response to oxygen-containing compound")

GO_pathways <- c("immune system process",
                 "response to oxygen-containing compound",
                 "actin filament binding",
                 "actin cytoskeleton",
                 "actin filament-based process",
                 "regulation of angiogenesis",
                 "reactive oxygen species metabolic process",
                 "positive regulation of programmed cell death",
                 "leukocyte activation",
                 "negative regulation of cellular metabolic process")

res_gse_devil <- res_gse_devil %>%
  filter(Description %in% GO_pathways) %>%
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

res_gse_nebula <- res_gse_nebula %>%
  filter(Description %in% GO_pathways) %>%
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

res_gse_glm <- res_gse_glm %>%
  filter(Description %in% GO_pathways) %>%
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

res_gse <- rbind(res_gse_devil, res_gse_glm, res_gse_nebula)
res_gse$method <- as.factor(res_gse$method)

saveRDS(res_gse, file = "results/res_gse_v2.RDS")


### Visualize enrichment results ###

#res_gse$log_padjust <- -log10(res_gse$p.adjust)

#plot1 <- ggpubr::ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    #color = "gene_clusters",
                    #palette = c("blue", "#FC4E07"),
                    #sorting = "descending",
                    #rotate = TRUE,
                    #group = "gene_clusters",
                    #dot.size = "setSize",
                    #add = "segments",
                    #title = "GO enrichment - Devil",
                    #xlab = "GO Pathways",
                    #ylab = "-log10(padjust)",
                    #ggtheme = theme_pubr()
#)
#plot1 + theme(legend.position = "right")+
  #theme_bw()+
  #font("xy.text", size = 12, color = "black", face = "plain")+
  #font("title", size = 10, color = "black", face = "bold")+
  #font("xlab", size = 10)+
  #font("ylab", size = 10)


# Barplot #

p_gse <- ggplot(res_gse, aes(x = log_padjust, y = Description, fill=method)) +
  geom_bar(stat="identity", position="dodge", width = 0.50)+
  labs(x = "-Log10(pvalue)",
       y = "GO Pathways",
       title = "GO analysis")+
  theme(strip.text.x = element_text(size=12, color="black", face="bold.italic"))+
  theme_bw()
p_gse <- p_gse + scale_fill_manual(values=c("#0072B2", "#CC79A7"))

p1 <- p1 + font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 10, color = "black", face = "bold")+
  font("xlab", size = 10)+
  font("ylab", size = 10)
p1


saveRDS(res_gse_devil, file = "multiomics_analysis/res_gse_devil.RDS")
saveRDS(res_gse_glm, file = "multiomics_analysis/res_gse_glm.RDS")
saveRDS(res_gse_nebula, file = "multiomics_analysis/res_gse_nebula.RDS")
