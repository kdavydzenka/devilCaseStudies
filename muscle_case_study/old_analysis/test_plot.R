setwd("~/GitHub/devilCaseStudies/multiomics_analysis/")
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork")
sapply(pkgs, require, character.only = TRUE)

# read results ####
grange_path <- "results/grange_annot_scADA.RDS"
#grange_path <- "multiomics_analysis/results/grange_annot.RDS"

grange <- readRDS(grange_path)
atac_scaDA_path <- "results/MuscleATAC/scADA_res.RDS"
atac_scaDA <- readRDS(atac_scaDA_path)
atac_scaDA <- cbind(grange, atac_scaDA)
atac_scaDA$adj_pval = p.adjust(atac_scaDA$pval, "BH")

atac_scaDA$annotation %>% unique()
atac_scaDA <- atac_scaDA %>%
  dplyr::filter(annotation == 'Promoter (<=1kb)') %>%
  dplyr::group_by(SYMBOL) %>%
  #dplyr::filter(abs(log2fc) == max(abs(log2fc))) %>% 
  dplyr::filter(adj_pval == min(adj_pval))

atac_scaDA$log2fc <- (-1 * atac_scaDA$log2fc)
atac_scaDA$geneID = atac_scaDA$SYMBOL

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))

RNAs <- list(
  "devil" = rna_devil,
  "glmGamPoi" = rna_glm,
  "NEBULA" = rna_nebula
)
log_intervals <- 10^seq(log10(1e-5), log10(1), length.out = 6)
rrr <- lapply(c('devil', 'glmGamPoi', 'NEBULA'), function(m) {
  print(m)
  rr <- lapply(seq(.0, 2, by = .1), function(lfc_cut) {
    print(lfc_cut)
    lfc_cut_atac <- lfc_cut / 2
    r <- lapply(log_intervals, function(pval_cut) {
      atac_deg <- atac_scaDA %>% dplyr::filter(FDR < pval_cut, abs(log2fc) > lfc_cut_atac)
      rna_deg <- RNAs[[m]] %>% dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut)
      d_corr <- rna_deg %>%
        dplyr::left_join(atac_deg, by="geneID") %>%
        dplyr::filter(!is.na(log2fc)) %>% dplyr::mutate(method = 'devil')
      
      if (nrow(d_corr) > 5) {
        c <- cor.test(d_corr$lfc, d_corr$log2fc)
        dplyr::tibble(method=m, corr = c$estimate, pval = c$p.value, lfc_cut=lfc_cut, pval_cut=pval_cut, n_genes = nrow(d_corr))
      } else {
        dplyr::tibble(method=m, corr = NA, pval = NA, lfc_cut=lfc_cut, pval_cut=pval_cut, n_genes = nrow(d_corr))
      }
    }) %>% do.call(bind_rows, .)
    r
  }) %>% do.call(bind_rows, .)
  rr
}) %>% do.call(bind_rows, .)

rrr %>% 
  dplyr::filter(pval_cut <= .05, pval_cut >= .01) %>% 
  dplyr::filter(lfc_cut <= 1) %>% 
  dplyr::arrange(-corr)

rrr %>% 
  dplyr::filter(pval_cut == .01) %>% 
  dplyr::filter(lfc_cut == 0) %>% 
  dplyr::arrange(-corr)

rrr %>% 
  na.omit() %>% 
  ggplot(mapping = aes(x=pval_cut, y=lfc_cut, fill=corr)) +
  geom_raster() +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  facet_wrap(~method)
