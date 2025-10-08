library(scaDA)
pkgs <- c("ChIPseeker","TxDb.Hsapiens.UCSC.hg38.knownGene", "ggplot2", "dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula")
sapply(pkgs, require, character.only = TRUE)
data_path <- "datasets/All_snATAC_95021.RDS"

atac <- readRDS(data_path)
counts <- atac@assays@data@listData[["PeakMatrix"]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
grange_annot <- ChIPseeker::annotatePeak(
  peak = atac@rowRanges,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  level = "gene",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                "Downstream", "Intergenic"),
  annoDb = "org.Hs.eg.db",
  addFlankGeneInfo = FALSE,
  flankDistance = 5000,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  overlap = "TSS",
  verbose = TRUE,
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
) %>% as.data.frame()

metadata <- atac@colData %>% as.data.frame()
tissue <- "muscle"

cell_filter <- (metadata$group %in% c("young", "old") & metadata$Annotation %in% c("Type I", "Type II"))


genes_of_interest <- c()
glm_rna <- readRDS("results/MuscleRNA/glmGamPoi_rna.RDS")
genes_of_interest <- c(genes_of_interest, glm_rna$name)
devil_rna <- readRDS("results/MuscleRNA/devil_rna.RDS")
genes_of_interest <- c(genes_of_interest, devil_rna$name)
nebula_rna <- readRDS("results/MuscleRNA/nebula_rna.RDS")
genes_of_interest <- c(genes_of_interest, nebula_rna$name)
genes_of_interest <- unique(genes_of_interest)

range_filter <- (grange_annot$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")) & (grange_annot$SYMBOL %in% genes_of_interest)

metadata <- metadata[cell_filter,]
counts <- counts[,cell_filter]
counts <- counts[range_filter,]
grange_annot <- grange_annot[range_filter,]


#
# coldata = metadata$group
# scaDA.obj <- scaDAdatasetFromMatrix(count = as.matrix(counts), colData = data.frame(coldata))
#
# scaDA.obj <- estParams(scaDA.obj, group.1 = "old", group.2 = "young")
# scaDA.obj <- shrinkDisp(scaDA.obj)
# scaDA.obj <- optParams(scaDA.obj)
# res = scaDA.obj@result
#
# saveRDS(res, "results/scADA_res.glm")
# saveRDS(grange_annot, "results/grange_annot")

counts_avg <- lapply(unique(grange_annot$SYMBOL), function(s) {
  print(s)
  if (length(which(grange_annot$SYMBOL == s)) > 1) {
    rowSums(counts[which(grange_annot$SYMBOL == s),])
  } else {
    counts[which(grange_annot$SYMBOL == s),]
  }
}) %>% do.call("rbind", .)

counts <- counts_avg
colnames(counts) <- NULL
sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=counts), colData=metadata)
sce.pb <- glmGamPoi::pseudobulk_sce(
  sce.obj,
  group_by=vars(sample, group),
  verbose=FALSE
)


design <- model.matrix(~1+group, data=sce.pb@colData)
edger.obj <- edgeR::DGEList(counts(sce.pb))
edger.obj <- edgeR::estimateDisp(edger.obj, design)
fit <- edgeR::glmQLFit(y=edger.obj, design=design)
test <- edgeR::glmTreat(fit, coef=2)

beta <- test$coefficients[,2]
pval <- test$table[,'PValue']
tval <- qnorm(1-pval/2) * sign(beta)
se <- beta/tval

result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')

saveRDS(result, "results/edgeR_res.rds")
saveRDS(grange_annot, "results/grange_annot_edgeR.rds")
