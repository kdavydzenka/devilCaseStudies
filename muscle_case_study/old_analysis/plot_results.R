# Plot results #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

load("~/metadata_atac_umap.Rdata")
load("~/metadata_rna_umap.Rdata")

### UMAP ###
metadata_rna <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]
umap_rna <- metadata_rna %>% dplyr::select(umap_1, umap_2, cell_type)
umap_rna$cell_type <- as.factor(umap_rna$cell_type)

p_umap_rna <- ggplot(umap_rna, aes(x = umap_1, y = umap_2, color = cell_type))+
  geom_point(size=0.2)+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snRNA (~ 212 000 nuclei)",
       color = "cell_type") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_classic() +
  scale_fill_manual("cell_type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)
  
p_umap_rna <- Seurat::LabelClusters(plot = p_umap_rna, id = 'cell_type', color="black")+
  theme(legend.position="none")


umap_atac <- metadata_atac %>% dplyr::select(UMAP_1, UMAP_2, cell_type)
umap_atac$cell_type <- as.factor(umap_atac$cell_type)

p_umap_atac <-ggplot(umap_atac, aes(x = UMAP_1, y = UMAP_2, color = cell_type))+
  geom_point(size=0.2)+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snATAC (~ 90 000 nuclei)")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_classic() +
  scale_fill_manual("cell type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)

p_umap_atac <- Seurat::LabelClusters(plot = p_umap_atac, id = 'cell_type', color="black")+
  theme(legend.position="none")

p_umap_rna + p_umap_atac


### Volcano plot ###
  
# Top cell type & muscle aging specific genes #
data_rna_devil <- "multiomics_analysis/results/MuscleRNA/devil_rna_new.RDS"
data_rna_glm <- "multiomics_analysis/results/MuscleRNA/glmGamPoi_rna_new.RDS"
data_rna_nebula <- "multiomics_analysis/results/MuscleRNA/nebula_rna.RDS"
#data_rna_edge <- "multiomics_analysis/results/MuscleRNA/edge_rna.RDS"
data_atac_scaDA <- "multiomics_analysis/results/atac_nodup_scaDA.RDS"


rna_devil <- readRDS(data_rna_devil)
rna_glm <- readRDS(data_rna_glm)
#rna_edge <- readRDS(data_rna_edge)
rna_nebula <- readRDS(data_rna_nebula)
atac <- readRDS(data_atac_scaDA)

#atac <- atac[!duplicated(atac$geneID), ]
#rna_edge$name <- rownames(rna_edge)

top_genes <- devil_rna %>% 
  dplyr::filter(name %in% c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1")) 
# Remove outliers
row.remove.neb <- c("C21orf91", "AL137246.2")
row.remove.devil <- c("CASP4", "PREPL", "SATB2")
rna_nebula <- rna_nebula[!(rna_nebula$name %in% row.remove.neb), ]
rna_devil <- rna_devil[!(rna_devil$name %in% row.remove.devil), ]

p3 <- EnhancedVolcano::EnhancedVolcano(rna_nebula,
                                       lab = rna_nebula$name,
                                       x = 'lfc',
                                       y = 'adj_pval',
                                       selectLab = c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                                                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1"),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       pCutoff = 0.05,
                                       FCcutoff = 1.0,
                                       pointSize = 1.0,
                                       labSize = 2.0,
                                       labCol = 'black',
                                       labFace = 'bold',
                                       boxedLabels = TRUE,
                                       col=c('gray', 'gray', 'gray', 'purple'),
                                       colAlpha = 4/5,
                                       legendPosition = 'bottom',
                                       legendLabSize = 10,
                                       legendIconSize = 3.0,
                                       subtitle = "",
                                       gridlines.major = F, gridlines.minor = F,
                                       border = 'full', borderWidth = 0.5, borderColour = 'black',
                                       titleLabSize = 10,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

plot1 <- p1 + ggplot2::labs(title="Devil: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot1

plot2 <- p2 + ggplot2::labs(title="glmGamPoi: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot2

plot3 <- p3 + ggplot2::labs(title="Nebula: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot3

plot4 <- p4 + ggplot2::labs(title="scaDA: snATAC") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot4

gridExtra::grid.arrange(plot1, plot2, plot3, nrow = 1)


### P-values and lfc comparison ###

rna_devil <- rna_devil[rna_devil$name %in% rna_nebula$name,]
rna_nebula <- rna_nebula[rna_nebula$name %in% rna_devil$name,]
rna_glm <- rna_glm[rna_glm$name %in% rna_devil$name,]

p_devil <- rna_devil %>% dplyr::select(name,adj_pval) %>% 
  dplyr::rename(geneID=name,pval_devil=adj_pval)

p_nebula <- rna_nebula %>% dplyr::select(name,adj_pval) %>% 
  dplyr::rename(geneID=name,pval_nebula=adj_pval)
p_pval_1 <- dplyr::full_join(p_devil, p_nebula, by = "geneID")

p_glm <- rna_glm %>% dplyr::select(name,adj_pval) %>% 
  dplyr::rename(geneID=name,pval_glm=adj_pval)
p_pval_2 <- dplyr::full_join(p_devil, p_glm, by = "geneID")

pval_rna_1 <- ggplot(p_pval_1, aes(x=-log10(pval_nebula), y=-log10(pval_devil))) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "p-value snRNA")+
  xlab("-log10(adjpval) snRNA from Nebula") +
  ylab ("-log10(adjpval) snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

pval_rna_2 <- ggplot(p_pval_2, aes(x=-log10(pval_glm), y=-log10(pval_devil))) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "p-value snRNA")+
  xlab("-log10(adjpval) snRNA from glmGamPoi") +
  ylab ("-log10(adjpval) snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

lfc_devil <- rna_devil %>% dplyr::select(name,lfc) %>% 
  dplyr::rename(geneID=name,lfc_devil=lfc)

lfc_nebula <- rna_nebula %>% dplyr::select(name,lfc) %>% 
  dplyr::rename(geneID=name,lfc_nebula=lfc)

lfc_glm <- rna_glm %>% dplyr::select(name,lfc) %>% 
  dplyr::rename(geneID=name,lfc_glm=lfc)

p_lfc_1 <- dplyr::full_join(lfc_devil, lfc_nebula, by = "geneID")
p_lfc_2 <- dplyr::full_join(lfc_devil, lfc_glm, by = "geneID")

lfc_rna_1 <- ggplot(p_lfc_1, aes(x=lfc_nebula, y=lfc_devil)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "log2FC snRNA")+
  xlab("log2FC snRNA from Nebula") +
  ylab ("log2FC snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

lfc_rna_2 <- ggplot(p_lfc_2, aes(x=lfc_glm, y=lfc_devil)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "log2FC snRNA")+
  xlab("log2FC snRNA from glmGamPoi") +
  ylab ("log2FC snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

#pval_rna_1+lfc_rna_1
#pval_rna_2+lfc_rna_2

gridExtra::grid.arrange(pval_rna_1,lfc_rna_1, pval_rna_2,lfc_rna_2, nrow = 2)
