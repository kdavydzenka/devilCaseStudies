
go_classification <- c(
  "oxidative phosphorylation"                                          = "Mitochondrial / Energy metabolism",
  "cellular respiration"                                               = "Mitochondrial / Energy metabolism",
  "ATP synthesis coupled electron transport"                           = "Mitochondrial / Energy metabolism",
  "mitochondrial ATP synthesis coupled electron transport"             = "Mitochondrial / Energy metabolism",
  "respiratory electron transport chain"                               = "Mitochondrial / Energy metabolism",
  "aerobic respiration"                                                = "Mitochondrial / Energy metabolism",
  "energy derivation by oxidation of organic compounds"                = "Mitochondrial / Energy metabolism",
  "aerobic electron transport chain"                                   = "Mitochondrial / Energy metabolism",
  
  "leukocyte migration"                                                = "Immune activation / Inflammation",
  "leukocyte apoptotic process"                                        = "Immune activation / Inflammation",
  "immune response-activating cell surface receptor signaling pathway" = "Immune activation / Inflammation",
  "immune response-activating signaling pathway"                       = "Immune activation / Inflammation",
  "T cell activation"                                                  = "Immune activation / Inflammation",
  "immune response-regulating cell surface receptor signaling pathway" = "Immune activation / Inflammation",
  "immune response-regulating signaling pathway"                       = "Immune activation / Inflammation",
  "activation of immune response"                                      = "Immune activation / Inflammation",
  "positive regulation of B cell activation"                           = "Immune activation / Inflammation",
  
  "regulation of apoptotic signaling pathway"                          = "Apoptosis / Cell death",
  "positive regulation of apoptotic process"                           = "Apoptosis / Cell death",
  "positive regulation of programmed cell death"                       = "Apoptosis / Cell death",
  "regulation of extrinsic apoptotic signaling pathway"                = "Apoptosis / Cell death",
  "negative regulation of apoptotic signaling pathway"                 = "Apoptosis / Cell death",
  "intrinsic apoptotic signaling pathway"                              = "Apoptosis / Cell death",
  
  "response to endoplasmic reticulum stress"                           = "Protein homeostasis / Stress response",
  "protein targeting"                                                  = "Protein homeostasis / Stress response",
  "vesicle organization"                                               = "Protein homeostasis / Stress response",
  
  "myofibril assembly"                                                 = "Muscle structural organization",
  
  "homophilic cell adhesion via plasma membrane adhesion molecules"    = "Cell adhesion / Communication",
  "cell-cell adhesion via plasma-membrane adhesion molecules"          = "Cell adhesion / Communication",
  
  "regulation of cell division"                                        = "Cell cycle / Proliferation"
)


ora_classification <- c(
  # Muscle structure, contraction, sarcomere, actin
  "muscle contraction"                                 = "Muscle structural organization / contraction",
  "muscle system process"                              = "Muscle structural organization / contraction",
  "actin-mediated cell contraction"                    = "Muscle structural organization / contraction",
  "skeletal system development"                        = "Muscle structural organization / contraction",
  "striated muscle cell differentiation"               = "Muscle structural organization / contraction",
  "myofibril assembly"                                 = "Muscle structural organization / contraction",
  "actin filament-based movement"                      = "Muscle structural organization / contraction",
  "muscle cell differentiation"                        = "Muscle structural organization / contraction",
  "sarcomere organization"                             = "Muscle structural organization / contraction",
  "actomyosin structure organization"                  = "Muscle structural organization / contraction",
  "striated muscle cell development"                   = "Muscle structural organization / contraction",
  "muscle cell development"                            = "Muscle structural organization / contraction",
  
  # Immune activation & inflammation
  "acute-phase response"                               = "Immune activation / inflammation",
  "humoral immune response"                            = "Immune activation / inflammation",
  "leukocyte mediated immunity"                        = "Immune activation / inflammation",
  "lymphocyte mediated immunity"                       = "Immune activation / inflammation",
  
  # Apoptosis & immune cell death regulation
  "lymphocyte apoptotic process"                       = "Apoptosis / cell death",
  "leukocyte apoptotic process"                        = "Apoptosis / cell death",
  "cell killing"                                       = "Apoptosis / cell death",
  "negative regulation of leukocyte apoptotic process" = "Apoptosis / cell death",
  "regulation of lymphocyte apoptotic process"         = "Apoptosis / cell death",
  
  # Connective tissue development / ECM changes
  "connective tissue development"                      = "Connective tissue / ECM organization"
)



### Results downstream analysis ###
# 
# rm(list = ls())
# pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
#           "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr',
#           "enrichplot", "clusterProfiler", "data.table", "reactome.db", "fgsea", "org.Hs.eg.db")
# sapply(pkgs, require, character.only = TRUE)
# #set.seed(1234)
# source("utils/utils_analysis.R")
# 
# methods <- c("devil", "glmGamPoi", "nebula")
# conditions <- c("age_only", "age_type1", "age_type2", "interaction")
# 
# c = conditions[1]
# TERMS = c()
# for (c in conditions) {
#   res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/", c)  
#   d1 = readRDS(file.path("results/MuscleRNA/per_contrast_vector_analysis/", c, "gsea_GO_list_df.RDS"))
#   d2 = readRDS(file.path("results/MuscleRNA/per_contrast_vector_analysis/",c,"/gsea_GO_list_df_simp.RDS"))
#   
#   terms = dplyr::bind_rows(d1, d2) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::select(Description) %>% dplyr::pull(Description) %>% unique()
#   
#   TERMS = c(TERMS, terms)
#   
# }
# 
# c = conditions[4]
# TERMS = c()
# for (c in conditions) {
#   res.dir = file.path("results/MuscleRNA/per_contrast_vector_analysis/", c)  
#   d = readRDS(file.path("results/MuscleRNA/per_contrast_vector_analysis/", c, "ORA_list_df.RDS"))
#   
#   terms = dplyr::bind_rows(d) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::select(Description) %>% dplyr::pull(Description) %>% unique()
#   
#   TERMS = c(TERMS, terms)
#   
# }
# 
# classify_go_term <- function(term) {
#   case_when(
#     # ---------- VERY SPECIFIC NEW BITS ----------
#     # Ferroptosis: non-apoptotic cell death
#     str_detect(term, regex("ferroptosis", TRUE)) ~
#       "Cell death – ferroptosis",
#     
#     # Carbohydrate / glucose / glycogen metabolism & homeostasis
#     str_detect(term, regex(
#       "glucose metabolic process|glucose catabolic process|gluconeogenesis|glycogen metabolic process|glycogen catabolic process|carbohydrate (meta|catabo|biosynthetic)|monosaccharide|hexose|polysaccharide|glucan|energy reserve metabolic process|carbohydrate homeostasis|glucose homeostasis",
#       TRUE
#     )) ~ "Metabolism – carbohydrate / glucose",
#     
#     # Lipids, membrane lipids, glyco-lipids
#     str_detect(term, regex(
#       "glycolipid|membrane lipid|phosphatidylserine|liposaccharide|glycosyl compound|phosphatidylinositol|prenylation|protein prenylation",
#       TRUE
#     )) ~ "Metabolism – lipid / membrane / glyco",
#     
#     # ---------- IMMUNE / INFLAMMATION ----------
#     str_detect(term, regex(
#       "T cell|B cell|lymphocyte|humoral immune response|type 2 immune response|immune effector process|leukocyte mediated immunity|lymphocyte mediated immunity",
#       TRUE
#     )) ~ "Immune – adaptive / lymphoid",
#     str_detect(term, regex(
#       "myeloid cell differentiation|reactive oxygen species metabolic process",
#       TRUE
#     )) ~ "Immune – innate / myeloid (or ROS-linked)",
#     
#     str_detect(term, regex(
#       "response to molecule of bacterial origin|response to lipopolysaccharide|cellular response to lipopolysaccharide|response to bacterium|cellular response to biotic stimulus|response to xenobiotic stimulus",
#       TRUE
#     )) ~ "Immune – host defence (bacteria / LPS / xenobiotic)",
#     str_detect(term, regex(
#       "interleukin-6 production|interleukin-8 production|tumor necrosis factor( |$)|TNF|cytokine production|type 2 immune response",
#       TRUE
#     )) ~ "Immune – cytokine signalling",
#     
#     str_detect(term, regex(
#       "inflammatory response|acute inflammatory response|acute-phase response|wound healing|response to wounding",
#       TRUE
#     )) ~ "Stress / inflammation / tissue damage",
#     
#     # ---------- APOPTOSIS & CELL DEATH ----------
#     str_detect(term, regex(
#       "apoptotic signaling pathway|apoptotic process|leukocyte apoptotic process|lymphocyte apoptotic process",
#       TRUE
#     )) ~ "Cell death – apoptosis",
#     
#     # ---------- DNA DAMAGE / P53 / CHECKPOINTS ----------
#     str_detect(term, regex(
#       "DNA damage response, signal transduction by p53 class mediator",
#       TRUE
#     )) ~ "DNA damage response / p53",
#     
#     # ---------- SIGNAL TRANSDUCTION & KINASES ----------
#     str_detect(term, regex(
#       "p38MAPK cascade|MAPK cascade|ERK1 and ERK2 cascade|phosphatidylinositol 3-kinase/protein kinase B",
#       TRUE
#     )) ~ "Signal transduction – MAPK / PI3K-AKT",
#     str_detect(term, regex(
#       "regulation of kinase activity|regulation of protein kinase activity|negative regulation of protein kinase activity",
#       TRUE
#     )) ~ "Signal transduction – kinases (regulation)",
#     str_detect(term, regex(
#       "regulation of transferase activity|negative regulation of catalytic activity",
#       TRUE
#     )) ~ "Enzyme activity regulation – general",
#     
#     # ---------- STRESS RESPONSES & DETOX ----------
#     str_detect(term, regex(
#       "reactive oxygen species|oxidative stress|hydrogen peroxide|cellular oxidant detoxification|detoxification|cellular detoxification|toxic substance",
#       TRUE
#     )) ~ "Stress response – oxidative / detox",
#     
#     str_detect(term, regex(
#       "response to (nutrient|glucose|hexose|monosaccharide|carbohydrate)",
#       TRUE
#     )) ~ "Stress / sensing – nutrient / carbohydrate",
#     
#     str_detect(term, regex(
#       "response to tumor necrosis factor|cellular response to tumor necrosis factor",
#       TRUE
#     )) ~ "Stress / signalling – TNF",
#     
#     str_detect(term, regex(
#       "response to corticosteroid|response to glucocorticoid|response to steroid hormone",
#       TRUE
#     )) ~ "Hormone signalling – steroid / glucocorticoid",
#     
#     # ---------- CELL ADHESION / MIGRATION / MOTILITY ----------
#     str_detect(term, regex(
#       "cell-cell adhesion via plasma-membrane adhesion molecules|homophilic cell adhesion via plasma membrane adhesion molecules|cell-substrate adhesion|leukocyte cell-cell adhesion",
#       TRUE
#     )) ~ "Cell adhesion / cell–cell interaction",
#     
#     str_detect(term, regex(
#       "regulation of cell-cell adhesion|positive regulation of cell adhesion|regulation of leukocyte cell-cell adhesion",
#       TRUE
#     )) ~ "Cell adhesion – regulation",
#     
#     str_detect(term, regex(
#       "cell migration|cell motility|locomotion|endothelial cell migration|positive regulation of cell migration|positive regulation of cell motility|positive regulation of locomotion",
#       TRUE
#     )) ~ "Cell migration / motility",
#     
#     # ---------- MUSCLE / CONTRACTILE / HEART ----------
#     str_detect(term, regex(
#       "muscle system process|muscle contraction|striated muscle contraction|skeletal muscle contraction|cardiac muscle contraction|heart contraction|muscle filament sliding|actin-myosin filament sliding|actomyosin structure organization|actin-mediated cell contraction|actin filament-based movement|regulation of muscle contraction|regulation of skeletal muscle contraction|regulation of heart contraction|regulation of system process",
#       TRUE
#     )) ~ "Muscle / contractile function",
#     
#     str_detect(term, regex(
#       "muscle cell differentiation|muscle cell development|striated muscle cell differentiation|striated muscle cell development|myotube differentiation|muscle tissue development|cardiac muscle tissue growth|muscle cell proliferation|striated muscle cell proliferation|muscle organ development|musculoskeletal movement|multicellular organismal movement",
#       TRUE
#     )) ~ "Muscle – development / growth",
#     
#     str_detect(term, regex(
#       "vascular associated smooth muscle cell proliferation|regulation of vascular associated smooth muscle cell proliferation",
#       TRUE
#     )) ~ "Muscle – vascular smooth muscle",
#     
#     # Heart & circulatory
#     str_detect(term, regex(
#       "heart morphogenesis|heart growth|heart process|cardiac chamber morphogenesis|ventricular septum morphogenesis|endocardial cushion development|endocardial cushion morphogenesis|ventricular cardiac muscle cell action potential|blood circulation|regulation of blood circulation",
#       TRUE
#     )) ~ "Cardiovascular – heart / conduction / circulation",
#     
#     # ---------- DEVELOPMENT, MORPHOGENESIS & ORGANOGENESIS ----------
#     str_detect(term, regex(
#       "morphogenesis of an epithelium|epithelial tube morphogenesis|branching morphogenesis of an epithelial tube|morphogenesis of a branching epithelium|morphogenesis of a branching structure|epithelial cell differentiation|regulation of epithelial cell differentiation|negative regulation of epithelial cell differentiation|epithelial cell proliferation|positive regulation of epithelial cell proliferation|regulation of epithelial cell proliferation",
#       TRUE
#     )) ~ "Development – epithelium / branching",
#     
#     str_detect(term, regex(
#       "renal system development|kidney development|ureteric bud development|ureteric bud morphogenesis|mesonephros development|mesonephric tubule development|mesonephric tubule morphogenesis|nephron development|nephron morphogenesis|nephron epithelium morphogenesis|nephron epithelium development|nephron tubule development|nephron tubule morphogenesis|renal tubule development|renal tubule morphogenesis|kidney epithelium development",
#       TRUE
#     )) ~ "Development – kidney / nephron",
#     
#     str_detect(term, regex(
#       "lung development|respiratory system development|respiratory tube development|lung epithelium development",
#       TRUE
#     )) ~ "Development – lung / respiratory",
#     
#     str_detect(term, regex(
#       "skeletal system development|bone morphogenesis|ossification|cartilage development",
#       TRUE
#     )) ~ "Development – bone / cartilage",
#     
#     str_detect(term, regex(
#       "mammary gland morphogenesis|gland morphogenesis|gland development|connective tissue development|mesenchyme development|appendage development|limb development|skin development|epidermis development|keratinocyte differentiation|odontogenesis",
#       TRUE
#     )) ~ "Development – skin / gland / limb / tooth",
#     
#     str_detect(term, regex(
#       "gliogenesis|regulation of nervous system development",
#       TRUE
#     )) ~ "Development – nervous system / glia",
#     
#     str_detect(term, regex(
#       "female pregnancy|multi-organism reproductive process|multi-multicellular organism process",
#       TRUE
#     )) ~ "Reproduction / pregnancy",
#     
#     str_detect(term, regex(
#       "organ growth|heart growth",
#       TRUE
#     )) ~ "Growth – organ / tissue",
#     
#     str_detect(term, regex(
#       "extracellular structure organization|extracellular matrix organization|external encapsulating structure organization",
#       TRUE
#     )) ~ "ECM / extracellular structure",
#     
#     str_detect(term, regex(
#       "angiogenesis|regulation of angiogenesis|positive regulation of angiogenesis|regulation of vasculature development|positive regulation of vasculature development|vascular process in circulatory system|regulation of vasoconstriction|regulation of vascular endothelial growth factor receptor signaling pathway",
#       TRUE
#     )) ~ "Vasculature / angiogenesis",
#     
#     # ---------- HORMONES, NUCLEOTIDES & MISC METABOLISM ----------
#     str_detect(term, regex(
#       "nucleoside metabolic process|ADP metabolic process|pyridine nucleotide metabolic process|nicotinamide nucleotide metabolic process|pyridine-containing compound metabolic process|purine nucleoside diphosphate metabolic process|purine ribonucleoside diphosphate metabolic process",
#       TRUE
#     )) ~ "Metabolism – nucleotides / energy carriers",
#     
#     str_detect(term, regex(
#       "hormone metabolic process",
#       TRUE
#     )) ~ "Metabolism – hormones",
#     
#     str_detect(term, regex(
#       "hemoglobin metabolic process",
#       TRUE
#     )) ~ "Metabolism – heme / hemoglobin",
#     
#     # ---------- TRANSPORT / ION HOMEOSTASIS ----------
#     str_detect(term, regex(
#       "regulation of calcium ion transport|positive regulation of calcium ion transmembrane transporter activity|regulation of metal ion transport|regulation of monoatomic ion transport",
#       TRUE
#     )) ~ "Ion transport / homeostasis – Ca2+ / metals",
#     
#     str_detect(term, regex(
#       "purine-containing compound transmembrane transport",
#       TRUE
#     )) ~ "Transport – nucleotides",
#     
#     # ---------- REGULATION CATCH-ALLS ----------
#     str_detect(term, regex(
#       "regulation of leukocyte activation|regulation of lymphocyte apoptotic process|regulation of leukocyte apoptotic process|negative regulation of leukocyte apoptotic process|negative regulation of lymphocyte apoptotic process",
#       TRUE
#     )) ~ "Immune – regulation of activation / death",
#     
#     str_detect(term, regex(
#       "regulation of system process|regulation of nervous system development",
#       TRUE
#     )) ~ "Regulation – system-level process",
#     
#     # ---------- DEFAULT ----------
#     TRUE ~ "Other / general"
#   )
# }