
go_classification <- c(
  # ---- GSEA NEBULA ----
  "rRNA metabolic process" = "Protein/RNA synthesis, processing & turnover",
  "ribosome biogenesis" = "Protein/RNA synthesis, processing & turnover",
  "leukocyte apoptotic process" = "Stress response, cell death & microenvironment",
  "RNA catabolic process" = "Protein/RNA synthesis, processing & turnover",
  "T cell activation" = "Immune & inflammatory response",
  "pigmentation" = "Development, differentiation & morphogenesis",
  "leukocyte migration" = "Immune & inflammatory response",
  "mitochondrial ATP synthesis coupled electron transport" = "Metabolism & mitochondrial/energy processes",
  "DNA damage response, signal transduction by p53 class mediator" = "Stress response, cell death & microenvironment",
  "myeloid cell homeostasis" = "Immune & inflammatory response",
  
  # ---- GSEA devil ----
  "homophilic cell adhesion via plasma membrane adhesion molecules" = "Development, differentiation & morphogenesis",
  "immune response-activating cell surface receptor signaling pathway" = "Immune & inflammatory response",
  "tRNA processing" = "Protein/RNA synthesis, processing & turnover",
  "positive regulation of cytokine production" = "Immune & inflammatory response",
  "regulation of apoptotic signaling pathway" = "Stress response, cell death & microenvironment",
  "protein localization to nucleus" = "Protein/RNA synthesis, processing & turnover",
  "response to endoplasmic reticulum stress" = "Stress response, cell death & microenvironment",
  "ribosome biogenesis" = "Protein/RNA synthesis, processing & turnover",
  "T cell activation" = "Immune & inflammatory response",
  "macroautophagy" = "Stress response, cell death & microenvironment",
  
  # ---- GSEA glmGamPoi ----
  "positive regulation of apoptotic process" = "Stress response, cell death & microenvironment",
  "positive regulation of response to external stimulus" = "Stress response, cell death & microenvironment",
  
  # ---- ORA NEBULA ----
  "muscle contraction" = "Muscle structure, development & physiology",
  "muscle system process" = "Muscle structure, development & physiology",
  "actin-mediated cell contraction" = "Muscle structure, development & physiology",
  "actin filament-based movement" = "Muscle structure, development & physiology",
  "muscle filament sliding" = "Muscle structure, development & physiology",
  "actin-myosin filament sliding" = "Muscle structure, development & physiology",
  "striated muscle contraction" = "Muscle structure, development & physiology",
  "skeletal system development" = "Development, differentiation & morphogenesis",
  "response to hypoxia" = "Stress response, cell death & microenvironment",
  "acute-phase response" = "Immune & inflammatory response",
  
  # ---- ORA devil ----
  "sarcomere organization" = "Muscle structure, development & physiology",
  "myofibril assembly" = "Muscle structure, development & physiology",
  "cell-cell adhesion via plasma-membrane adhesion molecules" = "Development, differentiation & morphogenesis",
  "regulation of blood circulation" = "Muscle structure, development & physiology",
  "regulation of system process" = "Development, differentiation & morphogenesis",
  "striated muscle cell differentiation" = "Muscle structure, development & physiology",
  "blood circulation" = "Muscle structure, development & physiology",
  "positive regulation of cell migration" = "Development, differentiation & morphogenesis",
  
  # ---- ORA glmGamPoi ----
  "lymphocyte apoptotic process" = "Stress response, cell death & microenvironment",
  "negative regulation of apoptotic signaling pathway" = "Stress response, cell death & microenvironment",
  "extrinsic apoptotic signaling pathway" = "Stress response, cell death & microenvironment",
  "blood vessel morphogenesis" = "Development, differentiation & morphogenesis",
  "glycolipid biosynthetic process" = "Metabolism & mitochondrial/energy processes",
  "gland morphogenesis" = "Development, differentiation & morphogenesis",
  "regulation of extrinsic apoptotic signaling pathway" = "Stress response, cell death & microenvironment",
  "cell killing" = "Stress response, cell death & microenvironment",
  
  # Metabolism & mitochondrial/energy processes
  "oxidative phosphorylation" = "Metabolism & mitochondrial/energy processes",
  "ATP biosynthetic process" = "Metabolism & mitochondrial/energy processes",
  "mitochondrial respirasome assembly" = "Metabolism & mitochondrial/energy processes",
  "ATP metabolic process" = "Metabolism & mitochondrial/energy processes",
  "glycerophospholipid metabolic process" = "Metabolism & mitochondrial/energy processes",
  
  # Immune & inflammatory response
  "leukocyte differentiation" = "Immune & inflammatory response",
  "immune response-regulating cell surface receptor signaling pathway" = "Immune & inflammatory response",
  "leukocyte cell-cell adhesion" = "Immune & inflammatory response",
  "T-helper 17 cell differentiation" = "Immune & inflammatory response",
  
  # Development, differentiation & morphogenesis
  "epithelial cell proliferation" = "Development, differentiation & morphogenesis",
  "angiogenesis" = "Development, differentiation & morphogenesis",
  "positive regulation of blood vessel endothelial cell migration" = "Development, differentiation & morphogenesis",
  "positive regulation of erythrocyte differentiation" = "Development, differentiation & morphogenesis",
  "myelin maintenance" = "Development, differentiation & morphogenesis",
  
  # Protein/RNA synthesis, processing & turnover
  "protein targeting" = "Protein/RNA synthesis, processing & turnover",
  "vesicle organization" = "Protein/RNA synthesis, processing & turnover",
  
  # Stress response, cell death & microenvironment
  "neuron apoptotic process" = "Stress response, cell death & microenvironment",
  
  # Muscle structure, development & physiology
  "actomyosin structure organization" = "Muscle structure, development & physiology",
  
  "regulation of transferase activity" = "Development, differentiation & morphogenesis",
  "regulation of protein kinase activity" = "Development, differentiation & morphogenesis",
  "regulation of receptor-mediated endocytosis" = "Protein/RNA synthesis, processing & turnover",
  "positive regulation of long-term synaptic potentiation" = "Development, differentiation & morphogenesis",

  # Organ-specific processes
  "renal system process" = "Development, differentiation & morphogenesis",
  
  "endothelial cell apoptotic process"                     = "Stress response, cell death & microenvironment",
  "proton motive force-driven ATP synthesis"               = "Metabolism & mitochondrial/energy processes",
  "response to xenobiotic stimulus"                        = "Stress response, cell death & microenvironment",
  "regulation of leukocyte migration"                      = "Immune & inflammatory response",
  "regulation of cell division"                            = "Development, differentiation & morphogenesis",
  "protein K48-linked ubiquitination"                      = "Protein/RNA synthesis, processing & turnover",
  "negative regulation of protein localization to nucleus" = "Protein/RNA synthesis, processing & turnover",
  "positive regulation of TORC1 signaling"                 = "Metabolism & mitochondrial/energy processes",
  
  "muscle cell differentiation"                           = "Muscle structure, development & physiology",
  "striated muscle cell development"                      = "Muscle structure, development & physiology",
  "striated muscle tissue development"                    = "Muscle structure, development & physiology",
  "striated muscle cell development"                      = "Muscle structure, development & physiology",
  "muscle cell development"                               = "Muscle structure, development & physiology",
  "cellular component assembly involved in morphogenesis" = "Development, differentiation & morphogenesis",
  "cellular anatomical entity morphogenesis"              = "Development, differentiation & morphogenesis",
  "response to molecule of bacterial origin"              = "Immune & inflammatory response",
  "response to bacterium"                                 = "Immune & inflammatory response",
  "leukocyte mediated immunity"                           = "Immune & inflammatory response",
  "humoral immune response"                               = "Immune & inflammatory response",
  "response to lipopolysaccharide"                        = "Immune & inflammatory response",
  "wound healing"                                         = "Immune & inflammatory response",
  "inner ear development"                                 = "Development, differentiation & morphogenesis",
  "tRNA metabolic process" = "Protein/RNA synthesis, processing & turnover",
  "glycoprotein metabolic process" = "Metabolism & mitochondrial/energy processes",
  "positive regulation of programmed cell death" = "Stress response, cell death & microenvironment",
  "ribosomal large subunit biogenesis" = "Protein/RNA synthesis, processing & turnover",
  "muscle tissue development"                     = "Muscle structure, development & physiology",
  "kidney development"                            = "Development, differentiation & morphogenesis",
  "renal system development"                      = "Development, differentiation & morphogenesis",
  "extracellular matrix organization"             = "Development, differentiation & morphogenesis",
  "external encapsulating structure organization" = "Development, differentiation & morphogenesis",
  "extracellular structure organization"          = "Development, differentiation & morphogenesis",
  "muscle tissue morphogenesis"                   = "Muscle structure, development & physiology",
  "sodium ion transport"                          = "Transport, ion homeostasis & vesicles",
  "chemotaxis"                                    = "Immune & inflammatory response",
  "muscle adaptation"                             = "Muscle structure, development & physiology",
  "striated muscle adaptation"                    = "Muscle structure, development & physiology",
  "cardiac muscle contraction"                    = "Muscle structure, development & physiology",
  "chemotaxis"                                    = "Immune & inflammatory response",
  "taxis"                                         = "Immune & inflammatory response",
  "regulation of chemotaxis"                      = "Immune & inflammatory response",
  "muscle tissue development"                     = "Muscle structure, development & physiology",
  "positive regulation of locomotion"             = "Immune & inflammatory response",
  "cell chemotaxis"                               = "Immune & inflammatory response"
)


go_program_5 <- c(
  ## 1. Mitochondrial metabolism & proteostasis
  "mitochondrial ATP synthesis coupled electron transport" = "Metabolism & proteostasis",
  "proton motive force-driven ATP synthesis"               = "Metabolism & proteostasis",
  "mitochondrial respirasome assembly"                     = "Metabolism & proteostasis",
  "glycerophospholipid metabolic process"                  = "Metabolism & proteostasis",
  "rRNA metabolic process"                                 = "Metabolism & proteostasis",
  "ribosome biogenesis"                                    = "Metabolism & proteostasis",
  "ribosomal large subunit biogenesis"                     = "Metabolism & proteostasis",
  "tRNA metabolic process"                                 = "Metabolism & proteostasis",
  "RNA catabolic process"                                  = "Metabolism & proteostasis",
  "protein K48-linked ubiquitination"                      = "Metabolism & proteostasis",
  "glycoprotein metabolic process"                         = "Metabolism & proteostasis",
  
  ## 2. Muscle structure & contractile apparatus
  "myofibril assembly"                                     = "Muscle structure & contractility",
  "muscle filament sliding"                                = "Muscle structure & contractility",
  "actomyosin structure organization"                      = "Muscle structure & contractility",
  
  ## 3. Cell stress, damage response & cell death
  "endothelial cell apoptotic process"                     = "Stress & cell death",
  "neuron apoptotic process"                               = "Stress & cell death",
  "leukocyte apoptotic process"                            = "Stress & cell death",
  "regulation of apoptotic signaling pathway"              = "Stress & cell death",
  "negative regulation of apoptotic signaling pathway"     = "Stress & cell death",
  "positive regulation of programmed cell death"           = "Stress & cell death",
  "DNA damage response, signal transduction by p53 class mediator" =
    "Stress & cell death",
  "response to endoplasmic reticulum stress"               = "Stress & cell death",
  "macroautophagy"                                         = "Stress & cell death",
  
  ## 4. Immune & inflammatory signaling
  "leukocyte migration"                                    = "Immune & inflammatory signaling",
  "regulation of leukocyte migration"                      = "Immune & inflammatory signaling",
  "immune response-activating cell surface receptor signaling pathway" =
    "Immune & inflammatory signaling",
  "T cell activation"                                      = "Immune & inflammatory signaling",
  "positive regulation of cytokine production"             = "Immune & inflammatory signaling",
  "acute-phase response"                                   = "Immune & inflammatory signaling",
  
  ## 5. Growth, signaling & tissue remodeling
  "epithelial cell proliferation"                          = "Growth, signaling & remodeling",
  "regulation of cell division"                            = "Growth, signaling & remodeling",
  "angiogenesis"                                           = "Growth, signaling & remodeling",
  "positive regulation of erythrocyte differentiation"     = "Growth, signaling & remodeling",
  "inner ear development"                                  = "Growth, signaling & remodeling",
  "protein targeting"                                      = "Growth, signaling & remodeling",
  "protein localization to nucleus"                        = "Growth, signaling & remodeling",
  "negative regulation of protein localization to nucleus" = "Growth, signaling & remodeling",
  "regulation of transferase activity"                     = "Growth, signaling & remodeling",
  "positive regulation of TORC1 signaling"                 = "Growth, signaling & remodeling",
  "homophilic cell adhesion via plasma membrane adhesion molecules" =
    "Growth, signaling & remodeling",
  "positive regulation of long-term synaptic potentiation" = "Growth, signaling & remodeling",
  "response to xenobiotic stimulus"                         = "Growth, signaling & remodeling",
  "positive regulation of response to external stimulus"   = "Growth, signaling & remodeling"
)
