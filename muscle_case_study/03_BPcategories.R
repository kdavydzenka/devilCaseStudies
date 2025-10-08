# Divide pathways into BP categories #

biological_processes <- list(
  immune_response_and_defense = c(
    "innate immune response activating cell surface receptor signaling pathway",
    "immune response-regulating cell surface receptor signaling pathway",
    "defense response",
    "inflammatory response",
    "positive regulation of immune system process",
    "regulation of response to external stimulus",
    "B cell activation",
    "cellular response to cytokine stimulus",
    "cytokine production",
    "negative regulation of immune system process",
    "positive regulation of response to external stimulus",
    "interleukin-8 production"
    
  ),
  cellular_processes_and_regulation = c(
    "regulation of cytokinesis",
    "positive regulation of cell cycle process",
    "positive regulation of cytokine production",
    "ncRNA processing",
    "cellular component assembly involved in morphogenesis",
    "regulation of protein secretion",
    "negative regulation of cell cycle phase transition",
    "extrinsic apoptotic signaling pathway",
    "extracellular matrix organization",
    "epithelial cell proliferation",
    "response to hydrogen peroxide",
    "regulation of receptor-mediated endocytosis",
    "regulation of RNA splicing",
    "activation of protein kinase activity"
  ),
  
  signaling = c("DNA damage checkpoint signaling",
                "hippo signaling"),
  
  cellular_adhesion_and_movement = c(
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "myeloid leukocyte migration",
    "regulation of cell adhesion",
    "positive regulation of multicellular organismal process",
    "leukocyte migration",
    "leukocyte cell-cell adhesion",
    "positive regulation of cell migration",
    "regulation of substrate adhesion-dependent cell spreading",
    
  ),
  
  metabolism = c("organic cyclic compound catabolic process",
                 "proteolysis",
                 "protein modification by small protein conjugation or removal",
                 "negative regulation of metabolic process",
                 "macromolecule modification"),
  
  transport = c("intracellular transport",
                "iron ion transport"
  ),
  
  development_and_differentiation = c(
    "osteoblast differentiation",
    "cartilage development"
    
  ),
  muscle_function = c(
    "actin filament-based movement",
    "muscle contraction",
    "actin filament-based process",
    "muscle system process"
  )
)

go_terms_levels <- unname(unlist(biological_processes))

data_join$Biological_process = lapply(data_join$Description, function(desc) {
  true_bp <- lapply(names(biological_processes), function(bp) {
    if (desc %in% biological_processes[[bp]]) {
      return(bp)
    }
  }) %>% unlist()
  true_bp
})
