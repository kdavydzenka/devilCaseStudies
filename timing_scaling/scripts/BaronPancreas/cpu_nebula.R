rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source("scripts/utils.R")
set.seed(123456)

MIN_ITER = 3
N_CELL_TYPES = 2

data = prep_data_small()

cnt <- data$cnt
design_matrix <- data$design_matrix
clusters = data$donors

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

# Scaling test
print("Scaling test...")
for (n_genes in c(100, 1000, 10000)) {
  for (n_cells in c(500, 1000, 4000)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))
    
    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells, clusters = clusters)
    c = as.matrix(input$c)
    d = input$d
    
    print("inference starting NEBULA ...")
    
    grouped_nebula = nebula::group_cell(c, id = input$clusters, pred = d)  
    b.nebula <- bench::mark(nebula::nebula(
      grouped_nebula$count,
      grouped_nebula$id,
      grouped_nebula$pred, 
      ncore = 1
    ), min_iterations = MIN_ITER, memory = T)
    b.nebula$result <- NULL
    res_name = paste0("results/baronPancreas/cpu/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.nebula, file = res_name)
  }
}

# Save beta and overdispersion for a specific object
n_genes <- 1000
n_cells <- 4000

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells, clusters = data$donors)
c = as.matrix(input$c)
d = input$d
clusters = input$clusters

print("Fitting NEBULA")

nebula_grouped = nebula::group_cell(as.matrix(c), id = clusters, pred = d)
fit.nebula = nebula::nebula(nebula_grouped$count, nebula_grouped$id, nebula_grouped$pred)

nebula.final.res = dplyr::tibble(
  lfc = fit.nebula$summary$logFC_labelbeta / log(2), 
  "1" = fit.nebula$summary$`logFC_(Intercept)`, 
  "2" = fit.nebula$summary$logFC_labelbeta, 
  theta = fit.nebula$overdispersion$Cell, 
  pval = fit.nebula$summary$p_labelbeta,
  adj_pval = p.adjust(fit.nebula$summary$p_labelbeta, "BH"))

print("Saving results...")

saveRDS(nebula.final.res, paste0("results/baronPancreas/fits/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))
