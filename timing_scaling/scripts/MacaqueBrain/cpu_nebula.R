rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source('scripts/utils.R')
set.seed(123456)

print("NEBULA")

MIN_ITER = 1
N_CELL_TYPES = 2

data = prep_MacaqueBrain_data(N_CELL_TYPES = N_CELL_TYPES)

cnt <- data$cnt
design_matrix <- data$design_matrix

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

print("Macaque brain dataset")
print(paste0(".  n genes = ", n.genes))
print(paste0(".  n cells = ", n.cells))

# Save beta and overdispersion for a specific object
print("Single fit test...")
n_genes <- 1000
n_cells <- 1e6

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.cells = n_cells, n.sub.genes = n_genes, clusters = clusters)
c = as.matrix(input$c)
d = input$d
cls = input$clusters

print(paste0("N genes = ", n_genes))
print(paste0("N cells = ", n_cells))

c = as.matrix(input$c)
d = input$d

print("NEBULA fitting...")

s <- Sys.time()
fit.nebula <- nebula::nebula(
  count = c, 
  id = cls, 
  pred = d, 
  ncore = 1
)
e <- Sys.time()
print(e-s)

nebula.final.res = dplyr::tibble(
  lfc = fit.nebula$summary$`logFC_cell_typeGABAergic neuron` / log(2), 
  "1" = fit.nebula$summary$`logFC_(Intercept)`, 
  "2" = fit.nebula$summary$`logFC_cell_typeGABAergic neuron`, 
  theta = fit.nebula$overdispersion$Cell, 
  pval = fit.nebula$summary$`p_cell_typeGABAergic neuron`,
  adj_pval = p.adjust(fit.nebula$summary$`p_cell_typeGABAergic neuron`, "BH"))
rm(glm.res, fit.nebula)

print("Saving results...")

saveRDS(nebula.final.res, paste0("results/MacaqueBrain/fits/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

print("Subsampling test...")
for (n_genes in c(100, 1000, 5000)) {
  for (n_cells in c(1000, 100000, 1e6)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))
    
    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells, clusters = clusters)
    c = as.matrix(input$c)
    d = input$d
    cls = input$clusters
    
    print("inference starting NEBULA ...")
    
    b.nebula <- bench::mark(
      nebula::nebula(count = c, id = cls, pred = d, ncore = 1), min_iterations = MIN_ITER, memory = T
    )
    b.nebula$result <- NULL
    b.nebula$memory <- NULL
    
    res_name = paste0("results/MacaqueBrain/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.glmGamPoi, file = res_name)
  }
}
