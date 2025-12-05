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

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

# Save beta and overdispersion for a specific object
n_genes <- 1000
n_cells <- 4000

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells)
c = as.matrix(input$c)
d = input$d

print("Devil fitting...")

s <- Sys.time()
fit.devil <- devil::fit_devil(
  as.matrix(c),
  d,
  overdispersion = "MOM",
  size_factors = NULL,
  verbose = T,
  parallel.cores = 1,
  offset = 1e-6,
  init_overdispersion = NULL, 
  init_beta_rough = TRUE,
  max_iter = 500,
  CUDA = TRUE,
  batch_size = min(1024L, n_genes)
)
e <- Sys.time()
print(e-s)

devil.res <- devil::test_de(fit.devil, contrast = c(0,1))
devil.final.res <- devil.res %>%
  cbind(fit.devil$beta) %>%
  cbind(dplyr::tibble(theta = fit.devil$overdispersion))
rm(fit.devil, devil.res)

saveRDS(devil.final.res, paste0("results/baronPancreas/fits/gpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

# Scaling test
print("Scaling test...")
for (n_genes in c(100, 1000, 10000)) {
  for (n_cells in c(500, 1000, 4000)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))

    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells)
    c = as.matrix(input$c)
    d = input$d

    print("inference starting devil ...")

    print("inference starting devil ...")

    b.devil <- bench::mark(devil::fit_devil(
      c,
      d,
      overdispersion = "MOM",
      size_factors = NULL,
      verbose = T,
      parallel.cores = 1,
      offset = 1e-6,
      init_overdispersion = NULL, 
      init_beta_rough = TRUE,
      max_iter = 500,
      CUDA = TRUE,
      batch_size = min(100, n_genes)
    ), min_iterations = MIN_ITER, memory = T)
    b.devil$result <- NULL

    res_name = paste0("results/baronPancreas/gpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil, file = res_name)
  }
}
