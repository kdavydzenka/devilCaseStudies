rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source('scripts/utils.R')
set.seed(123456)

print("DEVIL")

MIN_ITER = 1
N_CELL_TYPES = 2

data = prep_MacaqueBrain_data(N_CELL_TYPES = N_CELL_TYPES)

cnt <- data$cnt
design_matrix <- data$design_matrix
clusters = as.numeric(as.factor(data$clusters))

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

print("Macaque brain dataset")
print(paste0(".  n genes = ", n.genes))
print(paste0(".  n cells = ", n.cells))

# Single test
# Save beta and overdispersion for a specific object
n_genes <- 1000
n_cells <- 1e6

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.cells = n_cells, n.sub.genes = n_genes, clusters = clusters)

print(paste0("N genes = ", n_genes))
print(paste0("N cells = ", n_cells))

c = as.matrix(input$c)
d = input$d
cls = input$clusters

print("Devil fitting...")

s <- Sys.time()
fit.devil <- devil::fit_devil(
  c,
  d,
  overdispersion = "MOM",
  size_factors = NULL,
  verbose = T,
  parallel.cores = 1,
  offset = 1e-6,
  init_overdispersion = NULL, 
  init_beta_rough = TRUE,
  max_iter = 500
)
e <- Sys.time()
print(e-s)

devil.res <- devil::test_de(fit.devil, contrast = c(0,1))
devil.final.res <- devil.res %>%
  cbind(fit.devil$beta) %>%
  cbind(dplyr::tibble(theta = fit.devil$overdispersion))
rm(fit.devil, devil.res)

print("Saving results...")

saveRDS(devil.final.res, paste0("results/MacaqueBrain/fits/cpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

# Scaling test
print("Scaling test...")
for (n_genes in c(100, 1000, 5000)) {
  for (n_cells in c(1000, 100000, 1e6)) {

    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))

    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells)
    c = as.matrix(input$c)
    d = input$d

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
      max_iter = 500
    ), min_iterations = MIN_ITER, memory = T)
    b.devil$result <- NULL

    res_name = paste0("results/MacaqueBrain/cpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil, file = res_name)
  }
}

rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source('scripts/utils.R')
set.seed(123456)

print("GLMGAMPOI")

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

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.cells = n_cells, n.sub.genes = n_genes)

print(paste0("N genes = ", n_genes))
print(paste0("N cells = ", n_cells))

c = as.matrix(input$c)
d = input$d

print("glmGamPoi fitting...")

s <- Sys.time()
fit.glm <- glmGamPoi::glm_gp(
  c,
  d,
  overdispersion = T, 
  size_factors = FALSE, 
  offset = 1e-6
)
e <- Sys.time()
print(e-s)

glm.res <- glmGamPoi::test_de(fit.glm, contrast = c(0,1))
glm.final.res <- glm.res %>%
  cbind(fit.glm$Beta) %>%
  cbind(dplyr::tibble(theta = fit.glm$overdispersions))
rm(glm.res, fit.glm)

print("Saving results...")

saveRDS(glm.final.res, paste0("results/MacaqueBrain/fits/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

print("Subsampling test...")

for (n_genes in c(100, 1000, 5000)) {
  for (n_cells in c(1000, 100000, 1e6)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))

    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells)
    c = as.matrix(input$c)
    d = input$d

    print("inference starting glm ...")

    b.glmGamPoi <- bench::mark(
      glmGamPoi::glm_gp(c, d, overdispersion = T, size_factors = FALSE, offset = 1e-6), min_iterations = MIN_ITER, memory = T
    )
    b.glmGamPoi$result <- NULL

    res_name = paste0("results/MacaqueBrain/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.glmGamPoi, file = res_name)
  }
}
