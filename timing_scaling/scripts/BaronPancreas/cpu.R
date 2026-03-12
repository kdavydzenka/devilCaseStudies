rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source("scripts/utils.R")
set.seed(123456)

MIN_ITER = 1
N_CELL_TYPES = 2

data = prep_data_small()

cnt <- data$cnt
design_matrix <- data$design_matrix
clusters = data$donors

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

# Single full test ####
# Save beta and overdispersion for a specific object
n_genes <- 1000
n_cells <- 4000

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells, clusters = clusters)
c = as.matrix(input$c)
d = input$d
cls = input$clusters

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
  max_iter = 100
)
e <- Sys.time()
print(e-s)

devil.res <- devil::test_de(fit.devil, contrast = c(0,1))
devil.final.res <- devil.res %>%
  cbind(fit.devil$beta) %>%
  cbind(dplyr::tibble(theta = fit.devil$overdispersion))
rm(fit.devil, devil.res)

print("glmGamPoi fitting...")

s <- Sys.time()
fit.glm <- glmGamPoi::glm_gp(
  as.matrix(c),
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

print("Fitting NEBULA")

nebula_grouped = nebula::group_cell(as.matrix(c), id = cls, pred = d)
s <- Sys.time()
if (is.null(nebula_grouped)) {
  fit.nebula = nebula::nebula(count = as.matrix(c), id = cls, pred = d, ncore = 1)
} else {
  fit.nebula = nebula::nebula(nebula_grouped$count, nebula_grouped$id, nebula_grouped$pred, ncore = 1)
}
e <- Sys.time()
print(e-s)

nebula.final.res = dplyr::tibble(
  lfc = fit.nebula$summary$logFC_labelbeta / log(2), 
  "1" = fit.nebula$summary$`logFC_(Intercept)`, 
  "2" = fit.nebula$summary$logFC_labelbeta, 
  theta = fit.nebula$overdispersion$Cell, 
  pval = fit.nebula$summary$p_labelbeta,
  adj_pval = p.adjust(fit.nebula$summary$p_labelbeta, "BH"))

print("Saving results...")

saveRDS(devil.final.res, paste0("results/baronPancreas/fits/cpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))
saveRDS(glm.final.res, paste0("results/baronPancreas/fits/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))
saveRDS(nebula.final.res, paste0("results/baronPancreas/fits/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

# Scaling test ####
print("Scaling test...")
n_genes = 1000
n_cells = 4000
for (n_genes in c(100, 1000, 10000)) {
  for (n_cells in c(500, 1000, 4000)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))

    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells, clusters = clusters)
    c = as.matrix(input$c)
    d = input$d
    cls = input$clusters
    grouped_nebula = nebula::group_cell(c, id = cls, pred = d)
    if (!is.null(grouped_nebula)) {
      c = grouped_nebula$count
      cls = grouped_nebula$id
      d = grouped_nebula$pred 
    }
    
    print("inference starting glm ...")
    
    b.glmGamPoi <- bench::mark(
      glmGamPoi::glm_gp(c, d, overdispersion = T, size_factors = FALSE, offset = 1e-6), min_iterations = MIN_ITER, memory = T, iterations = MIN_ITER
    )
    b.glmGamPoi$result <- NULL
    
    # f = function() {glmGamPoi::glm_gp(c, d, overdispersion = T, size_factors = FALSE, offset = 1e-6)}
    # res = peakRAM::peakRAM(f)
    # b.glmGamPoi = dplyr::tibble(model = "glmGamPoi", time = res$Elapsed_Time_sec, mem_alloc = res$Peak_RAM_Used_MiB)
    # rm(res)
    
    res_name = paste0("results/baronPancreas/cpu/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.glmGamPoi, file = res_name)

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
    ), min_iterations = MIN_ITER, memory = T, iterations = MIN_ITER)
    b.devil$result <- NULL
    
    # f = function() {devil::fit_devil(
    #   c,
    #   d,
    #   overdispersion = "MOM",
    #   size_factors = NULL,
    #   verbose = F,
    #   parallel.cores = 1,
    #   offset = 1e-6,
    #   init_overdispersion = NULL, 
    #   init_beta_rough = TRUE,
    #   max_iter = 500
    # )}
    # res = peakRAM::peakRAM(f)
    # b.devil = dplyr::tibble(model = "devil", time = res$Elapsed_Time_sec, mem_alloc = res$Peak_RAM_Used_MiB)
    # rm(res)

    res_name = paste0("results/baronPancreas/cpu/cpu_devil_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil, file = res_name)
    
    print("inference starting NEBULA ...")
    b.nebula <- bench::mark(nebula::nebula(
      c,
      cls,
      d,
      ncore = 1
    ), min_iterations = MIN_ITER, memory = T, , iterations = MIN_ITER)
    b.nebula$result <- NULL
    b.nebula$memory <- NULL
    
    # f = function() {nebula::nebula(c, cls, d, ncore = 1)}
    # res = peakRAM::peakRAM(f)
    # b.nebula = dplyr::tibble(model = "nebula", time = res$Elapsed_Time_sec, mem_alloc = res$Peak_RAM_Used_MiB)
    
    res_name = paste0("results/baronPancreas/cpu/cpu_NEBULA_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.nebula, file = res_name)
  }
}

