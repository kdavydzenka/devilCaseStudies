
grid = expand_grid(
  ncells = c(10000),
  n_patients = c(8, 12, 16, 20, 30, 40, 50),
  n_covariates = c(1:5),
  n_genes = 1000,
  seed = 1:10,
  continuous = c(TRUE, FALSE)
)
saveRDS(grid, "data/param_grid.rds")
