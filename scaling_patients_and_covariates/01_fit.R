
rm(list = ls())
source("utils_sim.R")

param_grid = readRDS("data/param_grid.rds")
idx  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
row = param_grid[idx,]

ex <- simulate_splatter_simple(
  n_cells = row$ncells,
  n_patients = row$n_patients,
  n_covariates = row$n_covariates,  # -> 4 groups from two binary covariates
  n_genes = row$n_genes,
  seed = row$seed
)

input_matrix = as.matrix(ex$counts)
covariate_cols = colnames(ex$meta)[grepl("cov", colnames(ex$meta))]
form <- as.formula(paste("~", paste(covariate_cols, collapse = " + ")))
design_matrix <- model.matrix(form, data = ex$meta)

if (row$continuous) {
  for (i in 2:ncol(design_matrix)) {
    design_matrix[,i] = design_matrix[,i] + rnorm(nrow(design_matrix), 0, sd = .01)
  }
}

clusters = ex$meta$patient
contrast = rep(0, ncol(design_matrix))
contrast[1] = 1
input_matrix = input_matrix[rowSums(input_matrix) > 500, ]

# Fit devil
s = Sys.time()
devil.fit = devil::fit_devil(input_matrix, design_matrix, size_factors = NULL, init_overdispersion = 100,
                             parallel.cores = 1, overdispersion = "old", offset = 1e-6, max_iter = 500)
devil.test = devil::test_de(devil.fit, contrast = contrast, clusters = clusters, parallel.cores = 1)
e = Sys.time()
devil_time = e - s

# Fit glmGamPoi
s = Sys.time()
glm.fit = glmGamPoi::glm_gp(input_matrix, design_matrix, size_factors = FALSE)
glm.test = glmGamPoi::test_de(glm.fit, contrast = contrast)
e = Sys.time()
glmGamPoi_time = e - s

# Fit NEBULA
s = Sys.time()
nebula.fit = nebula::nebula(input_matrix, pred = design_matrix, id = clusters, ncore = 1, cpc = 0, mincp = 0)
e = Sys.time()
nebula_time = e - s

r = dplyr::bind_cols(row, 
                     dplyr::tibble(
                       NEBULA = as.numeric(nebula_time, units = "secs"), 
                       devil = as.numeric(devil_time, units = "secs"), 
                       glmGamPoi = as.numeric(glmGamPoi_time, units = "secs"))
                     )

saveRDS(r, paste0("results/sim_", idx, ".rds"))
