
Rcpp::sourceCpp("multiway_sandwich.cpp")

# design_matrix: n x p
# y: length n
# beta: length p
# overdispersion: scalar phi
# size_factors: length n
# clusters: list of 1–3 integer vectors (patient, batch, [site])
compute_sandwich_multiway <- function(design_matrix, y, beta, overdispersion, size_factors, clusters) {
  # B = ( -H )^{-1} already returned by your compute_hessian
  B <- compute_hessian(beta, 1 / overdispersion, y, design_matrix, size_factors)  # matches your current call
  
  M <- compute_multiway_meat(
    design_matrix = design_matrix,
    y = y,
    beta = beta,
    overdispersion = overdispersion,
    size_factors = size_factors,
    clusters_list = clusters,            # e.g., list(patient, batch)
    finite_adj = TRUE
  )
  
  # keep your scaling convention
  n <- length(y)
  V <- (B %*% M %*% B) * n
  return(V)
}
