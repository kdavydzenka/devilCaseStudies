#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd compute_scores(const Eigen::MatrixXd& design_matrix,
                               const Eigen::VectorXd& y,
                               const Eigen::VectorXd& beta,
                               const double overdispersion,
                               const Eigen::VectorXd& size_factors) {
  double alpha = 1.0 / overdispersion;
  
  // Vectorized computation
  VectorXd eta = design_matrix * beta;
  VectorXd mu = size_factors.array() * eta.array().exp();
  VectorXd residuals = (y.array() - mu.array()) / mu.array();
  VectorXd weights = mu.array() / (1.0 + mu.array() / alpha);
  VectorXd wr = residuals.array() * weights.array();
  
  return design_matrix.array().colwise() * wr.array();
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_hessian_expected(const Eigen::VectorXd& beta,
                                         const double overdispersion,
                                         const Eigen::VectorXd& y,
                                         const Eigen::MatrixXd& design_matrix,
                                         const Eigen::VectorXd& size_factors) {
  const double alpha = 1.0 / overdispersion;
  const int n = design_matrix.rows();
  // mu = s * exp(X beta)
  Eigen::VectorXd eta = design_matrix * beta;
  Eigen::VectorXd mu  = size_factors.array() * eta.array().exp();
  
  // w_i = mu_i / (1 + mu_i/alpha)
  Eigen::VectorXd w = mu.array() / (1.0 + mu.array() / alpha);
  
  // A = X^T W X (W diagonal)
  Eigen::MatrixXd XTW = design_matrix.transpose() * w.asDiagonal();
  Eigen::MatrixXd A   = XTW * design_matrix;   // this is the negative expected Hessian
  return A;                                    // return A, not A^{-1}
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_clustered_meat(const Eigen::MatrixXd& design_matrix,
                                       const Eigen::VectorXd& y,
                                       const Eigen::VectorXd& beta,
                                       const double overdispersion,
                                       const Eigen::VectorXd& size_factors,
                                       const Eigen::VectorXi& clusters) {
  const int n = design_matrix.rows();
  const int p = design_matrix.cols();
  
  // scores: n x p (row i is s_i^T)
  Eigen::MatrixXd S = compute_scores(design_matrix, y, beta, overdispersion, size_factors);
  
  // get unique cluster ids
  std::vector<int> ids;
  ids.reserve(n);
  std::unordered_set<int> seen;
  for (int i = 0; i < n; ++i) {
    if (!seen.count(clusters(i))) { seen.insert(clusters(i)); ids.push_back(clusters(i)); }
  }
  
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(p, p);
  for (int id : ids) {
    Eigen::VectorXd sum_g = Eigen::VectorXd::Zero(p);
    for (int i = 0; i < n; ++i) if (clusters(i) == id) sum_g += S.row(i).transpose();
    B.noalias() += sum_g * sum_g.transpose();
  }
  return B; // no 1/n here
}
