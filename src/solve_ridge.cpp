#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec ridge_solver_cpp(const arma::mat& X, const arma::vec& y, double lambda) {
    int n = X.n_cols;
    arma::mat I = arma::eye(n, n);
    return solve(X.t() * X + I * lambda, X.t() * y);
}
