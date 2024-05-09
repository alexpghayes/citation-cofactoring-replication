#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec lower_triangle_rowsums_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const int num_threads) {

  int n = U.n_rows;
  arma::mat DVt = arma::diagmat(d) * V.t();

  arma::vec rs = arma::zeros<arma::vec>(n);

  omp_set_num_threads(num_threads);

  #pragma omp parallel for
  for (int i = 1; i < n; i++) {
    rs(i) = arma::accu(U.row(i) * DVt.cols(0, i - 1));
  }

  return rs;
}

// [[Rcpp::export]]
arma::vec lower_triangle_abs_rowsums_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const int num_threads) {

  int n = U.n_rows;
  arma::mat DVt = arma::diagmat(d) * V.t();

  arma::vec rs = arma::zeros<arma::vec>(n);

  omp_set_num_threads(num_threads);

  #pragma omp parallel for
  for (int i = 1; i < n; i++) {
    rs(i) = arma::accu(arma::abs(U.row(i) * DVt.cols(0, i - 1)));
  }

  return rs;
}


// [[Rcpp::export]]
arma::vec lower_triangle_abs_colsums_impl(
    const arma::mat& U,
    const arma::rowvec& d,
    const arma::mat& V,
    const int num_threads) {

  int n = V.n_rows;
  arma::mat DVt = arma::diagmat(d) * V.t();

  arma::vec cs = arma::zeros<arma::vec>(n);

  omp_set_num_threads(num_threads);

  #pragma omp parallel for
  for (int i = 0; i < (n - 1); i++) {
    cs(i) = arma::accu(arma::abs(U.rows(i + 1, n - 1) * DVt.col(i)));
  }

  return cs;
}
