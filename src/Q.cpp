// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "rcpp_copula.h"

// [[Rcpp::export]]
arma::mat create_Q_grid(int n = 10) {
  arma::vec u = arma::linspace<arma::vec>(0.5/(n+1), (0.5 + n)/(n+1), n+1);
  arma::mat grid((n+1)*(n+1), 2);

  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= n; ++j) {
      grid(i * (n+1) + j, 0) = u(i);
      grid(i * (n+1) + j, 1) = u(j);
    }
  }

  return grid;
}

// [[Rcpp::export]]
Rcpp::DataFrame Q(const arma::mat& Q_grid, const arma::mat& C_grid) {
  int n = C_grid.n_rows - 1;
  Rcpp::NumericVector x = Rcpp::wrap(Q_grid.col(0));
  Rcpp::NumericVector y = Rcpp::wrap(Q_grid.col(1));
  Rcpp::NumericVector z(Q_grid.n_rows);

  for (arma::uword i = 0; i < Q_grid.n_rows; ++i) {
    int row_idx = std::ceil(Q_grid(i, 0) * (n + 1)) - 1;
    int col_idx = std::ceil(Q_grid(i, 1) * (n + 1)) - 1;
    if (row_idx > n) row_idx = n;
    if (col_idx > n) col_idx = n;
    z[i] = C_grid(row_idx, col_idx);
  }

  return Rcpp::DataFrame::create(Rcpp::Named("x") = x,
                                 Rcpp::Named("y") = y,
                                 Rcpp::Named("z") = z);
}
