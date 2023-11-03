#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]


arma::mat arma_copula(arma::vec Rx,arma::vec Ry);
arma::mat calculate_copula_mc_grid(arma::vec X, arma::vec Y, int MC = 100);

