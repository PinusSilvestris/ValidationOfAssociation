// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "rcpp_copula.h"
#include <RcppArmadilloExtensions/sample.h>
# include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec seq_int(long int a, long int b){
  long int d = std::abs(b-a)+1;
  return arma::linspace<arma::uvec>(a, b, d);
}

// [[Rcpp::export()]]
arma::uvec rnd_idx(int n = 10) {
  return  Rcpp::RcppArmadillo::sample(seq_int(0, n-1), n, true);
}

// [[Rcpp::export()]]
std::vector<arma::uvec> bootstrap_samples(int n, int MC=100){
  std::vector<arma::uvec> bootstrap = std::vector<arma::uvec>(MC);

  for (int t = 0; t < MC; ++t) {
    arma::uvec resample = rnd_idx(n);
    bootstrap[t] = resample;
  }
  return bootstrap;
}




// [[Rcpp::export()]]
arma::mat arma_copula_mc(arma::vec Rx, arma::vec Ry, int MC=100, int t=1) {
  int n = Rx.size();
  std::vector<arma::uvec> bootstrap_idx = bootstrap_samples(n, MC);
  //std::vector<arma::mat> copulas = std::vector<arma::mat>(MC);
  arma::mat sum = arma::zeros(n+1, n+1);

  omp_set_num_threads(t) ;
#pragma omp parallel for reduction(+:sum)
  for (int b = 0 ; b < MC ; b++) {
    arma::uvec idx = bootstrap_idx[b];

    arma::vec Rx_sample = Rx.elem(idx);
    arma::vec Ry_sample = Ry.elem(idx);

    arma::mat C_sample = arma_copula(Rx_sample, Ry_sample);
    sum += C_sample;
  }

  return sum/MC;

}

// [[Rcpp::export()]]
void omp2 (int t = 1) {
  omp_set_num_threads(t) ;
# pragma omp parallel for
  for (int i = 0 ; i < 10 ; i++) {
    Rcpp::Rcout << " " << i << " " ;
  }
  Rcpp::Rcout << std::endl ;
}

/*** R

*/
