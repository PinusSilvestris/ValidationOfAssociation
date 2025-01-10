// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <omp.h>
#include "rcpp_copula.h"
#include <iostream>
#include <atomic>
#include <chrono>
#include <thread>

// Simple function to display the progress
void print_progress(std::atomic<int>& counter, int total) {
  while (counter < total) {
    int progress = counter.load();
    std::cout << "\rProgress: " << progress << "/" << total << " (" << (100 * progress / total) << "%)";
    std::cout.flush();
    std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Update every 100ms
  }
  std::cout << "\rProgress: " << total << "/" << total << " (100%)\n";
}
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec swap_ranks(arma::vec R){
  return (R.size() + 1.0) - R;
}

double Wn(int n, double C, double i, double j)
{
  double u = (i+0.5)/(n+1);
  double v = (j+0.5)/(n+1);
  double w = pow(n, 0.5) * (C - u * v) / pow(u * v * (1.0 - u) * (1.0 - v), 0.5);
  return w;
}

arma::vec fill_T(arma::vec &Rx, arma::vec &Ry){
  int n = Rx.size();
  arma::vec T = arma::zeros(n);

  for (int i = 0; i < n; i++)
  {
    T((int)Rx[i] - 1) = Ry[i];
  }
  return T;
}

void calculate_copula_part(arma::vec &T, arma::mat &C){
  int n = T.size();
  double n_inv = 1.0 / n;
  C.zeros(); // Reset C to zero before reuse

  for (int i = 1; i <= n / 2 + 1; ++i) {
    for (int j = 0; j <= n / 2 + 1; ++j) {
      C(i, j) = C(i - 1, j) + (j >= T[i - 1]) * n_inv;
    }
  }}

void fill_matrix(int n, arma::mat &KS, arma::mat &C, bool x_prim, bool y_prim){

  int sign = pow(-1, x_prim+y_prim);

  //Rcpp::Rcout<<"x_prim: "<<x_prim<<" y_prim"<<y_prim<<std::endl;
  for (int i = 0; i <= (double)(n + 1) / 2.0; i++)
  {
    for (int j = 0; j <= (double)(n + 1) / 2.0; j++)
    {
      int x = x_prim ? n-i : i;
      int y = y_prim ? n-j : j;
      KS(x, y) = sign*Wn(n, C(i, j), i, j);
      //Rcpp::Rcout<<"i: "<<i<<" j:"<<" x:"<<x<<" y:"<<y<<" sign:"<<sign<<" Wn:"<<Wn(n, C(i, j), i, j)<<std::endl;
    }
  }}

// [[Rcpp::export]]
arma::mat arma_copula(arma::vec Rx, arma::vec Ry) {
  int n = Rx.size();
  // TODO: decrease memory footprint
  // int k = (int)(1.0*n/2.0);

  arma::mat Ctab = arma::zeros(n+1, n+1);
  arma::mat Ctabs22 = arma::zeros(n+1, n+1);
  arma::mat Ctabs12 = arma::zeros(n+1, n+1);
  arma::mat Ctabs21 = arma::zeros(n+1, n+1);
  arma::vec Rsx = arma::zeros(n);
  arma::vec Rsy = arma::zeros(n);
  arma::vec T;
  arma::vec Ts22;
  arma::vec Ts12;
  arma::vec Ts21;

  arma::mat KS = arma::zeros(n+1, n+1);

  Rsx = swap_ranks(Rx);
  Rsy = swap_ranks(Ry);

  T = fill_T(Rx, Ry);
  Ts22 = fill_T(Rsx, Rsy);
  Ts12 = fill_T(Rx, Rsy);
  Ts21 = fill_T(Rsx, Ry);

  // calculation of empirical copula on the grid, for ranks and transformed ranks
  calculate_copula_part(T, Ctab);
  calculate_copula_part(Ts22, Ctabs22);
  calculate_copula_part(Ts12, Ctabs12);
  calculate_copula_part(Ts21, Ctabs21);

  fill_matrix(n, KS, Ctab, false, false);
  fill_matrix(n, KS, Ctabs22, true, true);
  fill_matrix(n, KS, Ctabs12, false, true);
  fill_matrix(n, KS, Ctabs21, true, false);

  return KS;
}

arma::uvec rank(const arma::vec& da_ta) {
  return (arma::sort_index(arma::sort_index(da_ta, "ascend")) + 1);
}

// [[Rcpp::export]]
arma::mat calculate_copula_grid(const arma::vec X, const arma::vec Y) {
  arma::vec R = arma::conv_to<arma::vec>::from(rank(X));
  arma::vec S = arma::conv_to<arma::vec>::from(rank(Y));
  return arma_copula(R, S);
}

// [[Rcpp::export]]
arma::mat calculate_copula_mc_grid(const arma::vec& X, const arma::vec& Y, int MC) {
  int K = X.n_elem;
  arma::mat g = arma::zeros<arma::mat>(X.n_elem+1, Y.n_elem+1); // Adjust size as necessary
  std::atomic<int> progress_counter(0);
  std::thread progress_thread(print_progress, std::ref(progress_counter), MC);

  if(X.size() != Y.size()){
    std::cout<<"Size of X and Y do not match";
    throw std::invalid_argument( "Size of X and Y do not match" );
  }

#pragma omp parallel for
  for (int i = 0; i < MC; ++i) {
    arma::uvec sample_indices = arma::randi<arma::uvec>(K, arma::distr_param(0, K-1));
    arma::vec X_s = X.elem(sample_indices);
    arma::vec Y_s = Y.elem(sample_indices);
    arma::mat mat = calculate_copula_grid(X_s, Y_s);
    progress_counter++;
    //#pragma omp critical
    g += mat / MC;
  }
  progress_thread.join(); // Wait for the progress display to finish

  return g;
}
