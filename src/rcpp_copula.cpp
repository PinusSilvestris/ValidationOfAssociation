// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <atomic>
#include <chrono>
#include <thread>

// [[Rcpp::depends(RcppArmadillo)]]

// Simple function to display the progress
void print_progress(std::atomic<int>& counter, int total) {
  while (counter < total) {
    int progress = counter.load();
    Rcpp::Rcout << "\rProgress: " << progress << "/" << total << " ("
                << (100 * progress / total) << "%)";
    Rcpp::Rcout.flush();
    std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Update every 100ms
  }
  Rcpp::Rcout << "\rProgress: " << total << "/" << total << " (100%)\n";
}

// Swap ranks
arma::vec swap_ranks(const arma::vec& R){
  return (R.size() + 1.0) - R;
}

// The Wn function
double Wn(int n, double C, double i, double j)
{
  double u = (i + 0.5) / (n + 1);
  double v = (j + 0.5) / (n + 1);
  double w = std::sqrt(n) * (C - u * v) /
    std::sqrt(u * v * (1.0 - u) * (1.0 - v));
  return w;
}

// Fill vector T from ranks
arma::vec fill_T(const arma::vec &Rx, const arma::vec &Ry){
  int n = Rx.size();
  arma::vec T = arma::zeros(n);

  for (int i = 0; i < n; i++){
    T((int)Rx[i] - 1) = Ry[i];
  }
  return T;
}

// Calculate part of the empirical copula
void calculate_copula_part(const arma::vec &T, arma::mat &C) {
  int n = T.size();
  double n_inv = 1.0 / n;
  C.zeros(); // Reset C to zero before reuse

  // Fill upper-left quadrant of the C matrix
  for (int i = 1; i <= n / 2 + 1; ++i) {
    for (int j = 0; j <= n / 2 + 1; ++j) {
      C(i, j) = C(i - 1, j) + (j >= T[i - 1]) * n_inv;
    }
  }
}

// Fill KS matrix with Wn values
void fill_matrix(int n, arma::mat &KS, const arma::mat &C,
                 bool x_prim, bool y_prim)
{
  int sign = std::pow(-1, x_prim + y_prim);

  for (int i = 0; i <= (n + 1) / 2; i++){
    for (int j = 0; j <= (n + 1) / 2; j++){
      int x = x_prim ? n - i : i;
      int y = y_prim ? n - j : j;
      KS(x, y) = sign * Wn(n, C(i, j), i, j);
    }
  }
}

// Core function that constructs the copula KS matrix
// [[Rcpp::export]]
arma::mat arma_copula(const arma::vec &Rx, const arma::vec &Ry) {
  int n = Rx.size();

  arma::mat Ctab      = arma::zeros(n+1, n+1);
  arma::mat Ctabs22   = arma::zeros(n+1, n+1);
  arma::mat Ctabs12   = arma::zeros(n+1, n+1);
  arma::mat Ctabs21   = arma::zeros(n+1, n+1);

  arma::mat KS        = arma::zeros(n+1, n+1);

  // Swapped versions of Rx, Ry
  arma::vec Rsx = swap_ranks(Rx);
  arma::vec Rsy = swap_ranks(Ry);

  // Fill T
  arma::vec T    = fill_T(Rx,  Ry);
  arma::vec Ts22 = fill_T(Rsx, Rsy);
  arma::vec Ts12 = fill_T(Rx,  Rsy);
  arma::vec Ts21 = fill_T(Rsx, Ry);

  // Calculation of empirical copulas
  calculate_copula_part(T,    Ctab);
  calculate_copula_part(Ts22, Ctabs22);
  calculate_copula_part(Ts12, Ctabs12);
  calculate_copula_part(Ts21, Ctabs21);

  // Fill the KS matrix with all 4 quadrants
  fill_matrix(n, KS, Ctab,    false, false);
  fill_matrix(n, KS, Ctabs22, true,  true );
  fill_matrix(n, KS, Ctabs12, false, true );
  fill_matrix(n, KS, Ctabs21, true,  false);

  return KS;
}

// Rank function
arma::uvec rank(const arma::vec& data) {
  // 1-based ranking
  return (arma::sort_index(arma::sort_index(data, "ascend")) + 1);
}

// [[Rcpp::export]]
arma::mat calculate_copula_grid(const arma::vec &X, const arma::vec &Y) {
  // Simple dimension check
  if (X.size() != Y.size()) {
    Rcpp::stop("Size of X and Y do not match");
  }

  arma::vec R = arma::conv_to<arma::vec>::from(rank(X));
  arma::vec S = arma::conv_to<arma::vec>::from(rank(Y));
  return arma_copula(R, S);
}

// [[Rcpp::export]]
arma::mat calculate_copula_mc_grid(const arma::vec &X, const arma::vec &Y, int MC) {
  // Ensure same lengths
  if (X.size() != Y.size()) {
    Rcpp::stop("Size of X and Y do not match");
  }

  int K = X.n_elem;
  arma::mat g = arma::zeros<arma::mat>(K + 1, K + 1);

  // Create atomic counter for progress
  std::atomic<int> progress_counter(0);

  // Start a thread to print progress
  std::thread progress_thread(print_progress, std::ref(progress_counter), MC);

  // Variables to capture exceptions from parallel region
  bool error_occurred = false;
  std::string error_msg;

  // OpenMP parallel section
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < MC; ++i) {
    // If an error is already recorded, skip further work
    if (error_occurred) {
      continue;
    }
    try {
      // Draw random indices
      arma::uvec sample_indices = arma::randi<arma::uvec>(K, arma::distr_param(0, K-1));
      arma::vec X_s = X.elem(sample_indices);
      arma::vec Y_s = Y.elem(sample_indices);

      // Compute one realization
      arma::mat mat = calculate_copula_grid(X_s, Y_s);

      // Accumulate results
#ifdef _OPENMP
#pragma omp critical
#endif
{
  g += mat / MC;
}

// Update progress
progress_counter++;
    }
    catch (std::exception &ex) {
      // Record first error in a thread-safe manner
#ifdef _OPENMP
#pragma omp critical
#endif
{
  if (!error_occurred) {
    error_occurred = true;
    error_msg = ex.what();
    // Force progress to end so thread can join quickly
    progress_counter = MC;
  }
}
    }
    catch (...) {
      // Catch-all for non-std::exception
#ifdef _OPENMP
#pragma omp critical
#endif
{
  if (!error_occurred) {
    error_occurred = true;
    error_msg = "Unknown error in calculate_copula_mc_grid.";
    progress_counter = MC;
  }
}
    }
  } // end parallel for

  // Join the progress thread
  progress_thread.join();

  // If any thread hit an error, stop() now
  if (error_occurred) {
    Rcpp::stop(error_msg);
  }

  return g;
}
