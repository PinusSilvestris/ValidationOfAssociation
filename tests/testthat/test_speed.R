context("Speed tests")
library(microbenchmark)

D <- get_sample_data("SR5", 500)
R <- rank(D$X, ties.method = "random")
S <- rank(D$Y, ties.method = "random")


microbenchmark("Rcpp_mc_100" = {
  d <- calculate_copula_mc_grid(D$X, D$Y, 100)
},
"Rcpp_mc_100_t_1" = {
  d <- arma_copula_mc(R, S, 100, 8)
},
times = 5)


D <- get_sample_data("SR5", 1000)
R <- rank(D$X, ties.method = "random")
S <- rank(D$Y, ties.method = "random")

microbenchmark("Rcpp_mc_10000" = {
  d <- calculate_copula_mc_grid(D$X, D$Y, 10000)
},
"Rcpp_mc_10000_t_8" = {
  d <- arma_copula_mc(R, S, 10000, 8)
},
times = 5)

microbenchmark(
"Rcpp_mc_10000_t_8" = {
  d <- arma_copula_mc(R, S, 10000, 16)
},
times = 1)
