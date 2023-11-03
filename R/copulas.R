#' if (!require(parallel)) install.packages("parallel")
#' library(parallel)
#'
#' if (!require(pbmcapply)) {
#'   install.packages("pbmcapply")
#' }
#' library(pbmcapply)
#'
#'
#' Calculates copula based on a given sample
#'
#' #' @description Calculates copula based on a given sample.
#' #'
#' #' @param X numerical vector with first random variable e.g. c(1.1, 2.2, 1.73).
#' #' @param Y numerical vector with second random variable e.g. c(3.1, 1.2, 1.93).
#' #'
#' #' @return 2d matrix of Copula function calculated on a grid
#' #' @export
#' #'
#' #' @examples
#' #' 'calculate_copula_grid(c(1, 3, 2), c(3, 2, 1)) '
#' calculate_copula_grid <- function(X, Y) {
#'   R <- rank(X, ties.method = "random")
#'   S <- rank(Y, ties.method = "random")
#'   return(arma_copula(R, S))
#' }
#'
#' #' Calculate copula based on monte carlo samples
#' #'
#' #' @description Calculates copula based on a given sample using bootstrap.
#' #'
#' #' @param X numerical vector with first random variable e.g. c(1.1, 2.2, 1.73).
#' #' @param Y numerical vector with second random variable e.g. c(3.1, 1.2, 1.93).
#' #' @param MC number of replications e.g. 1000, defaults to 100.
#' #'
#' #' @return 2d matrix of Copula function calculated on a grid
#' #' @export
#' #'
#' #' @examples
#' #' 'calculate_copula_mc_grid(c(1, 3, 2), c(3, 2, 1))'
#' # Ensure the 'pbapply' package is installed and loaded for the progress bar
#' calculate_copula_mc_grid <- function(X, Y, MC = 100, no_cores = 4) {
#'   K <- length(X)
#'
#'   # Pre-allocate the matrix for the result
#'   if (!exists("calculate_copula_grid")) {
#'     stop("The function 'calculate_copula_grid' does not exist in the current environment.")
#'   }
#'   result_size <- dim(calculate_copula_grid(X, Y))
#'   g <- matrix(0, nrow = result_size[1], ncol = result_size[2])
#'
#'   # Efficient sampling
#'   # sample_indices <- matrix(sample.int(K, K * MC, replace = TRUE), ncol = K, byrow = TRUE)
#'
#'   # Set up the parallel backend to use multiple cores
#'   if(no_cores == -1) { no_cores <- parallel::detectCores() - 1 }
#'
#'   # Run the copula calculations in parallel with a progress bar
#'   time_taken <- system.time({
#'     # Use mclapply for better memory efficiency
#'     pbmclapply(1:MC, function(i) {
#'       # Sample indices for the current iteration
#'       sample_indices <- sample.int(length(X), length(X), replace = TRUE)
#'       X_s <- X[sample_indices]
#'       Y_s <- Y[sample_indices]
#'       mat <- calculate_copula_grid(X_s, Y_s)
#'
#'       # Update the result matrix within the loop to avoid storing all matrices
#'       parallel::mccollect(g <<- g + mat / MC, wait = FALSE)
#'     }, mc.cores = no_cores)
#'   })[3]  # Extract the elapsed time
#'
#'   # Print the time taken
#'   cat("Time taken for calculations:", time_taken, "seconds\n")
#'
#'   return(g)
#' }
