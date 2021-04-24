library(ggplot2)
#' Plot Q function based on mone carlo estimation
#'
#' @description Greet a person and appropriately capitalize their name.
#'
#' @param X numerical vector with first random variable e.g. c(1.1, 2.2, 1.73).
#' @param Y numerical vector with second random variable e.g. c(3.1, 1.2, 1.93).
#' @param k_plot_grid number of grid points for plot e.g. 1000, defaults to 100.
#' @param MC number of replications e.g. 1000, defaults to 100.
#' @param print boolean value, if true a plot will be displayed, defaults to TRUE.
#'
#' @return A character string, capitalized to title case.
#' @export
#'
#' @examples
#' 'calculate_copula_mc_grid(c(1, 3, 2), c(3, 2, 1)) '
create_Q_plot <- function(X, Y, k_plot_grid = 100, MC = 100, print = TRUE){
  C_grid<- calculate_copula_mc_grid(X, Y, MC)

  Q_grid <- create_Q_grid(k_plot_grid)
  plot_points <- Q(Q_grid, C_grid)
  plot <- ggplot(plot_points, aes(x, y, z = z)) + geom_contour_filled()

  if(print) {
    print(plot)
  }

  return(plot)
}
