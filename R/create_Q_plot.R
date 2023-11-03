#' Plot Q function based on mone carlo estimation
#'
#' @description Plot Q function based on mone carlo estimation
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
  time_taken <- system.time({
    C_grid<- calculate_copula_mc_grid(X, Y, MC)

    Q_grid <- create_Q_grid(k_plot_grid)
    plot_points <- Q(Q_grid, C_grid)
    plot <- ggplot(plot_points, aes(x, y, z = z)) + geom_contour_filled()

    if(print) {
      print(plot)
    }
  })[3]

  cat("Time taken for calculations:", time_taken, "seconds\n")

  return(list(Q_plot = plot, C_grid = C_grid, Q_grid = Q_grid, plot_points = plot_points))
}
