library(ggplot2)
#' Calculate copula based on monte carlo samples
#'
#' @description Greet a person and appropriately capitalize their name.
#'
#' @param X Your name (character string; e.g. "john doe").
#' @param Y Your name (character string; e.g. "john doe").
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
