#' Calculate copula based on single sample
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
#' 'calculate_copula_grid(c(1, 3, 2), c(3, 2, 1)) '
calculate_copula_grid <- function(X, Y) {
  R <- rank(X, ties.method = "random")
  S <- rank(Y, ties.method = "random")
  return(arma_copula(R, S))
}



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
#' 'calculate_copula_mc_grid(c(1, 3, 2), c(3, 2, 1))'
calculate_copula_mc_grid <- function(X, Y, MC=100){

  g <- NULL
  K <- length(X)
  pb <- txtProgressBar(min=0, max=MC, style=3)

  for(i in 1:MC){

    s <- sample(seq(1:K), K, TRUE)
    X_s <- X[s]
    Y_s <- Y[s]

    g_tmp <- calculate_copula_grid(X_s, Y_s)

    if(is.null(g)) {
      g <- g_tmp
    } else {
      g <- g+g_tmp
    }
    setTxtProgressBar(pb, i)
  }

  g <- g/MC
  return(g)

}
