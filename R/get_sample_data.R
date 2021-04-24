#' Calculate copula based on monte carlo samples
#'
#' @description Greet a person and appropriately capitalize their name.
#'
#' @param type One of the examples from the paper i.e. SR1-SR5
#' @param n number of samples to generate
#'
#' @return list(X=..., Y=...) List of two vectors X and Y generatted accordingly to the type
#' @export
#'
#' @examples
#' 'calculate_copula_mc_grid(c(1, 3, 2), c(3, 2, 1)) '
get_sample_data <- function(type, n) {
  data <- NULL
  switch(type,
         SR1 = {  X<-runif(n)
                  Y<-2+X+rnorm(n)},
         SR2 = {  X<-runif(n)
                  Y<-X**0.25+rnorm(n, 0, 0.25)},
         SR3 = {  X<-runif(n)
                  Y<-(X>=0.5)+rnorm(n, 0,2)},
         SR4 = {  X<-rnorm(n)
                  Y<-log(1+abs(X))+rnorm(n, 0,1)},
         SR5 = {  X<-runif(n)
                  Y<-4*((2*X-1)**2-0.5)**2+rnorm(n, 0, 0.5)}
  )
  data <- list(X=X, Y=Y)
  return(data)
}
