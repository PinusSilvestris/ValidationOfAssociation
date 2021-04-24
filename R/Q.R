create_Q_grid <- function(n=10){
  u <- (0.5+seq(0, n))/(n+1)
  v <- (0.5+seq(0, n))/(n+1)
  return(expand.grid(x=u, y=v))
}

Q <- function(Q_grid, C_grid) {

  n <- nrow(C_grid)-1

  get_val <- function(x, y){
    i <- ceiling(x*(n+1))
    j <- ceiling(y*(n+1))

    return(C_grid[i,j])
  }
  Q_grid$z <- mapply(get_val, Q_grid$x, Q_grid$y)
  return(Q_grid)
}
