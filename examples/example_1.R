library(VoA)
library(ggplot2)

R <- seq(1,100)
S <- R^2

ret <- create_Q_plot(R, S)

ret$Q_plot
