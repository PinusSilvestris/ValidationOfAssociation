install.packages("devtools")
devtools::install_github("PinusSilvestris/ValidationOfAssociation")

library(VoA)
library(ggplot2)

R <- seq(1,10)
S <- R^2

ret <- create_Q_plot(R, S, k_plot_grid = 10, MC = 1)

ret$Q_plot

SR5 <- get_sample_data('SR5', 1000)

result <- create_Q_plot(SR5$X, SR5$Y, MC=100)
