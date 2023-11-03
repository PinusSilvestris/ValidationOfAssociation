library(ggplot2)
library(VoA)
SR5 <- get_sample_data('SR5', 5000)

result <- create_Q_plot(SR5$X, SR5$Y, MC=500)

