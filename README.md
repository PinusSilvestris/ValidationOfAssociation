# Install

```R
install.packages("devtools")
devtools::install_github("PinusSilvestris/ValidationOfAssociation")
```

# Usage

```R
library(ggplot2)
library(VoA)
SR5 <- get_sample_data('SR5', 1000)

result <- create_Q_plot(SR5$X, SR5$Y, MC=100)
```
