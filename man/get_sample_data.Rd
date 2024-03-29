% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sample_data.R
\name{get_sample_data}
\alias{get_sample_data}
\title{Calculate copula based on monte carlo samples}
\usage{
get_sample_data(type, n)
}
\arguments{
\item{type}{One of the examples from the paper i.e. SR1-SR5}

\item{n}{number of samples to generate}
}
\value{
list(X=..., Y=...) List of two vectors X and Y generatted accordingly to the type
}
\description{
Greet a person and appropriately capitalize their name.
}
\examples{
'calculate_copula_mc_grid(c(1, 3, 2), c(3, 2, 1)) '
}
