% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AJDN_Functions.R
\name{generate_upsilon_matrix}
\alias{generate_upsilon_matrix}
\title{Generate bootstrap matrix}
\usage{
generate_upsilon_matrix(x_matrix, s_prime)
}
\arguments{
\item{x_matrix}{\eqn{p} x \eqn{p} (\eqn{p} = # dimensions, \eqn{n} = # obsrevations) matrix of time series}

\item{s_prime}{block size parameter}
}
\value{
\eqn{p} x \eqn{n} matrix to be subsequently used in the block multiplier bootstrap
}
\description{
Generates the matrix used in the high dimensional block multiplier bootstrap
}
\examples{
x_matrix = matrix(stats::rnorm(1000),nrow=10)
upsilon_matrix = ajdn::generate_upsilon_matrix(x_matrix,1/100)
}
