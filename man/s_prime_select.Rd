% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim_Functions.R
\name{s_prime_select}
\alias{s_prime_select}
\title{Choose \eqn{s'}}
\usage{
s_prime_select(dgp, n, p, max_s_prime = NULL)
}
\arguments{
\item{dgp}{data generating process, one of: IID, IIDT, GS, GST, PS, PST, LS, LST, PLS, PLST}

\item{n}{length of time series}

\item{p}{dimension of time series}

\item{max_s_prime}{OPTIONAL argument, the maximium \eqn{ns'} to consider (default is \eqn{n^{1/3}})}
}
\value{
value of \eqn{s'}
}
\description{
This function computes \eqn{s'} based on the desired LRV ratio from Appendix B.1.2, for data generating processes used in the simulations.
}
\examples{
ajdn::s_prime_select('GS',1000,100)
}
