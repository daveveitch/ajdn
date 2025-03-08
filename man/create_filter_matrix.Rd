% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AJDN_Functions.R
\name{create_filter_matrix}
\alias{create_filter_matrix}
\title{Create filter matrix}
\usage{
create_filter_matrix(n, unique_scale_counts)
}
\arguments{
\item{n}{\eqn{n} = # obsrevations of time series}

\item{unique_scale_counts}{vector of \eqn{s_jn} for all unique \eqn{s_j} that are being used across all dimensions}
}
\value{
List of two items
1st entry - sparse matrix that represents
2nd entry - filter matrix dictionary which specifies what scale_count (i.e. \eqn{sn}) and \eqn{t_j n} each
row of the filter matrix corresponds to.
}
\description{
This function creates a filter matrix, as well as a filter matrix dictionary
which is then used to calculate both test statistics and bootstrap statistics for AJDN. The entires of this matrix
are equal to the \eqn{W((j/n-t)/s)} components of \eqn{H(t,s,r)} for various values of \eqn{j,t,s}.
}
\examples{
filter_matrix = ajdn::create_filter_matrix(1000,c(50,100))
plot(filter_matrix[[1]][100,])
print(filter_matrix[[2]][100,])
}
