% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rfilter}
\alias{rfilter}
\title{Time Series Recursive Filters}
\usage{
rfilter(x, filter, init)
}
\arguments{
\item{x}{A \code{column vector} of length T}

\item{filter}{A \code{column vector} of length f}

\item{init}{A \code{column vector} of length f that contains the initial values of the time series in reverse.}
}
\value{
x A \code{column vector} with its contents reversed.
}
\description{
Applies a recursive filter to a univariate time series.
}
\details{
Note: The length of 'init' must be equal to the length of 'filter'.
This is a port of the rfilter function harnessed by the filter function in stats. 
It is about 6-7 times faster than R's base function. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
}
\examples{
x = 1:100
# 
rfilter(x, rep(1, 3), rep(1, 3))
# Using R's function
filter(x, rep(1, 3), method="recursive", init=rep(1, 3))
}
\author{
R Core Team and JJB
}
