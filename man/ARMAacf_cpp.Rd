% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ARMAacf_cpp}
\alias{ARMAacf_cpp}
\title{Compute Theoretical ACF for an ARMA Process}
\usage{
ARMAacf_cpp(ar,ma,lag_max)
}
\arguments{
\item{ar}{A \code{vector} of length p containing AR coefficients}

\item{ma}{A \code{vector} of length q containing MA coefficients}

\item{lag_max}{A \code{unsigned integer} indicating the maximum lag necessary}
}
\value{
x A \code{matrix} listing values from 1...nx in one column and 1...1, 2...2,....,n...n, in the other
}
\description{
Compute the theoretical autocorrelation function for an ARMA process.
}
\details{
This is an implementaiton of the ARMAacf function in R. It is approximately 40x times faster. The benchmark was done on iMac Late 2013 using vecLib as the BLAS.
}
\examples{
# ARMA(2,1)
ARMAacf_cpp(c(1.0, -0.25), 1.0, lag_max = 10)
# ARMA(0,1)
ARMAacf_cpp(numeric(0), .35, lag_max = 10)
}
\author{
R Core Team and JJB
}
