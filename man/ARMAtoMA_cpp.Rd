% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ARMAtoMA_cpp}
\alias{ARMAtoMA_cpp}
\title{Converting an ARMA Process to an Infinite MA Process}
\usage{
ARMAtoMA_cpp(ar, ma, lag_max)
}
\arguments{
\item{ar}{A \code{column vector} of length p}

\item{ma}{A \code{column vector} of length q}

\item{lag_max}{A \code{int} of the largest MA(Inf) coefficient required.}
}
\value{
A \code{column vector} containing coefficients
}
\description{
Takes an ARMA function and converts it to an infinite MA process.
}
\details{
This function is a port of the base stats package's ARMAtoMA. There is no significant speed difference between the two.
}
\examples{
# ARMA(2,1)
ARMAtoMA_cpp(c(1.0, -0.25), 1.0, 10)
# ARMA(0,1)
ARMAtoMA_cpp(numeric(0), 1.0, 10)
}
\author{
R Core Team and JJB
}
