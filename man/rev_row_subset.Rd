% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rev_row_subset}
\alias{rev_row_subset}
\title{Reverse Subset Row}
\usage{
rev_row_subset(x, start, end)
}
\arguments{
\item{x}{A \code{matrix} of dimensions M x N}

\item{start}{A \code{unsigned int} that indicates the starting row.}

\item{end}{A \code{unsigned int} that indicates the ending row.}
}
\value{
x A \code{matrix} with matrix rows displayed in reversed order
}
\description{
Subsets the row by going from high indices to low (the reverse of the supported practice)
}
\details{
Consider a vector x=[[1,2],[3,4]], the function would output x=[[3,4],[1,2]].
Start and end must be valid C++ matrix locations. (e.g. matrix rows start at 0 and not 1)
}
\examples{
x = matrix(c(1,2,3,4), nrow=2,byrow=TRUE)
rev_row_subset(x, 1, 0)
}
\author{
JJB
}
