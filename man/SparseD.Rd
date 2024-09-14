\name{SparseD}

\alias{SparseD}

\title{
General difference matrices for ordinary B-splines
}

\description{
Compute general difference matrices of order 1 to \eqn{d - 1} for ordinary B-splines of order \eqn{d}.
}

\usage{
SparseD(xt, d)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
}

\value{
A list of sparse matrices (of "dgCMatrix" class), giving general difference matrices of order 1 to \eqn{d - 1}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps.mgcv)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)

SparseD(xt, d = 4)
}
