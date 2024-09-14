\name{SandBar}

\alias{SandBar}

\title{
Build derivative penalty matrix or its sparse "root"
}

\description{
Compute order-\eqn{m} derivative penalty matrix \eqn{\bm{S}} or its sparse "root" \eqn{\bm{K}} such that \eqn{\bm{K'K} = \bm{S}}.
}

\usage{
SandBar(xt, d, m)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the derivative penalty (\eqn{0 \le m \le d - 1} required).}
}

\details{
Build \eqn{\bm{\bar{S}}} in the sandwich formula \eqn{\bm{S} = \bm{D'\bar{S}D}}.
}

\value{
A sparse matrice of "dsCMatrix" (when \eqn{m < d - 1}) or "ddiMatrix" (when \eqn{m = d - 1}) class.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps.mgcv)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)

SandBar(xt, d = 4, m = 2)
}
