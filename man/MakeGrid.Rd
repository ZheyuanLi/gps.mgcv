\name{MakeGrid}

\alias{MakeGrid}

\title{
Make a grid on the domain
}

\description{
Place equidistant grid points on each knot span to make a grid on the domain, suitable for evaluating B-splines.
}

\usage{
MakeGrid(xd, n, rm.dup = FALSE)
}

\arguments{
\item{xd}{domain knot sequence.}
\item{n}{number of equidistant grid points on each knot span.}
\item{rm.dup}{if FALSE, interior knots will appear twice on the grid; if TRUE, they will appear only once.}
}

\details{
Denote the domain knot sequence by \eqn{s_0, s_1, s_2, \ldots, s_k, s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain endpoints. A knot span is the interval between two successive knots, i.e., \eqn{[s_j, s_{j + 1}]}.

To make a suitable grid on \eqn{[a, b]} for evaluating B-splines, we can place \eqn{n} equidistant grid points on each knot span, ending up with \eqn{n(k + 1)} grid points in total. Interior knots will show up twice in this grid. To keep only one instance, set \code{rm.dup = TRUE}.
}

\value{
A vector of grid points.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps.mgcv)

## 4 domain knots: two interior knots 0.5 and 1.5 in domain [0, 3]
xd <- c(0, 0.5, 1.5, 3)

## interior knots will appear twice
MakeGrid(xd, 5, rm.dup = FALSE)

## interior knots only appear once
MakeGrid(xd, 5, rm.dup = TRUE)
}
