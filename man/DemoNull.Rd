\name{DemoNull}

\alias{DemoNull}

\title{
Demonstrate null space of a P-spline
}

\description{
Demonstrate the limiting behaviour of a P-spline at \eqn{\lambda = +\infty}. A cubic P-spline set up with non-uniform B-splines and 2nd order difference penalty is fitted to observations simulated from \eqn{y = x}. If this P-spline has the correct null space, its limiting fit will be a straight line regardless of knot locations. In this demo, non-uniform knots from different distributions (primarily Beta distributions with varying shape parameters) are attempted. The result shows that when handling non-uniform B-splines, standard P-spline has an incorrect and unpredictable limiting behavior that is sensitive to knot locations, whereas general P-spline has a correct and consistent limiting behavior.
}

\usage{
DemoNull(n, k, gps = FALSE)
}

\arguments{
\item{n}{number of simulated observations from \eqn{y = x}.}
\item{k}{number of interior knots to place.}
\item{gps}{TRUE to fit general P-spline; FALSE to fit standard P-spline.}
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps.mgcv)

## standard P-spline
DemoNull(n = 100, k = 10, gps = FALSE)

## general P-spline
DemoNull(n = 100, k = 10, gps = TRUE)
}
