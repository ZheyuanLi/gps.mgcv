\name{DemoSpl}

\alias{DemoSpl}

\title{
Demonstrate polynomial spline and its B-spline representation
}

\description{
Demonstrate cubic spline and its B-spline representation.
}

\usage{
DemoSpl(uniform = TRUE)
}

\arguments{
\item{uniform}{TRUE for uniform knots and FALSE for non-uniform knots.}
}

\value{
A list giving domain knots, B-spline coefficients and piecewise polynomial coefficients.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps.mgcv)

## cubic spline with uniform knots
DemoSpl(uniform = TRUE)

## cubic spline with non-uniform knots
DemoSpl(uniform = FALSE)
}
