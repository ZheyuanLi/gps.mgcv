\name{RandomSpl}

\alias{RandomSpl}

\title{
A random cubic spline on [0, 1].
}

\description{
Generate a random cubic spline \eqn{g(x)} on [0, 1] and sample its values at \eqn{n} unevenly spaced \eqn{x}-values.
}

\usage{
RandomSpl(n, periodic = FALSE, plot = TRUE)
}

\arguments{
\item{n}{number of unevenly spaced \eqn{x}-values.}
\item{periodic}{TRUE to generate a periodic spline; FALSE to generate an ordinary spline.}
\item{plot}{TRUE to plot the generated random spline.}
}

\details{
A random spline is generated as a linear combination of 8 uniform B-splines, whose coefficients are a realization of an AR(1) process with coefficient -1/3. In case \eqn{g(x)} needs be periodic, these B-spline coefficients are further transformed to satisfy the periodic constraints, i.e., \eqn{g^{(m)}(0) = g^{(m)}(1)}, \eqn{m} = 0, 1, 2.
}

\value{
A list of the following components:
\itemize{
\item \code{x}, \eqn{n} unevenly spaced \eqn{x}-values drawn from a "tent" distribution, whose density peaks at \eqn{g(x)}'s local extrema;
\item \code{y}, \eqn{g(x)} evaluated at \code{x};
\item \code{xg}, 101 equidistant \eqn{x}-values (i.e., stepsize = 0.01);
\item \code{yg}, \eqn{g(x)} evaluated at \code{xg};
}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\references{
Li, Z. and Cao, J. (2022). General P-Splines for Non-Uniform B-Splines.
}

\examples{
require(gps.mgcv)

op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
spl.ordinary <- RandomSpl(200)
title("a random cubic spline")
spl.periodic <- RandomSpl(200, periodic = TRUE)
title("a random periodic cubic spline")
par(op)
}
