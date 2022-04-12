\name{Fit4BS}

\alias{Fit4BS}
\alias{MSE4BS}

\title{
Experiment 4 variants of penalized B-splines
}

\description{
Fit 4 variants of penalized B-splines to noisy \eqn{(x, y)} data simulated from \eqn{g(x)} and compare their MSE performance under repeated simulations.
}

\usage{
Fit4BS(x, y, k, knots = NULL, domain = NULL, periodic = FALSE,
       select = c("os", "sps", "nps", "gps"))

MSE4BS(x, g, k, knots = NULL, domain = NULL, periodic = FALSE,
       n2s.ratio = 0.2, n.sim = 100, nc = 1,
       select = c("os", "sps", "nps", "gps"))
}

\arguments{
\item{x, y}{\eqn{(x, y)} data.}
\item{g}{the true generative function \eqn{g(x)} from which data are simulated.}
\item{k}{number of interior knots to place. This gives \eqn{k + 4} ordinary cubic B-splines or \eqn{k + 1} periodic cubic B-splines.}
\item{knots}{a sequence of interior knots on which non-uniform B-splines are constructed for "os", "nps" and "gps" (see Details). Will be automatically placed at equal quantiles of \eqn{x}-values if not specified. If specified, must match \code{k}.}
\item{domain}{a vector of two values giving domain interval \eqn{[a, b]}. Will use \code{min(x)} and \code{max(x)} if not specified.}
\item{periodic}{TRUE to fit a periodic spline; FALSE to fit an ordinary spline.}
\item{n2s.ratio}{a noise-to-signal ratio that is used when simulating noisy \eqn{y}-values from \eqn{g(x)}.}
\item{n.sim}{number of simulations.}
\item{nc}{number of CPU cores to use for parallel simulations.}
\item{select}{a character vector giving names (see Details) of the variants to fit.}
}

\details{
The 4 variants of penalized B-splines are:
\itemize{
\item O-spline ("os"): non-uniform B-splines + derivative penalty;
\item standard P-spline ("sps"): uniform B-splines + standard difference penalty;
\item naive P-spline ("nps"): non-uniform B-splines + standard difference penalty;
\item general P-spline ("gps"): non-uniform B-splines + general difference penalty.
}
They are constructed with cubic B-splines on \eqn{k} interior knots, penalized with a 2nd order derivative or difference penalty. By default, quantile knots are placed automatically for "os", "nps" and "gps". But this can be overridden by user-specified interior knots through argument \code{knots}.
}

\value{
A list of 4 fitted "gam" models, with components "os", "sps", "nps" and "gps". A component may be empty (with value \code{NULL}) if it is excluded from \code{select}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\references{
Li, Z. and Cao, J. (2022). General P-Splines for Non-Uniform B-Splines.
}

\examples{
\dontrun{

require(mgcv)
require(gps.mgcv)

## a U-shaped curve
g <- function (x) (1 / 8) * abs(x) ^ 3

## number of data
n <- 400
## noise-to-signal ratio
n2s.ratio <- 0.4
## number of interior knots
k <- 40

## (x, y) data
x <- qnorm(seq.int(pnorm(-3), pnorm(3), length.out = n))
gx <- g(x)
y <- rnorm(length(gx), gx, n2s.ratio * sd(gx))

## fit 4 variants of penalized B-splines
fit <- Fit4BS(x, y, k)

## plot a fitted penalized B-spline against observed (x, y) data and true function g(x)
## gamfit: a fitted GAM model, should only have a single s(x) term of B-spline class
## x, y: observed (x, y) data
## g: the true function g(x)
PlotXYFit <- function (gamfit, x, y, g) {
  ## g must be a function
  if (!is.function(g)) stop("g is not a function!")
  ## plot (x, y) data
  plot(x, y, col = 8, ann = FALSE)
  ## observed x-values may not be evenly spaced
  ## to plot g(x), we'd better use evenly spaced x-values
  nx <- length(x)
  xg <- seq.int(min(x), max(x), length.out = nx)
  lines(xg, g(xg), col = 2, lwd = 2)
  ## show knot location
  abline(v = gamfit$smooth[[1]]$knots, lty = 2, col = 8)
  ## plot estimated f(x) at xg
  yh <- predict(gamfit, newdata = data.frame(x = xg))
  lines(xg, yh, lwd = 2)
}

## compare 4 variants of penalized B-splines
op <- par(no.readonly = TRUE)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2))
par(mar = c(2, 2, 1.5, 0.5))
## O-spline
PlotXYFit(fit$os, x, y, g)
title("O-spline")
## standard P-spline
PlotXYFit(fit$sps, x, y, g)
title("standard P-spline")
## naive P-spline
PlotXYFit(fit$nps, x, y, g)
title("naive P-spline")
## general P-spline
PlotXYFit(fit$gps, x, y, g)
title("general P-spline")
## boxplot of MSE
mse <- MSE4BS(x, g, k = k, n2s.ratio = n2s.ratio, n.sim = 100, nc = 1)
boxplot(mse, main = "MSE (100 simulations)")
par(op)
}
}
