# gps.mgcv (version 1.0)

**General P-Splines for Package 'mgcv'**

The package constructs and predicts the general P-splines of [Li and Cao (2022)](https://arxiv.org/abs/2201.06808) in **mgcv** by defining a new 'gps' smooth class (see `?gps.smooth`). A general P-spline *f(x)* is specified as `s(x, bs = 'gps', ...)` in a formula and estimated using **mgcv**'s model fitting functions, namely `gam` (generalized additive models, or GAMs), `bam` (GAMs for big data) and `gamm` (GAMs as mixed-effect models). General P-splines are state-of-the-art penalized B-splines. Unlike the standard P-splines of [Eilers and Marx (1996)](https://doi.org/10.1214/ss/1038425655) that only make sense for uniform B-splines on equidistant knots, they are properly defined for non-uniform B-splines on irregularly spaced knots, thanks to their powerful general difference penalty that accounts for uneven knot spacing. The package also contains functions for fitting and benchmarking different penalized B-splines (see `?Fit4BS` and `?SimStudy`) and a case study of smoothing Bone Mineral Content longitudinal data (see `?BMC`).

# Installation from CRAN

The package is not on CRAN yet, so you need to install it from this GitHub repository for now.

# Installation from GitHub

The package requires **R** (>= 4.0.0), **gps** (>= 1.1) and **mgcv** (>= 1.8-40).

Package **mgcv** comes with an **R** distribution, but you may need to update it if your existing version is older than the requirement. You can always install its latest version from CRAN using `install.packages("mgcv")`.

<!--Package **gps** is also [on CRAN](https://CRAN.R-project.org/package=gps), but [its GitHub version](https://github.com/Zheyuan/gps) may be more up to date. If both repositories have the same version-->

Install both **gps** and **gps.mgcv** from GitHub:

```r
## you may need to first install package 'devtools' from CRAN
devtools::install_github("ZheyuanLi/gps")
devtools::install_github("ZheyuanLi/gps.mgcv")
```

# Vignette

R code for [Li and Cao (2022)](https://arxiv.org/abs/2201.06808) is at: https://github.com/ZheyuanLi/gps-vignettes/blob/main/gps1.pdf
