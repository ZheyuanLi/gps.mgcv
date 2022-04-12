synBMC <- function (gender) {
  if (gender == "male") {
    with(synMBMC, data.frame(id = id, age = age, bmc = bmc, log.bmc = log(bmc)))
  } else if (gender == "female") {
    with(synFBMC, data.frame(id = id, age = age, bmc = bmc, log.bmc = log(bmc)))
  } else {
    stop("unknown 'gender'!")
  }
}

FitBMC <- function (gender, p, gps = TRUE, subset = NULL) {
  BMC <- synBMC(gender)
  if (p < 5) stop("p >= 5 required!")
  bs <- if (gps) "gps" else "ps"
  form <- sprintf("log.bmc ~ s(age, id, bs = 'fs', k = %d, m = c(%d, 1), xt = list(bs = '%s')) + s(age, bs = '%s', k = %d)",
                  p, 3 - !gps, bs, bs, p)
  form <- as.formula(form)
  suppressWarnings(mgcv::gamm(form, data = BMC, subset = subset))
}

PSE.serial <- function (n.sim = 100, gender, p, gps = TRUE) {
  BMC <- synBMC(gender)
  n <- nrow(BMC)
  ind <- split(1:n, BMC$id)
  id.age.min <- with(BMC, id[which.min(age)])
  id.age.max <- with(BMC, id[which.max(age)])
  ind[[id.age.min]] <- ind[[id.age.min]][-1]
  ind[[id.age.max]] <- ind[[id.age.max]][-length(ind[[id.age.max]])]
  ni <- unname(lengths(ind))
  ind <- ind[ni > 3]
  nSubjects <- length(ind)
  nTest <- round(0.1 * n)
  nf <- floor(nTest / nSubjects)
  nr <- nTest - nf * nSubjects
  PSE <- numeric(n.sim)
  for (i in 1:n.sim) {
    nd <- rep.int(nf, nSubjects)
    nd[sample.int(nSubjects, nr)] <- nf + 1
    TestInd <- sort(unlist(Map(sample, ind, nd), use.names = FALSE))
    TrainInd <- setdiff(1:n, TestInd)
    test <- BMC[TestInd, ]
    GAM <- FitBMC(gender, p, gps, TrainInd)
    log.bmc.pred <- predict(GAM$gam, test)
    PSE[i] <- mean((log.bmc.pred - test$log.bmc) ^ 2)
  }
  PSE
}

PSE.parallel <- function (gender, p, gps = TRUE, n.sim = 100, nc = 1) {
  if (nc == 1) {
    pse <- PSE.serial(n.sim, gender, p, gps)
  } else {
    n.sim.i <- n.sim / nc
    if (n.sim.i != round(n.sim.i)) {
      stop("n.sim is not a multiple of nc!")
    }
    n.sim.i <- rep.int(n.sim.i, nc)
    cl <- parallel::makeCluster(nc)
    pse <- parallel::clusterApply(cl, n.sim.i, PSE.serial, gender, p, gps)
    parallel::stopCluster(cl)
    pse <- unlist(pse)
  }
  cvm <- mean(pse)
  cve <- sd(pse) / sqrt(n.sim)
  c(cvm, cvm - 2 * cve, cvm + 2 * cve)
}

cvBMC <- function (gender, gps = TRUE, n.sim = 100, nc = 1) {
  if (!(gender %in% c("male", "female"))) stop("unknown 'gender'!")
  p <- 5:20
  pse <- matrix(0, length(p), 3, dimnames = list(p, c("mean", "lwr", "upr")))
  for (i in 1:length(p)) {
    cat(sprintf("-> doing p = %d\n", p[i]))
    pse[i, ] <- PSE.parallel(gender, p[i], gps, n.sim, nc)
  }
  pse
}

Fit4BS <- function (x, y, k, knots = NULL, domain = NULL, periodic = FALSE,
                    select = c("os", "sps", "nps", "gps")) {
  if (anyNA(match(select, c("os", "sps", "nps", "gps")))) {
    stop("invalid values in 'select'!")
  }
  os <- sps <- nps <- gps <- NULL
  degree <- 3
  d <- degree + 1
  p <- k + d
  Iknots <- knots
  if (is.null(Iknots)) {
    Aknots <- gps::PlaceKnots(x, d, k, domain, FALSE, FALSE)
    Nknots <- Aknots[seq.int(d, length(Aknots) - degree)]
    Iknots <- Nknots[seq.int(2, length.out = k)]
  } else if (length(Iknots) != k) {
    stop(sprintf("%d interior knots expected!", k))
  } else {
    x.min <- min(x)
    x.max <- max(x)
    if (is.null(domain)) {
      domain <- c(x.min, x.max)
    } else if (domain[1] > x.min || domain[2] < x.max) {
      stop("domain does not contain all x-values!")
    }
    if (min(Iknots) <= x.min || max(Iknots) >= x.max) {
      stop("interior knots must be strictly within range of x-values!")
    }
    Nknots <- c(domain[1], Iknots, domain[2])
    Aknots <- c(rep.int(domain[1], d), Iknots, rep.int(domain[2], d))
  }
  if ("gps" %in% select) {
    xt.arg <- list(domain = domain, periodic = periodic, derivative = FALSE)
    gam.form <- y ~ s(x, bs = 'gps', k = p, xt = xt.arg)
    gam.knots <- list(x = Iknots)
    gps <- mgcv::gam(gam.form, knots = gam.knots)
  }
  if ("os" %in% select) {
    xt.arg <- list(domain = domain, periodic = periodic, derivative = TRUE)
    gam.form <- y ~ s(x, bs = 'gps', k = p, xt = xt.arg)
    gam.knots <- list(x = Iknots)
    os <- mgcv::gam(gam.form, knots = gam.knots)
  }
  if (periodic) {
    pp <- p - degree
    if ("sps" %in% select) {
      Nknots.uniform <- gps::PlaceKnots(x, d, k, domain, TRUE, TRUE)
      gam.form <- y ~ s(x, bs = 'cp', k = pp)
      gam.knots <- list(x = Nknots.uniform)
      sps <- suppressWarnings(mgcv::gam(gam.form, knots = gam.knots))
    }
    if ("nps" %in% select) {
      gam.form <- y ~ s(x, bs = 'cp', k = pp)
      gam.knots <- list(x = Nknots)
      nps <- mgcv::gam(gam.form, knots = gam.knots)
    }
  } else {
    if ("sps" %in% select) {
      Aknots.uniform <- gps::PlaceKnots(x, d, k, domain, TRUE, FALSE)
      gam.form <- y ~ s(x, bs = 'ps', k = p)
      gam.knots <- list(x = Aknots.uniform)
      sps <- suppressWarnings(mgcv::gam(gam.form, knots = gam.knots))
    }
    if ("nps" %in% select) {
      gam.form <- y ~ s(x, bs = 'ps', k = p)
      gam.knots <- list(x = Aknots)
      nps <- mgcv::gam(gam.form, knots = gam.knots)
    }
  }
  list(os = os, sps = sps, nps = nps, gps = gps)
}

MSE4BS.serial <- function (n.sim = 100, x, g, k, knots = NULL, domain = NULL,
                           periodic = FALSE, n2s.ratio = 0.2,
                           select = c("os", "sps", "nps", "gps")) {
  n <- length(x)
  gx <- g(x)
  xg <- x
  yg <- gx
  newx <- data.frame(x = xg)
  mse <- matrix(0, n.sim, 4, dimnames = list(NULL, c("os", "sps", "nps", "gps")))
  for (i in 1:n.sim) {
    y <- rnorm(n, gx, n2s.ratio * sd(gx))
    fit <- Fit4BS(x, y, k, knots, domain, periodic, select)
    if (!is.null(fit$os)) {
      yh.os <- predict(fit$os, newdata = newx)
      mse[i, 1] <- c(base::crossprod(yh.os - yg)) / n
    }
    if (!is.null(fit$sps)) {
      yh.sps <- predict(fit$sps, newdata = newx)
      mse[i, 2] <- c(base::crossprod(yh.sps - yg)) / n
    }
    if (!is.null(fit$nps)) {
      yh.nps <- predict(fit$nps, newdata = newx)
      mse[i, 3] <- c(base::crossprod(yh.nps - yg)) / n
    }
    if (!is.null(fit$gps)) {
      yh.gps <- predict(fit$gps, newdata = newx)
      mse[i, 4] <- c(base::crossprod(yh.gps - yg)) / n
    }
  }
  mse[, select]
}

MSE4BS.parallel <- function (x, g, k, knots = NULL, domain = NULL, periodic = FALSE,
                             n2s.ratio = 0.2, n.sim = 100, nc = 2,
                             select = c("os", "sps", "nps", "gps")) {
  n.sim.i <- n.sim / nc
  if (n.sim.i != round(n.sim.i)) {
    stop("n.sim is not a multiple of nc!", call. = FALSE)
  }
  n.sim.i <- rep.int(n.sim.i, nc)
  cl <- parallel::makeCluster(nc)
  mse <- parallel::clusterApply(cl, n.sim.i, MSE4BS.serial, x, g, k, knots, domain, 
                                periodic, n2s.ratio, select)
  parallel::stopCluster(cl)
  do.call("rbind", mse)
}

MSE4BS <- function (x, g, k, knots = NULL, domain = NULL, periodic = FALSE,
                    n2s.ratio = 0.2, n.sim = 100, nc = 1,
                    select = c("os", "sps", "nps", "gps")) {
  if (anyNA(match(select, c("os", "sps", "nps", "gps")))) {
    stop("invalid values in 'select'!")
  }
  if (!is.function(g)) stop("g is not a function!")
  if (nc == 1) {
    MSE4BS.serial(n.sim, x, g, k, knots, domain, periodic, n2s.ratio, select)
  } else {
    MSE4BS.parallel(x, g, k, knots, domain, periodic, n2s.ratio, n.sim, nc, select)
  }
}

NullC_periodic <- function (xt, d) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  degree <- d - 1
  K <- length(xt)
  if (K < 2 * d + degree) {
    stop("length(xt) >= 3 * d - 1 required", call. = FALSE)
  }
  a <- xt[d]
  b <- xt[K - degree]
  local.knots.a <- xt[seq.int(from = 1, length.out = 2 * d)]
  local.knots.b <- xt[seq.int(to = K, length.out = 2 * d)]
  qs <- seq.int(0, degree - 1)
  Ca <- splines::splineDesign(local.knots.a, rep.int(a, degree), d, qs)
  Cb <- splines::splineDesign(local.knots.b, rep.int(b, degree), d, qs)
  Ca <- Ca[, seq.int(1, degree)]
  Cb <- Cb[, seq.int(2, d)]
  C <- cbind(Ca, -Cb, deparse.level = 0L)
  QR <- qr.default(t.default(C))
  Q <- qr.Q(QR, complete = TRUE)
  Q[, seq.int(d, 2 * degree), drop = FALSE]
}

BQ_periodic <- function (B, Q) {
  p <- ncol(B)
  degree <- ncol(Q)
  B.boundary <- B[, c(seq.int(1L, degree), seq.int(p - degree + 1L, p))]
  B.transformed <- B.boundary %*% Q
  B.unchanged <- B[, seq.int(degree + 1L, p - degree), drop = FALSE]
  cbind(B.unchanged, B.transformed)
}

smooth.construct.gps.smooth.spec <- function (object, data, knots) {
  if (length(object$term) != 1) {
    stop("can not construct multivariate splines!")
  }
  m <- object$p.order
  if (length(m) == 1 && is.na(m)) m <- c(3L, 2L)
  m <- as.integer(m)
  degree <- m[1]
  if (is.na(degree)) {
    stop("missing polynomial degree!")
  }
  if (degree < 1) {
    stop("can not fit piecewise constant!")
  }
  nDiff <- unique(sort(m[2:length(m)], decreasing = TRUE, na.last = TRUE))
  NumPen <- length(nDiff)
  nDiff.max <- nDiff[1]
  nDiff.min <- nDiff[NumPen]
  if (is.na(nDiff.min)) {
    stop("missing penalty order!")
  }
  if (nDiff.min <= 0 || nDiff.max > degree) {
    stop("penalty order must be between 1 and polynomial degree!")
  }
  object$m <- object$p.order <- m <- c(degree, nDiff)
  derivative <- object$xt$derivative
  if (is.null(derivative)) derivative <- FALSE
  else if (!is.logical(derivative) || is.na(derivative)) {
    stop("'xt$derivative' must be TRUE or FALSE!")
  }
  d <- degree + 1
  if (object$bs.dim < 0) object$bs.dim <- 2 * d
  p <- object$bs.dim
  x <- data[[object$term]]
  Iknots <- knots[[object$term]]
  domain <- object$xt$domain
  xu <- sort(unique.default(x))
  nx <- length(xu)
  periodic <- object$xt$periodic
  if (is.null(periodic)) periodic <- FALSE
  else if (!is.logical(periodic) || is.na(periodic)) {
    stop("'xt$periodic' must be TRUE or FALSE!")
  }
  object$xt <- list(derivative = derivative, domain = domain, periodic = periodic)
  min.p <- if (periodic) d + degree else d + 1
  if (nx < min.p) {
    stop("too few unique x-values to construct a spline!")
  }
  if (p < min.p || p > nx) {
    stop(sprintf("%d <= k <= %d required!", min.p, nx))
  }
  if (is.null(domain)) {
    a <- xu[1]
    b <- xu[nx]
  } else {
    a <- domain[1]
    b <- domain[2]
    if (xu[1] < a || xu[nx] > b) {
      stop("xt$domain does not contain all x-values!")
    }
  }
  k <- p - d
  if (is.null(Iknots)) {
    prob <- seq.int(0, 1, length.out = k + 2)
    Iknots <- quantile(xu, prob = prob, names = FALSE)[-c(1, k + 2)]
  } else if (length(Iknots) != k) {
    stop(sprintf("%d interior knots expected for basis dimension k = %d!", k, p))
  } else if (is.unsorted(Iknots, strictly = TRUE)) {
    stop("interior knots must be strictly ascending!")
  } else if (Iknots[1] <= xu[1] || Iknots[k] >= xu[nx]) {
    stop("interior knots must be strictly within range of x-values!")
  }
  Aknots <- c(rep.int(a, d), Iknots, rep.int(b, d))
  object$xt$domain <- c(a, b)
  object$interior.knots <- Iknots
  object$knots <- Aknots
  X <- splines::splineDesign(Aknots, x, d)
  D <- gps::SparseD(Aknots, d)[nDiff]
  if (derivative) {
    for (i in 1:NumPen) {
      S.bar <- gps::SandBar(Aknots, d, nDiff[i])
      U <- Matrix::chol(S.bar)
      D[[i]] <- U %*% D[[i]]
    }
  }
  if (periodic) {
    object$xt$Q <- Q <- NullC_periodic(Aknots, d)
    X <- BQ_periodic(X, Q)
    for (i in 1:NumPen) D[[i]] <- BQ_periodic(D[[i]], Q)
    p <- p - degree
  }
  S <- vector("list", NumPen)
  for (i in 1:NumPen) {
    S[[i]] <- as.matrix(Matrix::crossprod(D[[i]]))
    attr(S[[i]], "dimnames") <- NULL
  }
  names(S) <- paste0("order.", nDiff)
  object$X <- X
  object$S <- S
  if (periodic) {
    object$rank <- rep.int(p - 1, NumPen)
    object$null.space.dim <- 1
  } else {
    object$rank <- p - nDiff
    object$null.space.dim <- nDiff.min
  }
  structure(object, class = "gps.smooth")
}

Predict.matrix.gps.smooth <- function (object, data) {
  d <- object$m[1] + 1
  Aknots <- object$knots
  periodic <- object$xt$periodic
  x <- data[[object$term]]
  domain <- object$xt$domain
  a <- domain[1]
  b <- domain[2]
  if (periodic) {
    period <- b - a
    if (max(x) > b) {
      ind <- (x > b)
      x[ind] <- a + (x[ind] - b) %% period
    }
    if (min(x) < a) {
      ind <- (x < a)
      x[ind] <- b - (a - x[ind]) %% period
    }
  } else if (min(x) < a || max(x) > b) {
    stop("out-of-domain prediction is not possible for non-periodic splines!")
  }
  X <- splines::splineDesign(Aknots, x, d)
  if (periodic) X <- BQ_periodic(X, object$xt$Q)
  X
}

RandomSpl <- function (n, degree, periodic = FALSE, plot = TRUE) {
  if (degree < 1) stop("degree >= 1 required!")
  d <- degree + 1
  k <- 2 * d
  p <- k + d
  min.n <- 20 * p
  if (n < min.n) stop(sprintf("n >= %d required for this example!", min.n))
  h <- 1 / (k + 1)
  Nknots <- seq.int(0, 1, length.out = k + 2)
  laux <- seq.int(to = -h, by = h, length.out = degree)
  raux <- seq.int(from = 1 + h, by = h, length.out = degree)
  Aknots <- c(laux, Nknots, raux)
  K <- length(Aknots)
  xg <- seq.int(0, 1, by = 0.01)
  Bg <- splines::splineDesign(Aknots, xg, d, sparse = TRUE)
  if (periodic) {
    qs <- seq.int(0, degree - 1)
    Ca <- splines::splineDesign(Aknots, numeric(degree), d, qs)
    Cb <- splines::splineDesign(Aknots, rep.int(1, degree), d, qs)
    C <- Ca - Cb
    Q <- qr.Q(qr.default(t.default(C)))
  }
  repeat {
    b <- cumsum(rnorm(p))
    b <- b - mean.default(b)
    if (periodic) b <- b - c(Q %*% crossprod(Q, b))
    yg <- (Bg %*% b)@x
    ip <- which(diff.default(sign(diff.default(yg))) != 0) + 1L
    ip01 <- c(1L, ip, 101L)
    gap <- diff.default(ip01)
    if (min(gap) >= 10L && max(gap) < 100L) break
  }
  shift <- min(yg)
  yg <- yg - shift
  b <- b - shift
  xp <- xg[ip]
  yp <- yg[ip]
  np <- length(ip)
  x <- xtent(n, 0, 1, xp, h = 0.2)
  B <- splines::splineDesign(Aknots, x, d, sparse = TRUE)
  y <- (B %*% b)@x
  if (plot) {
    plot.default(xg, yg, type = "l", ann = FALSE)
    points.default(xp, yp, pch = 19)
    tent.x <- c(0, (xp[-np] + xp[-1]) / 2, 1, xp)
    tent.y <- c(numeric(np + 1), rep.int(0.2 * max(yg), np))
    ind <- order(tent.x)
    tent.x <- tent.x[ind]
    tent.y <- tent.y[ind]
    lines.default(tent.x, tent.y, lty = 2)
  }
  list(x = x, y = y, xg = xg, yg = yg)
}

SimStudy.serial <- function (n.sim = 100, n, degree, periodic = FALSE,
                             n2s.ratio = 0.2) {
  mse.os <- matrix(0, n.sim, degree)
  mse.sps <- matrix(0, n.sim, degree)
  mse.gps <- matrix(0, n.sim, degree)
  for (j in 1:n.sim) {
    spl <- RandomSpl(n, degree, periodic, plot = FALSE)
    x <- spl$x
    gx <- spl$y
    newx <- data.frame(x = spl$xg)
    sig <- n2s.ratio * sd(gx)
    sig2 <- sig * sig
    y <- rnorm(n, mean = gx, sd = sig)
    d <- degree + 1
    p <- 15 * d
    k <- p - d
    prob <- seq.int(0, 1, length.out = k + 2)
    Iknots <- quantile(x, prob = prob, names = FALSE)[-c(1, k + 2)]
    Nknots <- c(0, Iknots, 1)
    Aknots <- c(rep.int(0, d), Iknots, rep.int(1, d))
    for (m in 1:degree) {
      fit.os <- mgcv::gam(y ~ s(x, bs = "gps", k = p, m = c(degree, m),
                                xt = list(derivative = TRUE, periodic = periodic)),
                          knots = list(x = Iknots))
      fit.gps <- mgcv::gam(y ~ s(x, bs = "gps", k = p, m = c(degree, m),
                                 xt = list(derivative = FALSE, periodic = periodic)),
                           knots = list(x = Iknots))
      if (periodic) {
        pp <- p - degree
        fit.sps <- mgcv::gam(y ~ s(x, bs = "cp", k = pp, m = c(degree - 1, m)))
      } else {
        fit.sps <- mgcv::gam(y ~ s(x, bs = "ps", k = p, m = c(degree - 1, m)))
      }
      yh.os <- predict(fit.os, newdata = newx)
      yh.sps <- predict(fit.sps, newdata = newx)
      yh.gps <- predict(fit.gps, newdata = newx)
      mse.os[j, m] <- c(base::crossprod(yh.os - spl$yg)) / (101 * sig2)
      mse.sps[j, m] <- c(base::crossprod(yh.sps - spl$yg)) / (101 * sig2)
      mse.gps[j, m] <- c(base::crossprod(yh.gps - spl$yg)) / (101 * sig2)
    }
  }
  nm <- paste(rep.int(c("os", "sps", "gps"), degree),
              rep(1:degree, each = 3), sep = ".")
  matrix(rbind(mse.os, mse.sps, mse.gps), nrow = n.sim,
         dimnames = list(NULL, nm))
}

SimStudy.parallel <- function (n, degree, periodic = FALSE, n2s.ratio = 0.2,
                               n.sim = 100, nc = 2) {
  n.sim.i <- n.sim / nc
  if (n.sim.i != round(n.sim.i)) {
    stop("n.sim is not a multiple of nc!", call. = FALSE)
  }
  n.sim.i <- rep.int(n.sim.i, nc)
  cl <- parallel::makeCluster(nc)
  mse <- parallel::clusterApply(cl, n.sim.i, SimStudy.serial, n, degree,
                                periodic, n2s.ratio)
  parallel::stopCluster(cl)
  do.call("rbind", mse)
}

SimStudy <- function (n, degree, periodic = FALSE, n2s.ratio = 0.2,
                      n.sim = 100, nc = 1) {
  if (nc == 1) {
    SimStudy.serial(n.sim, n, degree, periodic, n2s.ratio)
  } else {
    SimStudy.parallel(n, degree, periodic, n2s.ratio, n.sim, nc)
  }
} 

xtent <- function (n, a, b, xp, h = 0) {
  np <- length(xp)
  breaks <- (xp[-length(xp)] + xp[-1]) / 2
  lp <- c(a, breaks)
  rp <- c(breaks, b)
  l1 <- xp - lp
  l2 <- rp - xp
  N <- round(c(rbind(l1, l2) / (b - a)) * (n - 1))
  N[1] <- (n - 1) - sum(N[-1])
  if (N[1] <= 0) stop("xtent failure!")
  x <- vector("list", 2 * np)
  for (i in 1:np) {
    h <- 0.2
    c1 <- 1 - h
    c2 <- 4 - 3 * h
    p1 <- seq.int(0, 0.5, length.out = N[2 * i - 1] + 1)[-1]
    p2 <- seq.int(0.5, 1, length.out = N[2 * i] + 1)[-1]
    q1 <- (sqrt(h * h + 8 * c1 * p1) - h) / (4 * c1)
    q2 <- (c2 - sqrt(c2 * c2 - 8 * c1 * (c1 + p2))) / (4 * c1)
    q1 <- 2 * q1
    q2 <- 2 * q2 - 1
    x[[2 * i - 1]] <- q1 * l1[i] + lp[i]
    x[[2 * i]] <- q2 * l2[i] + xp[i]
  }
  x <- c(a, unlist(x))
  x[n] <- b
  x
}

