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

Fit4BS <- function (x, y, k, m = 2, knots = NULL, domain = NULL, periodic = FALSE,
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
    Aknots <- PlaceKnots(x, d, k, domain, FALSE, FALSE)
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
    gam.form <- y ~ s(x, bs = 'gps', k = p, m = c(3, m), xt = xt.arg)
    gam.knots <- list(x = Iknots)
    gps <- mgcv::gam(gam.form, knots = gam.knots)
  }
  if ("os" %in% select) {
    xt.arg <- list(domain = domain, periodic = periodic, derivative = TRUE)
    gam.form <- y ~ s(x, bs = 'gps', k = p, m = c(3, m), xt = xt.arg)
    gam.knots <- list(x = Iknots)
    os <- mgcv::gam(gam.form, knots = gam.knots)
  }
  if (periodic) {
    pp <- p - degree
    if ("sps" %in% select) {
      Nknots.uniform <- PlaceKnots(x, d, k, domain, TRUE, TRUE)
      gam.form <- y ~ s(x, bs = 'cp', k = pp, m = c(2, m))
      gam.knots <- list(x = Nknots.uniform)
      sps <- suppressWarnings(mgcv::gam(gam.form, knots = gam.knots))
    }
    if ("nps" %in% select) {
      gam.form <- y ~ s(x, bs = 'cp', k = pp, m = c(2, m))
      gam.knots <- list(x = Nknots)
      nps <- mgcv::gam(gam.form, knots = gam.knots)
    }
  } else {
    if ("sps" %in% select) {
      Aknots.uniform <- PlaceKnots(x, d, k, domain, TRUE, FALSE)
      gam.form <- y ~ s(x, bs = 'ps', k = p, m = c(2, m))
      gam.knots <- list(x = Aknots.uniform)
      sps <- suppressWarnings(mgcv::gam(gam.form, knots = gam.knots))
    }
    if ("nps" %in% select) {
      gam.form <- y ~ s(x, bs = 'ps', k = p, m = c(2, m))
      gam.knots <- list(x = Aknots)
      nps <- mgcv::gam(gam.form, knots = gam.knots)
    }
  }
  list(os = os, sps = sps, nps = nps, gps = gps)
}

MSE4BS.serial <- function (n.sim = 100, x, g, k, m = 2, knots = NULL, domain = NULL,
                           periodic = FALSE, n2s.ratio = 0.2,
                           select = c("os", "sps", "nps", "gps")) {
  n <- length(x)
  gx <- g(x)
  xg <- seq.int(min(x), max(x), length.out = n)  ## xg <- x
  yg <- g(xg)  ## yg <- gx
  newx <- data.frame(x = xg)
  mse <- matrix(0, n.sim, 4, dimnames = list(NULL, c("os", "sps", "nps", "gps")))
  for (i in 1:n.sim) {
    y <- rnorm(n, gx, n2s.ratio * sd(gx))
    fit <- Fit4BS(x, y, k, m, knots, domain, periodic, select)
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

MSE4BS.parallel <- function (x, g, k, m = 2, knots = NULL, domain = NULL,
                             periodic = FALSE, n2s.ratio = 0.2, n.sim = 100, nc = 2,
                             select = c("os", "sps", "nps", "gps")) {
  n.sim.i <- n.sim / nc
  if (n.sim.i != round(n.sim.i)) {
    stop("n.sim is not a multiple of nc!", call. = FALSE)
  }
  n.sim.i <- rep.int(n.sim.i, nc)
  cl <- parallel::makeCluster(nc)
  mse <- parallel::clusterApply(cl, n.sim.i, MSE4BS.serial, x, g, k, m, knots, domain, 
                                periodic, n2s.ratio, select)
  parallel::stopCluster(cl)
  do.call("rbind", mse)
}

MSE4BS <- function (x, g, k, m = 2, knots = NULL, domain = NULL, periodic = FALSE,
                    n2s.ratio = 0.2, n.sim = 100, nc = 1,
                    select = c("os", "sps", "nps", "gps")) {
  if (anyNA(match(select, c("os", "sps", "nps", "gps")))) {
    stop("invalid values in 'select'!")
  }
  if (!is.function(g)) stop("g is not a function!")
  if (nc == 1) {
    MSE4BS.serial(n.sim, x, g, k, m, knots, domain, periodic, n2s.ratio, select)
  } else {
    MSE4BS.parallel(x, g, k, m, knots, domain, periodic, n2s.ratio, n.sim, nc, select)
  }
}


PlaceKnots <- function (x, d, k, domain = NULL, uniform = FALSE, periodic = FALSE) {
  xu <- sort.int(unique.default(x))
  nx <- length(xu)
  if (is.null(domain)) {
    a <- xu[1L]
    b <- xu[nx]
    domain <- c(a, b)
  } else {
    a <- domain[1L]
    b <- domain[2L]
    if (xu[1L] < a || xu[nx] > b) {
      stop("'domain' does not contain all x-values!")
    }
  }
  degree <- d - 1
  if (uniform) {
    xd <- seq.int(a, b, length.out = k + 2)
    h <- xd[2L] - xd[1L]
    laux <- seq.int(to = a - h, by = h, length.out = degree)
    raux <- seq.int(from = b + h, by = h, length.out = degree)
  } else {
    prob <- seq.int(0, 1, length.out = k + 2)
    xd <- quantile(xu, prob, names = FALSE)
    xd[c(1, k + 2)] <- domain
    laux <- rep.int(a, degree)
    raux <- rep.int(b, degree)
  }
  if (periodic) xd else c(laux, xd, raux)
}

IsAscending <- function (x, n = length(x), xi = 1L) {
  if (!is.double(x)) stop("'x' is not in double-precision mode!")
  .Call("C_IsAscending", x, n, xi, PACKAGE = "gps.mgcv") > 0L
}

MonoKnots <- function (xt, d) {
  xt <- as.double(xt)
  K <- length(xt)
  if (K < 2 * d) stop("length(xt) >= 2 * d required!", call. = FALSE)
  flag <- IsAscending(xt, n = K - 2 * (d - 1), xi = d)
  if (!flag) stop("Domain knots are not strictly ascending!", call. = FALSE)
  lbnd <- xt[seq.int(from = 1L, length.out = d)]
  rbnd <- xt[seq.int(to = K, length.out = d)]
  if (is.unsorted(lbnd) || is.unsorted(rbnd)) {
    stop("Boundary knots are not non-decreasing!", call. = FALSE)
  }
  xt
}

Diff <- function (x, k = 1L, n = length(x), xi = 1L) {
  if (!is.double(x)) stop("'x' is not in double-precision mode!")
  .Call("C_Diff", x, k, n, xi, PACKAGE = "gps.mgcv")
}

SparseWtDelta <- function (h) {
  r <- length(h)
  x <- rep.int(c(-1, 1), r) * rep(1 / h, each = 2)
  i <- rep(seq.int(0L, r - 1L), each = 2)
  p <- c(0L, seq.int(1L, length.out = r, by = 2L), 2L * r)
  methods::new("dgCMatrix", i = i, p = p, Dim = c(r, r + 1L), x = x)
}

SparseD <- function (xt, d) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  K <- length(xt)
  D <- vector("list", d - 1)
  h <- Diff(x = xt, k = d - 1, n = K - 2, xi = 2)
  D[[1]] <- SparseWtDelta(h)
  i <- 2
  while (i < d) {
    h <- Diff(x = xt, k = d - i, n = K - 2 * i, xi = i + 1)
    D[[i]] <- SparseWtDelta(h) %*% D[[i - 1]]
    i <- i + 1
  }
  D
}

MakeGrid <- function (xd, n, rm.dup = FALSE) {
  xd <- as.double(xd)
  if (!IsAscending(xd)) stop("'xd' is not strictly ascending!")
  if (n == 1) {
    lp <- xd[-length(xd)]
    rp <- xd[-1]
    mp <- 0.5 * (lp + rp)
    return(mp)
  }
  if (rm.dup && (n == 2)) return(xd)
  .Call("C_MakeGrid", xd, n, rm.dup, PACKAGE = "gps.mgcv")
}

QuadWts <- function (ord) {
  if (ord == 1) {
    return(2)
  }
  p <- seq.int(0, ord - 1)
  P <- outer(seq.int(-1, 1, length.out = ord), p, "^")
  Pinv <- solve.default(P)
  pow <- outer(1:ord, p, "+")
  H <- (1 - (-1) ^ pow) / pow
  base::crossprod(Pinv, H %*% Pinv)
}

SbarBlocks <- function (xt, d, m) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  if (m < 0 || m >= d) stop("0 <= m <= d - 1 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  K <- length(xt)
  ord <- d - m
  if (ord == 1) {
    h <- Diff(x = xt, n = K - 2 * (d - 1), xi = d)
    return(h)
  }
  W <- QuadWts(ord)
  xd <- xt[seq.int(d, K - d + 1)]
  xg <- MakeGrid(xd, ord)
  xt.local <- xt[seq.int(1 + m, K - m)]
  B <- splines::splineDesign(xt.local, xg, ord, sparse = TRUE)
  .Call("C_SbarBlocks", xd, W, B@x, PACKAGE = "gps.mgcv")
}

LTB2Csparse <- function (L, k = n, symmetric = FALSE) {
  b1 <- nrow(L); n <- ncol(L); b <- b1 - 1L
  k <- as.integer(k)
  if (symmetric) {
    k <- n
  } else if (k < n || k > n + b) {
    stop(sprintf("%d <= k <= %d expected!", n, n + b))
  }
  r <- n + b - k
  p <- rep.int(seq.int(b1, by = -1L, length.out = r + 1L), c(n - r, rep.int(1, r)))
  i <- sequence.default(nvec = p, from = seq.int(0L, n - 1L))
  p <- c(0L, cumsum(p))
  if (r) {
    x <- numeric(p[length(p)])
    x[seq_len(b1 * (k - b))] <- L[, seq_len(k - b)]
    columns <- seq.int(to = n, length.out = r)
    for (j in columns) {
      ind <- seq.int(p[j] + 1L, p[j + 1L])
      x[ind] <- L[seq_len(length(ind)), j]
    }
  } else {
    x <- c(L)
  }
  if (symmetric) {
    methods::new("dsCMatrix", i = i, p = p, x = x, Dim = c(k, n), uplo = "L")
  } else {
    methods::new("dgCMatrix", i = i, p = p, x = x, Dim = c(k, n))
  }
}

SparseSbar <- function (xt, d, m, do.chol = FALSE) {
  blocks <- SbarBlocks(xt, d, m)
  if (d - m == 1) {
    if (do.chol) blocks <- sqrt(blocks)
    return(blocks)
  }
  LTB <- .Call("C_SbarLTB", blocks, do.chol, PACKAGE = "gps.mgcv")
  LTB2Csparse(LTB, symmetric = !do.chol)
}

SandBar <- function (xt, d, m) {
  Sbar <- SparseSbar(xt, d, m)
  if (d - m == 1) {
    q <- length(Sbar)
    Sbar <- methods::new("ddiMatrix", diag = "N", x = Sbar, Dim = c(q, q))
  }
  Sbar
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
  D <- SparseD(Aknots, d)[nDiff]
  if (derivative) {
    for (i in 1:NumPen) {
      S.bar <- SandBar(Aknots, d, nDiff[i])
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

ComputeLD <- function (xt, d) {
  if (d < 2) stop("d >= 2 required!", call. = FALSE)
  xt <- MonoKnots(xt, d)
  .Call("C_ComputeLD", xt, d, PACKAGE = "gps.mgcv")
}

OrthNullD <- function (ld, m = 1, orthogonal = TRUE) {
  if (m < 1 || m > ncol(ld)) stop("1 <= m <= d - 1 required!", call. = FALSE)
  basis <- .Call("C_NullD", ld, m, PACKAGE = "gps.mgcv")
  if (orthogonal && (m > 1)) {
    Q <- qr.Q(qr.default(basis[, m:1]))[, m:1]
    i <- sequence.default(1:(m - 1))
    j <- rep.int(2:m, 1:(m - 1))
    Q[cbind(i, j)] <- 0
    basis <- Q
  }
  basis
}

NullD <- function (xt, d, m) {
  ld <- ComputeLD(xt, d)
  OrthNullD(ld, m, orthogonal = FALSE)
}

as_matrix <- function (A) {
  if (is.matrix(A)) return(A)
  sparse <- inherits(A, "dsparseMatrix")
  dense <- inherits(A, "ddenseMatrix")
  if (!sparse && !dense) {
    stop("'A' is not a \"dsparseMatrix\" or \"ddenseMatrix\"!")
  }
  nnz <- length(A@x)
  nr <- A@Dim[1]
  nc <- A@Dim[2]
  if (nnz == nr * nc) {
    denseA <- matrix(A@x, nr, nc)
  } else if (inherits(A, "dCsparseMatrix")) {
    i <- A@i
    j <- rep.int(seq.int(0L, nc - 1L), diff.default(A@p))
    denseA <- matrix(0, nr, nc)
    denseA[j * nr + (i + 1L)] <- A@x
    if (inherits(A, "dsCMatrix")) denseA[i * nc + (j + 1L)] <- A@x
  } else {
    stop("Not implemented yet!")
  }
  denseA
}

PlotNull <- function (x, y, k, shape1 = 3, shape2 = 3, gps = FALSE) {
  d <- 4
  m <- 2
  if (shape1 < 1 && shape2 < 1) {
    xd <- seq.int(-1/8, 1/8, length.out = k + 2)
    xd <- sign(xd) * abs(xd) ^ (1 / 3) + 0.5
    distr <- "U-quadratic(0, 1)"
  } else {
    xd <- qbeta(seq.int(0, 1, length.out = k + 2), shape1, shape2)
    distr <- sprintf("Beta(%d, %d)", shape1, shape2)
  }
  xt <- c(numeric(d - 1), xd, rep.int(1, d - 1))
  K <- length(xt)
  xg <- seq.int(0, 1, length.out = 101)
  B <- splines::splineDesign(xt, x, d, sparse = TRUE)
  Bg <- splines::splineDesign(xt, xg, d, sparse = TRUE)
  if (gps) {
    H <- NullD(xt, d, m)
  } else {
    H <- NullD(seq.int(0, 1, length.out = K), d, m)
  }
  X <- as_matrix(B %*% H)
  Xg <- as_matrix(Bg %*% H)
  XtX <- base::crossprod(X)
  Xty <- base::crossprod(X, y)
  U <- chol.default(XtX)
  f <- forwardsolve(U, Xty, upper.tri = TRUE, transpose = TRUE)
  b <- backsolve(U, f)
  yg <- c(Xg %*% b)
  ylim <- range(y, yg)
  plot.default(x, y, col = 8, ylim = ylim, ann = FALSE, xaxt = "n", yaxt = "n")
  abline(v = xd, lty = 2, col = 8)
  lines.default(xg, yg, lwd = 2)
  title(distr)
}

DemoNull <- function (n, k, gps = FALSE) {
  x <- seq.int(0, 1, length.out = n)
  y <- rnorm(n, mean = x, sd = 0.2 * sd(x))
  op <- par(mfrow = c(2, 3), mar = c(0.25, 0.25, 1.5, 0.25), oma = c(0, 0, 1.5, 0))
  on.exit(par(op))
  PlotNull(x, y, k, 3, 5, gps)
  PlotNull(x, y, k, 5, 5, gps)
  PlotNull(x, y, k, 5, 3, gps)
  PlotNull(x, y, k, 1, 3, gps)
  PlotNull(x, y, k, 0.5, 0.5, gps)
  PlotNull(x, y, k, 3, 1, gps)
  label <- sprintf("limiting %s P-spline fit at ", c("standard", "general")[gps + 1L])
  expr <- substitute(expression(paste(label, lambda == +infinity)), list(label = label))
  title(eval(expr), outer = TRUE)
}

Zero2NA <- function (Bsparse) {
  if (!inherits(Bsparse, "dgCMatrix")) {
    stop("'Bsparse' is not a \"dgCMatrix\"!")
  }
  B <- matrix(NA_real_, Bsparse@Dim[1], Bsparse@Dim[2])
  i <- Bsparse@i + 1L
  j <- rep.int(seq_len(Bsparse@Dim[2]), diff.default(Bsparse@p))
  B[cbind(i, j)] <- Bsparse@x
  B
}

DemoSpl <- function (uniform = TRUE) {
  if (uniform) {
    xd <- 1:6
    xt <- seq.int(-2, 9)
    xg <- MakeGrid(xt, 21)
    b <- c(0.44, 1.11, 1.66, 0.25, 1.6, 1.43, 1.49, 2.52)
  } else {
    xd <- c(1, 2.06, 2.65, 3.22, 3.91, 6)
    xt <- xd[c(1, 1, 1, 1, 2:5, 6, 6, 6, 6)]
    xg <- MakeGrid(xd, 21)
    b <- c(1.9, 1.84, 1.35, 1.67, 0.48, 1.85, 1.59, 1.35)
  }
  Bg <- splines::splineDesign(xt, xg, outer.ok = uniform, sparse = TRUE)
  Bgdense <- Zero2NA(Bg)
  yg <- (Bg %*% b)@x
  if (uniform) {
    i1 <- 21 * 3 + 1
    i2 <- 21 * 8
    x <- xg[i1:i2]
    y <- yg[i1:i2]
  } else {
    x <- xg
    y <- yg
  }
  ymax <- max(y)
  yd <- y[c(1, 1:5 * 21)]
  op <- par(xpd = TRUE, mar = c(2, 2, 1.5, 0.5))
  on.exit(par(op))
  plot.default(x, y, type = "n", xlim = c(-2, 9), ylim = c(0, ymax),
               xaxs = "i", yaxs = "i", ann = FALSE, bty = "n")
  polygon(c(1, 6, 6, 1), c(0, 0, ymax, ymax), col = gray(0.3, 0.2), border = NA)
  lines.default(x, y, lwd = 2)
  points.default(xd, yd, pch = c(17, 19, 19, 19, 19, 17))
  matlines(xg, Bgdense, type = "l", lty = 1, lwd = 2, col = c(1:6, 8, 7))
  points.default(xd[2:5], numeric(4), pch = 19)
  if (uniform) {
    points.default(c(-2:0, 7:9), numeric(6), pch = 15)
  } else {
    segments(-2, 0, 1, 0, lwd = 2, col = 1)
    segments(6, 0, 9, 0, lwd = 2, col = 7)
    segments(1, 1, 1, 0, lwd = 2, col = 1, lty = 2)
    segments(6, 1, 6, 0, lwd = 2, col = 7, lty = 2)
    ystack <- seq.int(0.1, by = 0.1, length.out = 3)
    points.default(c(1, 1, 1, 6, 6, 6), c(ystack, ystack), pch = 15)
  }
  points.default(c(1, 6), numeric(2), pch = 17)
  type <- c("non-uniform", "uniform")[uniform + 1L]
  title(sprintf("cubic spline with %s knots", type))
  dB <- splines::splineDesign(xt, rep.int(xd[1:5], 3), derivs = rep(1:3, each = 5))
  pc <- c(dB %*% b) * rep(1 / factorial(1:3), each = 5)
  lab1 <- paste0("f", 1:5)
  lab2 <- c("constant", "linear", "quadratic", "cubic")
  pc <- matrix(c(yd[1:5], pc), ncol = 4, dimnames = list(lab1, lab2))
  list(xd = xd, bspl.coef = b, piecepoly.coef = pc)
}

RandomSpl <- function (n, periodic = FALSE, plot = TRUE) {
  p <- 8
  if (n < 200) stop("n >= 200 required!")
  Aknots <- c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)
  xg <- seq.int(0, 1, by = 0.01)
  Bg <- splines::splineDesign(Aknots, xg, sparse = TRUE)
  dBg <- splines::splineDesign(Aknots, xg, derivs = 1, sparse = TRUE)
  if (periodic) {
    qs <- seq.int(0, 2)
    Ca <- splines::splineDesign(Aknots, numeric(3), derivs = qs)
    Cb <- splines::splineDesign(Aknots, rep.int(1, 3), derivs = qs)
    C <- Ca - Cb
    Q <- qr.Q(qr.default(t.default(C)))
  }
  repeat {
    b <- c(arima.sim(list(ar = -1/3), n = p))
    if (periodic) b <- b - c(Q %*% crossprod(Q, b))
    yg <- (Bg %*% b)@x
    dyg <- (dBg %*% b)@x
    ip.yg <- which(diff.default(sign(diff.default(yg))) != 0) + 1L
    ip.dyg <- which(diff.default(sign(diff.default(dyg))) != 0) + 1L
    ip <- sort(c(ip.yg, ip.dyg))
    x.gap <- diff.default(c(1L, ip, 101L))
    ye <- yg[c(1L, ip.yg, 101L)]
    y.gap <- abs(diff.default(ye))
    cond1 <- length(ip.yg) >= 2  ## 2+ local extrema wanted
    cond2 <- min(x.gap) > 5L  ## not too close in x-direction
    cond3 <- min(y.gap) > 0.2 * (max(ye) - min(ye))  ## not too close in y-direction
    if (cond1 && cond2 && cond3) break
  }
  shift <- min(ye)
  yg <- yg - shift
  b <- b - shift
  xp <- xg[ip]
  yp <- yg[ip]
  np <- length(ip)
  x <- xtent(n, 0, 1, xp, h = 0.1)
  B <- splines::splineDesign(Aknots, x, sparse = TRUE)
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
    spl <- RandomSpl(n, periodic, plot = FALSE)
    x <- spl$x
    gx <- spl$y
    newx <- data.frame(x = spl$xg)
    sig <- n2s.ratio * sd(gx)
    sig2 <- sig * sig
    y <- rnorm(n, mean = gx, sd = sig)
    d <- degree + 1
    p <- 40
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

.onAttach <- function (libname, pkgname) {
  packageStartupMessage(sprintf("--------------------\n< Urgent Notice on 2024-09-14 >\nPackage gps.mgcv no longer depends on package gps!\nThe up-to-date R code for our paper is available at:\n%s/gps.mgcv/code.html\n--------------------", .libPaths()[1]))
}

