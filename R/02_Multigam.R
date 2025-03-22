# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multigam Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Multigam",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Gamma Distribution
#' @name Multigam
#'
#' @param x an object of class `Multigam`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Multigam`.
#' @param shape,scale numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
#' @param log logical. Should the log of the density be returned?
#' @param n numeric. The sample size.
#'
#' @inherit Distributions return
#'
#' @export
Multigam <- function(shape = 1, scale = 1) {
  new("Multigam", shape = shape, scale = scale)
}

setValidity("Multigam", function(object) {
  if(any(object@shape <= 0)) {
    stop("shape has to be a vector with positive elements")
  }
  if(length(object@scale) != 1) {
    stop("scale has to be a numeric of length 1")
  }
  if(object@scale <= 0) {
    stop("scale has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
#' @export
dmultigam <- function(x, shape, scale, log = FALSE) {

  if (length(x) != length(shape)) {
    stop("The lengths of x (", length(x), ") and shape (",
         length(shape), ") must be equal.")
  }

  z <- fd(x)
  xk <- x[length(x)]
  a0 <- sum(shape)

  ld <- sum(shape * log(z)) - a0 * log(scale) - sum(lgamma(shape)) - xk / scale

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname Multigam
#' @export
rmultigam <- function(n, shape, scale) {

  k <- length(shape)
  x <- matrix(nrow = n, ncol = k)
  for (j in 1:k) {
    x[, j] <- stats::rgamma(n, shape[j], scale = scale)
  }

  t(apply(x, 1, cumsum))

}

#' @rdname Multigam
setMethod("d", signature = c(x = "Multigam"),
          function(x) {
            function(y, log = FALSE) {
              dmultigam(y, shape = x@shape, scale = x@scale, log = log)
            }
          })

#' @rdname Multigam
setMethod("r", signature = c(x = "Multigam"),
          function(x) {
            function(n) {
              rmultigam(n, shape = x@shape, scale = x@scale)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ATTENTION: some moments were copy-pasted from Gam
# mean is ok

#' @rdname Multigam
setMethod("mean",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  cumsum(x@shape * x@scale)

})

#' @rdname Multigam
setMethod("var",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  x@shape * x@scale ^ 2

})

#' @rdname Multigam
setMethod("entro",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  a <- x@shape
  a + log(x@scale) + lgamma(a) + (1 - a) * digamma(a)

})

#' @rdname Multigam
setMethod("finf",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  # Preliminaries
  a <- x@shape
  b <- x@scale
  k <- length(a)
  a0 <- sum(a)

  D <- rbind(cbind(diag(trigamma(a)),
                   matrix(1 / b, nrow = k, ncol = 1)),
             c(rep(1 / b, k), a0 / b ^ 2))

  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llmultigam <- function(x, shape, scale) {
  ll(x, prm = c(shape, scale), distr = Multigam())
}

#' @rdname Multigam
setMethod("ll",
          signature  = c(x = "matrix", prm = "numeric", distr = "Multigam"),
          definition = function(x, prm, distr) {

  k <- length(prm)
  sum(apply(x, MARGIN = 1, FUN = dmultigam,
            shape = prm[1:(k - 1)], scale = prm[k], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Multigam"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz - log(b))

  - sum(a) * log(b) - sum(lgamma(a)) - xk / b + sum((a - 1) * logz)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Multigam"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz -log(b))

  db <- - xk / par ^ 2
  da <- 1 / (par * trigamma(a))

  - sum(da) * log(b) - sum(a) * db / b - sum(digamma(a) * da) +
    xk * db / b ^ 2 + sum(logz * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
emultigam <- function(x, type = "mle", ...) {

  estim(x, Multigam(), type, ...)

}

#' @rdname Multigam
setMethod("mle",
          signature  = c(x = "matrix", distr = "Multigam"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  k <- ncol(x)
  logz <- colMeans(log(fd(x)))
  xk <- mean(x[, k])
  tx <- c(logz, xk)

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))[1:k]),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  b <- xk / par
  a <- idigamma(logz - log(b))

  par <- c(a, b)

  names(par) <- c(paste0("shape", 1:k), "scale")
  par

})

#' @rdname Multigam
setMethod("me",
          signature  = c(x = "matrix", distr = "Multigam"),
          definition = function(x, distr) {

  z <- fd(x)
  mz <- colMeans(z)
  scale <- mean(colVar(z) / mz)
  shape <- mz / scale

  c(shape = shape, scale = scale)

})

me2 <- function(x, distr) {

  w <- gendir(x)
  shape <- unname(me(w, Dir()))
  xk <- mean(x[, ncol(x)])
  scale <- xk / sum(shape)

  c(shape = shape, scale = scale)

}

#' @rdname Multigam
setMethod("same",
          signature  = c(x = "matrix", distr = "Multigam"),
          definition = function(x, distr) {

  z <- fd(x)
  scale <- mean(diag(stats::cov(z, log(z))))
  shape <- colMeans(z) / scale

  c(shape = shape, scale = scale)

})

same2 <- function(x, distr) {

  w <- gendir(x)
  shape <- unname(same(w, Dir()))
  xk <- mean(x[, ncol(x)])
  scale <- xk / sum(shape)

  c(shape = shape, scale = scale)

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vmultigam <- function(shape, scale, type = "mle") {

  avar(Multigam(shape = shape, scale = scale), type = type)

}

#' @rdname Multigam
setMethod("avar_mle",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  trinv <- 1 / trigamma(a)
  cons <- a0 - sum(trinv)

  D <- diag(cons * trinv) + Matrix(trinv, k, 1) %*% Matrix(trinv, 1, k)
  D <- cbind(D, - b * trinv)
  D <- rbind(D, t(c(- b * trinv, b ^ 2))) / cons

  D <- as.matrix(nearPD(D))
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

#' @rdname Multigam
setMethod("avar_me",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A11 <- (Matrix(a, k, 1) %*% Matrix(2 + 1 / a, 1, k) / k + diag(1, k, k)) / b
  A12 <- - Matrix(a, k, 1) %*% Matrix(1 / a, 1, k) / (k * b ^ 2)
  A21 <- - (2 + 1 / a) / k
  A22 <- 1 / (a * k * b)
  A <- rbind(cbind(A11, A12), Matrix(c(A21, A22), 1, 2 * k))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- 2 * a * (a + 1) * b ^ 4 * ( 2 * a + 3)
  B12 <- 2 * a * (a + 1) * b ^ 3
  B <- rbind(cbind(diag(B11), diag(B12)),
             cbind(diag(B12), diag(B22)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

avar_me2 <- function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  S11 <- avar_me(Dir(alpha = a))
  S21 <- - matrix(colSums(S11) * b / a0, nrow = 1)

  # Matrix D
  D <- cbind(rbind(S11, S21), c(t(S21), (sum(S11) + a0) * (b / a0) ^ 2))
  D <- as.matrix(nearPD(D))
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

}

#' @rdname Multigam
setMethod("avar_same",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A12 <- Matrix(a, k, 1) %*% Matrix(a, 1, k) / k
  A13 <- - Matrix(a, k, 1) %*% Matrix(1, 1, k) / (k * b)
  A21 <- - (digamma(a) + log(b)) / k
  A22 <- - a * b / k
  A23 <- rep(1 / k, k)
  A11 <- (- Matrix(a, k, 1) %*% Matrix(A21, 1, k) + diag(1, k, k)) / b
  A <- rbind(cbind(A11, A12, A13), Matrix(c(A21, A22, A23), 1, 3 * k))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- trigamma(a)
  B33 <- a * (a + 1) * b ^ 2 *
    (trigamma(a + 2) + (digamma(a + 2) + log(b)) ^ 2) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b)) ^ 2
  B12 <- rep(b, k)
  B13 <- a * (a + 1) * b ^ 2 * (digamma(a + 2) + log(b)) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b))
  B23 <- a * b * (trigamma(a + 1) + (digamma(a + 1) + log(b)) ^ 2) -
    a * b * (digamma(a) + log(b)) * (digamma(a + 1) + log(b))
  B <- rbind(cbind(diag(B11), diag(B12), diag(B13)),
             cbind(diag(B12), diag(B22), diag(B23)),
             cbind(diag(B13), diag(B23), diag(B33)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

avar_same2 <- function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  S11 <- avar_same(Dir(alpha = a))
  S21 <- - matrix(colSums(S11) * b / a0, nrow = 1)

  # Matrix D
  D <- cbind(rbind(S11, S21), c(t(S21), (sum(S11) + a0) * (b / a0) ^ 2))
  D <- as.matrix(nearPD(D))
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

}
