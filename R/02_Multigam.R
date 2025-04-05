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

#' @title Multivariate Gamma Distribution
#' @name Multigam
#'
#' @description
#' The multivariate gamma distribution is a multivariate absolute continuous
#' probability distribution, defined as the cumulative sum of independent
#' gamma random variables with possibly different shape parameters
#' \eqn{\alpha_i > 0, i\in\[k \]} and the same scale \eqn{\beta > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Gamma` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Gamma` instead.
#' @param shape,scale numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#' @param log logical. Should the logarithm of the density be returned?
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization.
#'
#' @details
#' The probability density function (PDF) of the multivariate gamma distribution
#' is given by:
#' \deqn{ f(x; \alpha, \beta) = \frac{\beta^{-\alpha_0}}{\prod_{i=1}^k\Gamma(\alpha_i)}, e^{-x_k/\beta} x_1^{\alpha_1-1}\prod_{i=1}^k (x_i - x_{i-1})^{(\alpha_i-1)} \quad x > 0. }
#'
#' @inherit Distributions return
#'
#' @references
#'
#' - Oikonomidis, I. & Trevezas, S. (2025), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Multivariate Gamma Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- c(0.5, 3, 5) ; b <- 5
#' D <- Multigam(a, b)
#' x <- c(0.3, 2, 10)
#' n <- 100
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, x) # density function
#' x <- r(D, n) # random generator function
#'
#' # alternative way to use the function
#' df <- d(D) ; df(x) # df is a function itself
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
#' var(D) # Variance
#' finf(D) # Fisher Information Matrix
#'
#' # List of all available moments
#' mom <- moments(D)
#' mom$mean # expectation
#'
#' # ------------------
#' # Point Estimation
#' # ------------------
#'
#' ll(D, x)
#' llmultigam(x, a, b)
#'
#' emultigam(x, type = "mle")
#' emultigam(x, type = "me")
#' emultigam(x, type = "same")
#'
#' mle(D, x)
#' me(D, x)
#' same(D, x)
#' e(D, x, type = "mle")
#'
#' mle("multigam", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vmultigam(a, b, type = "mle")
#' vmultigam(a, b, type = "me")
#' vmultigam(a, b, type = "same")
#'
#' avar_mle(D)
#' avar_me(D)
#' avar_same(D)
#'
#' avar(D, type = "mle")
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

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != length(shape)) {
    stop("the columns of x must be equal to the length of shape")
  }
  if (any(shape <= 0) || scale <= 0) {
    stop("shape and scale must be positive")
  }

  y <- apply(x, 1,
             FUN = function(x) {
               if (any(x <= 0)) {
                 -Inf
               } else {
                 sum(shape * log(fd(x))) - sum(shape) * log(scale) -
                   sum(lgamma(shape)) - x[length(x)] / scale
               }
             })

  if (!log) {
    return(exp(y))
  } else {
    return(y)
  }

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
setMethod("d", signature = c(distr = "Multigam", x = "numeric"),
          function(distr, x) {
            dmultigam(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Multigam
setMethod("d", signature = c(distr = "Multigam", x = "matrix"),
          function(distr, x) {
            dmultigam(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Multigam
setMethod("r", signature = c(distr = "Multigam", n = "numeric"),
          function(distr, n) {
            rmultigam(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

  cumsum(x@shape * x@scale ^ 2)

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

#' @rdname Multigam
#' @export
llmultigam <- function(x, shape, scale) {
  ll(Multigam(shape, scale), x)
}

#' @rdname Multigam
setMethod("ll",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x) {

  k <- length(distr@shape)
  sum(apply(x, MARGIN = 1, FUN = dmultigam,
            shape = distr@shape, scale = distr@scale, log = TRUE))

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
  a <- idigamma(logz - log(b))

  db <- - xk / par ^ 2
  da <- 1 / (par * trigamma(a))

  - sum(da) * log(b) - sum(a) * db / b - sum(digamma(a) * da) +
    xk * db / b ^ 2 + sum(logz * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
#' @export
emultigam <- function(x, type = "mle", ...) {

  e(Multigam(), x, type, ...)

}

#' @rdname Multigam
setMethod("mle",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  k <- ncol(x)
  logz <- colMeans(log(fd(x)))
  xk <- mean(x[, k])
  tx <- c(logz, xk)

  par <- optim(par = sum(do.call(par0, list(distr = distr, x = x))$shape),
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

  list(shape = a, scale = b)

})

#' @rdname Multigam
setMethod("me",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x) {

  z <- fd(x)
  mz <- colMeans(z)
  scale <- mean(colVar(z) / mz)
  shape <- mz / scale

  list(shape = shape, scale = scale)

})

me2 <- function(distr, x) {

  w <- gendir(x)
  shape <- unname(me(Dir(), w))
  xk <- mean(x[, ncol(x)])
  scale <- xk / sum(shape)

  list(shape = shape, scale = scale)

}

#' @rdname Multigam
setMethod("same",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x) {

  z <- fd(x)
  scale <- mean(diag(stats::cov(z, log(z))))
  shape <- colMeans(z) / scale

  list(shape = shape, scale = scale)

})

same2 <- function(distr, x) {

  w <- gendir(x)
  shape <- unname(same(Dir(), w))
  xk <- mean(x[, ncol(x)])
  scale <- xk / sum(shape)

  list(shape = shape, scale = scale)

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
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
