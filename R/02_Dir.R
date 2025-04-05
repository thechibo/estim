# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dir Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Dir",
         contains = "Distribution",
         slots = c(alpha = "numeric"),
         prototype = list(alpha = c(1, 1)))

#' @title Dirichlet Distribution
#' @name Dir
#'
#' @description
#' The Dirichlet distribution is an absolute continuous probability,
#' specifically a multivariate generalization of the beta distribution,
#' parameterized by a vector \eqn{\boldsymbol{\alpha} =
#' (\alpha_1, \alpha_2, ..., \alpha_k)} with \eqn{\alpha_i > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Dir` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Dir` instead.
#' @param alpha numeric. The distribution parameter vector.
#' @param log logical. Should the logarithm of the density be returned?
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization.
#'
#' @details
#' The probability density function (PDF) of the Dirichlet distribution is given
#' by:
#' \deqn{ f(x_1, ..., x_k; \alpha_1, ..., \alpha_k) = \frac{1}{B(\boldsymbol{\alpha})} \prod_{i=1}^k x_i^{\alpha_i - 1}, }
#' where \eqn{B(\boldsymbol{\alpha})} is the multivariate Beta function:
#' \deqn{ B(\boldsymbol{\alpha}) = \frac{\prod_{i=1}^k \Gamma(\alpha_i)}{\Gamma\left(\sum_{i=1}^k \alpha_i\right)} }
#' and \eqn{\sum_{i=1}^k x_i = 1}, \eqn{x_i > 0}.
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
#' # Dir Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- c(0.5, 2, 5)
#' D <- Dir(a)
#' x <- c(0.3, 0.2, 0.5)
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
#' mode(D) # Mode
#' var(D) # Variance
#' entro(D) # Entropy
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
#' lldirichlet(x, a)
#'
#' edirichlet(x, type = "mle")
#' edirichlet(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("dir", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vdirichlet(a, type = "mle")
#' vdirichlet(a, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Dir <- function(alpha = c(1, 1)) {
  new("Dir", alpha = alpha)
}

setValidity("Dir", function(object) {
  if(length(object@alpha) <= 1) {
    stop("alpha has to be a numeric of length 2 or more")
  }
  if(any(object@alpha <= 0)) {
    stop("alpha has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
ddir <- function(x, alpha, log = FALSE) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != length(alpha)) {
    stop("the columns of x must be equal to the length of alpha")
  }
  if (any(alpha <= 0)) {
    stop("alpha must be positive")
  }

  y <- apply(x, 1,
        FUN = function(x) {
          if (any(x <= 0) || abs(sum(x) - 1) > 1e-10) {
            -Inf
          } else {
            lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
          }
        })

  if (!log) {
    return(exp(y))
  } else {
    return(y)
  }

}

#' @rdname Dir
#' @export
rdir <- function(n, alpha) {

  if (length(n) > 1 || !is.numeric(n) || n < 0) {
    stop("n must be a positive integer")
  }
  if (any(alpha <= 0)) {
    stop("alpha must be positive")
  }

  k <- length(alpha)
  m <- matrix(rgamma(n * k, shape = alpha), nrow = n, ncol = k, byrow = TRUE)

  m / rowSums(m)

}

#' @rdname Dir
setMethod("d", signature = c(distr = "Dir", x = "numeric"),
          function(distr, x) {
            ddir(x, alpha = distr@alpha)
          })

#' @rdname Dir
setMethod("d", signature = c(distr = "Dir", x = "matrix"),
          function(distr, x) {
            ddir(x, alpha = distr@alpha)
          })

#' @rdname Dir
setMethod("r", signature = c(distr = "Dir", n = "numeric"),
          function(distr, n) {
            rdir(n, alpha = distr@alpha)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
setMethod("mean",
          signature  = c(x = "Dir"),
          definition = function(x) {

  x@alpha / sum(x@alpha)

})

#' @rdname Dir
setMethod("mode",
          signature  = c(x = "Dir"),
          definition = function(x) {

  (x@alpha - 1) / (sum(x@alpha) - length(x@alpha))

})

#' @rdname Dir
setMethod("var",
          signature  = c(x = "Dir"),
          definition = function(x) {

  # Required variables
  a <- x@alpha
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)
  Ik <- diag(k)

  y <- - Matrix(a, k, 1) %*% Matrix(a, 1, k)
  diag(y) <- a * b
  y <- y / (a0 ^ 2 * (a0 + 1))
  as.matrix(nearPD(y))

})

#' @rdname Dir
setMethod("entro",
          signature  = c(x = "Dir"),
          definition = function(x) {

  a <- x@alpha
  a0 <- sum(a)
  ba <- sum(lgamma(a)) - lgamma(a0)

  ba + (a0 - length(a)) * digamma(a0) - sum((a - 1) * digamma(a))

})

#' @rdname Dir
setMethod("finf",
          signature  = c(x = "Dir"),
          definition = function(x) {

  a <- x@alpha
  k <- length(a)

  D <- diag(trigamma(a)) - matrix(trigamma(sum(a)), k, k)

  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
lldirichlet <- function(x, alpha) {
  ll(distr = Dir(alpha), x)
}

#' @rdname Dir
setMethod("ll",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x) {

  a <- distr@alpha
  nrow(x) * (lgamma(sum(a)) - sum(lgamma(a))) + sum(log(x) %*% diag(a - 1))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dir"),
          definition = function(par, tx, distr) {

  a <- idigamma(digamma(par) + tx)

  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * tx)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dir"),
          definition = function(par, tx, distr) {

  # Shape parameters (a_i) as a function of a0
  a <- idigamma(digamma(par) + tx)

  # a_i derivative wrt a0
  da <- trigamma(par) / trigamma(a)

  # lloptim derivative wrt a0 (par)
  digamma(sum(a)) * sum(da) - sum(digamma(a) * da) + sum(tx * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
edirichlet <- function(x, type = "mle", ...) {

  e(Dir(), x, type, ...)

}

#' @rdname Dir
setMethod("mle",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- colMeans(log(x))

  par <- optim(par = sum(unlist(do.call(par0, list(distr = distr, x = x)))),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  list(alpha = idigamma(digamma(par) + tx))

})

#' @rdname Dir
setMethod("me",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x) {

  m  <- colMeans(x)
  m2  <- colMeans(x ^ 2)

  a0 <- (1 - sum(m2)) / (sum(m2) - sum(m ^ 2))

  list(alpha = a0 * m)

})

#' @rdname Dir
setMethod("same",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x) {

  m  <- colMeans(x)
  logm  <- colMeans(log(x))
  mlogm <- colMeans(x * log(x))

  list(alpha  = (length(m) - 1) * m / sum(mlogm - m * logm))

})

me1dir <- function(x) {

  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)

  list(alpha = m * (m - m2) / (m2 - m ^ 2))

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
vdirichlet <- function(alpha, type = "mle") {

  avar(Dir(alpha = alpha), type = type)

}

#' @rdname Dir
setMethod("avar_mle",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  a <- distr@alpha
  k <- length(a)
  a0 <- sum(a)
  trig <- 1 / trigamma(a)
  cons <- trigamma(a0) / (1 - trigamma(a0) * sum(trig))

  D <- diag(trig) + cons * Matrix(trig, k, 1) %*% Matrix(trig, 1, k)

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

#' @rdname Dir
setMethod("avar_me",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@alpha
  a0 <- sum(a)
  k <- length(a)
  dn <- a0 ^ 2 - sum(a ^ 2)
  a1 <- a0 + 1
  a2 <- a0 + 2
  a3 <- a0 + 3
  s2 <- sum(a^2)
  s3 <- sum(a^3)

  c0 <- ((- 4 * a0 * (a0 - 1) * a1 ^ 2 * s3 +
         (2 * a0 ^ 3 + a0 ^ 2 + a0) * s2 ^ 2 +
         (2 * a0 ^ 5 + 2 * a0 ^ 4 - 6 * a0 ^ 3 - 4 * a0 ^ 2 - 2 * a0) * s2 +
         a0 ^ 6 + a0 ^ 5 + 2 * a0 ^ 3) / (dn ^ 2 * a1 * a2 * a3))

  D <- (a0 / a1) * diag(a) + 2*a0 / (dn*a2) * (a %*% t(a^2) + a^2 %*% t(a)) +
    c0 * a %*% t(a)

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

#' @rdname Dir
setMethod("avar_same",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  # Required variables
  a <- distr@alpha
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)
  Ik <- diag(k)
  Amat <- Matrix(a, k, 1) %*% Matrix(a, 1, k)
  mat2 <- Matrix(1/a, k, 1) %*% Matrix(1, 1, k)

  par1 <- 1 / ((k - 1) * (a0 + 1))
  par2 <- 1 / ((k - 1) ^ 2 * (a0 + 1))

  c <- - par2 *  sum(a^2*trigamma(a)) + a0 * par2 * sum(a*trigamma(a)) +
    (a0+2) * par1

  D <- Amat * (c - (mat2 + Matrix::t(mat2)) * a0 * par1  +
                 a0 * diag(1/ a) / (a0+1))

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

avar_me1dir <- function(distr) {

  a <- distr@alpha
  k <- length(a)
  a0 <- sum(a)
  b <- a0 - a

  matai <- Matrix(a, k, 1) %*% Matrix(1, 1, k)
  mataj <- Matrix(1, k, 1) %*% Matrix(a, 1, k)

  com <- (Matrix((a + 1) / b, k, 1) %*% Matrix((a + 1) / b, 1, k)) *
    (diag(a0, k, k) - matai) * mataj * a0 / (a0 + 2)

  A <- diag(1 / (a + 1)) * 2 * (a0 + 1) ^ 2 / (a0 + 3)
  B <- (2 * a0 ^ 2 + a0 + 1) / ((a0 + 1) * (a0 + 3))

  D <- com * (A - B)

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

}
