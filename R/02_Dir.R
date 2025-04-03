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
#' @param x an object of class `Dir`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Dir`.
#' @param alpha numeric. The distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @importFrom extraDistr ddirichlet rdirichlet
#' @export
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
setMethod("d", signature = c(distr = "Dir", x = "numeric"),
          function(distr, x) {
            extraDistr::ddirichlet(x, alpha = distr@alpha)
          })

#' @rdname Dir
setMethod("r", signature = c(distr = "Dir", n = "numeric"),
          function(distr, n) {
            extraDistr::rdirichlet(n, alpha = distr@alpha)
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
