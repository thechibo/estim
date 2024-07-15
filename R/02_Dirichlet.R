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
#' @param distr an object of class `Dir`.
#' @param alpha numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
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
setMethod("d", signature = c(x = "Dir"),
          function(x) {
            function(y, log = FALSE) {
              extraDistr::ddirichlet(y, alpha = x@alpha, log = log)
            }
          })

#' @rdname Dir
setMethod("r", signature = c(x = "Dir"),
          function(x) {
            function(n) {
              extraDistr::rdirichlet(n, alpha = x@alpha)
            }
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
  nearPD(y)

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

#' @rdname ll
#' @export
lldirichlet <- function(x, alpha) {
  ll(x, prm = alpha, distr = Dir())
}

#' @rdname Dir
setMethod("ll",
          signature  = c(x = "matrix", prm = "numeric", distr = "Dir"),
          definition = function(x, prm, distr) {

  nrow(x) * (lgamma(sum(prm)) - sum(lgamma(prm))) + sum(log(x) %*% diag(prm - 1))

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

#' @rdname estimation
#' @export
edirichlet <- function(x, type = "mle", ...) {

  estim(x, Dir(), type, ...)

}

#' @rdname Dir
setMethod("mle",
          signature  = c(x = "matrix", distr = "Dir"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- colMeans(log(x))

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  alpha <- idigamma(digamma(par) + tx)

  names(alpha) <- paste0("alpha", seq_along(alpha))
  alpha

})

#' @rdname Dir
setMethod("me",
          signature  = c(x = "matrix", distr = "Dir"),
          definition = function(x, distr) {

  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  alpha  <- m * (m - m2) / (m2 - m ^ 2)

  names(alpha) <- paste0("alpha", seq_along(alpha))
  alpha

})

#' @rdname Dir
setMethod("same",
          signature  = c(x = "matrix", distr = "Dir"),
          definition = function(x, distr) {

  m  <- colMeans(x)
  logm  <- colMeans(log(x))
  mlogm <- colMeans(x * log(x))

  alpha  <- (length(m) - 1) * m / sum(mlogm - m * logm)

  names(alpha) <- paste0("alpha", seq_along(alpha))
  alpha

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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

  D <- solve(diag(trigamma(a)) - matrix(trigamma(sum(a)), k, k))
  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

#' @rdname Dir
setMethod("avar_me",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  a <- distr@alpha
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)

  A1 <- diag(a0 * (2 * a0 + 1) * (a + 1) / b)
  A2 <- diag(- a0 * (a0 + 1) ^ 2 / b)
  A <- Matrix(cbind(A1, A2))

  B11 <- - Matrix(a, k, 1) %*% Matrix(a, 1, k)
  diag(B11) <- a * b
  B11 <- B11 / (a0 ^ 2 * (a0 + 1))
  B11 <- nearPD(B11)

  B12 <- - Matrix(a, k, 1) %*% Matrix(a * (a + 1), 1, k)
  diag(B12) <- a * (a + 1) * b
  B12 <- B12 * 2 / (a0 ^ 2 * (a0 + 1) * (a0 + 2))

  c22 <- - 2 * (2 * a0 + 3) / (a0 ^ 2 * (a0 + 1) ^ 2 * (a0 + 2) * (a0 + 3))
  B22 <- c22 * Matrix(a * (a + 1), k, 1) %*% Matrix(a * (a + 1), 1, k)
  diag(B22) <- (a * (a + 1) * (a + 2) * (a + 3)) /
    (a0 * (a0 + 1) * (a0 + 2) * (a0 + 3)) -
    (a * (a + 1) / (a0 * (a0 + 1))) ^ 2
  B22 <- nearPD(B22)

  B <- rbind(cbind(B11, B12),
             cbind(Matrix::t(B12), B22))
  B <- nearPD(B)

  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
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

  # Matrix A

  A1 <- (a0 / (k - 1)) * Matrix(a, k, 1) %*% Matrix(Ddigamma(a, a0), 1, k) +
    a0 * diag(k)
  A2 <- (1 / (k - 1)) * Matrix(a, k, 1) %*% Matrix(a, 1, k)
  A3 <- - (a0 / (k - 1)) * Matrix(a, k, 1) %*% Matrix(1, 1, k)
  A <- cbind(A1, A2, A3)

  # Matrix B

  B11 <- - Amat
  diag(B11) <- a * b
  B11 <- B11 / (a0 ^ 2 * (a0 + 1))
  B11 <- nearPD(B11)

  B22 <- trigamma(a) * Ik - Matrix(trigamma(a0), k, k)
  B22 <- nearPD(B22)

  c331 <- Ddigamma(a + 1, a0 + 2)
  B331 <- Matrix(c331, k, 1) %*% Matrix(c331, 1, k)
  c332 <- Ddigamma(a + 1, a0 + 1)
  B332 <- Matrix(c332, k, 1) %*% Matrix(c332, 1, k)

  B33 <- Amat * (B331 - trigamma(a0 + 2)) / (a0 * (a0 + 1)) -
    Amat * B332 / (a0 ^ 2)
  diag(B33) <- (Ddigamma(a + 2, a0 + 2) ^ 2 + Dtrigamma(a + 2, a0 + 2))  *
    a * (a + 1) / (a0 * (a0 + 1)) - (Ddigamma(a + 1, a0 + 1) * a / a0) ^ 2
  B33 <- nearPD(B33)

  B12 <- Ik / a0 - Matrix(a, k, 1) %*% Matrix(1, 1, k) / a0 ^ 2

  c13 <- a * (Ddigamma(a + 1, a0 + 2) + 1) / (a0 ^ 2 * (a0 + 1))
  B13 <- - Matrix(a, k, 1) %*% Matrix(c13, 1, k)
  diag(B13) <- - Matrix::diag(B13) * b / a

  c231 <- (Ddigamma(a + 1, a0 + 1) + a * trigamma(a + 1)) / a0
  c232 <- (Ddigamma(a + 1, a0 + 1) / a0 + trigamma(a0 + 1)) * a / a0
  B23 <- c231 * Ik - Matrix(1, k, 1) %*% Matrix(c232, 1, k)

  B <- rbind(cbind(B11, B12, B13),
             cbind(Matrix::t(B12), B22, B23),
             cbind(Matrix::t(B13), Matrix::t(B23), B33))
  B <- nearPD(B)

  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})
