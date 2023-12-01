# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llbeta <- function(x, shape1, shape2) {
  ll(x, prm = c(shape1, shape2), distr = distr::Beta())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Beta"),
          definition = function(x, prm, distr) {

  sum(dbeta(x = x, shape1 = prm[1], shape2 = prm[2], log = TRUE))

})

# Bias Corrected log-likelihood
# (Firth, 1993, Cribari-Neto and Vasconcellos, 2010)
#ll = function(prm, x) {
#  p1a = trigamma(prm[1])
#  p1b = trigamma(prm[2])
#  p1  = trigamma(sum(prm))
#  d   = p1a * p1b - p1 * (p1a + p1b)
#  ld = log((length(x) ^ 2) * d) / 2
#  sum(do.call(dbeta, c(list(x = x, log = TRUE), prm))) + ld
#}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Beta"),
          definition = function(par, tx, distr) {

  a <- idigamma(digamma(par) + tx)
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * tx)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Beta"),
          definition = function(par, tx, distr) {

  # Shape parameters (a, b) as a function of a0
  a <- idigamma(digamma(par) + tx)

  # a_i derivative wrt a0
  da <- trigamma(par) / trigamma(a)

  # lloptim derivative wrt a0 (par)
  digamma(sum(a)) * sum(da) - sum(digamma(a) * da) + sum(tx * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
ebeta <- function(x, type = "mle", ...) {

  estim(x, distr::Beta(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- c(mean(log(x)), mean(log(1 - x)))

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  shape <- idigamma(digamma(par) + tx)

  names(shape) <- paste0("shape", seq_along(shape))
  shape

})

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  d  <- (m - m2) / (m2 - m ^ 2)

  c(shape1 = d * m, shape2 = d * (1 - m))

})

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  mx <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  my <- 1 - mx
  mly <- mean(log(1 - x))
  myly <- mean((1 - x) * log(1 - x))
  s <- mxlx - mx * mlx + myly - my * mly

  c(shape1 = mx / s, shape2 = my / s)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vbeta <- function(shape1, shape2, type = "mle") {

  avar(distr::Beta(shape1 = shape1, shape2 = shape2), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  p1a  <- trigamma(a)
  p1b  <- trigamma(b)
  p1   <- trigamma(a + b)
  D    <- 1 / (p1a * p1b - (p1a + p1b) * p1)

  D <- matrix(D * c(p1b - p1, p1, p1, p1a - p1), nrow = 2, ncol = 2)

  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})

#' @rdname avar_me
setMethod("avar_me",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  s4  <- s2 ^ 2
  m3  <- 2 * (b - a) * s2 / (th * (th + 2))
  m4  <- 3 * prd * (prd * (th + 2) + 2 * (b - a) ^2) /
    ((th ^ 4) * (th + 1) * (th + 2) * (th + 3))
  d   <- (th + 1) ^ 2 * (th + 2) ^ 2 * s2
  e   <- (th + 1) ^ 3 * (m4 - s4 - m3 ^ 2 / s2) / s2

  s11 <- (a * (a + 1)) ^ 2 / d + a * e / b
  s22 <- (b * (b + 1)) ^ 2 / d + b * e / a
  s12 <- - a * (a + 1) * b * (b + 1) / d + e

  D <- matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)
  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})

#' @rdname avar_same
setMethod("avar_same",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  p1a <- trigamma(a)
  p1b <- trigamma(b)
  m1  <- matrix(c(a ^ 2, prd, prd, b ^ 2), nrow = 2, ncol = 2)
  m2  <- matrix(c(prd, th2 - prd, th2 - prd, prd), nrow = 2, ncol = 2)

  D <- (s2 * th2 * (p1a + p1b) + 1) * m1 - m2 / (th + 1)
  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})
