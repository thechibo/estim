# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "numeric", distr = "Beta"),
          definition = function(prm, x, distr) {

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
## MLE                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  x  <- as.matrix(x)
  mle(x, distr)[ , 1]

})

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  dn <- list(prm = c("shape1", "shape2"), sam = seqcol(x))
  y  <- matrix(0, nrow = 2, ncol = ncol(x), dimnames = dn)

  for (j in seqcol(x)) {
    y[, j] <- optim(par = do.call(par0, list(x = x[ , j], distr = distr)),
                    fn = ll,
                    x = x[ , j],
                    distr = distr,
                    method = method,
                    lower = lower,
                    upper = upper,
                    control = list(fnscale = -1))$par
  }

  y

})

#' @rdname acov_mle
setMethod("acov_mle",
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  me(x, distr)[ , 1]

})

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  d  <- (m - m2) / (m2 - m ^ 2)
  a  <- d * m
  b  <- d * (1 - m)
  dn <- list(prm = c("shape1", "shape2"), sam = seqcol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_me
setMethod("acov_me",
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  same(x, distr)[ , 1]

})

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr) {

  x <- as.matrix(x)
  mx <- colMeans(x)
  mlx <- colMeans(log(x))
  mxlx <- colMeans(x * log(x))
  my <- 1 - mx
  mly <- colMeans(log(1 - x))
  myly <- colMeans((1 - x) * log(1 - x))

  sx <- mxlx - mx * mlx
  sy <- myly - my * mly

  a <- mx / (sx + sy)
  b <- my / (sx + sy)

  dn <- list(prm = c("shape1", "shape2"), sam = seqcol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_same
setMethod("acov_same",
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
