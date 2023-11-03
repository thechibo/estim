# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "numeric", distr = "Gammad"),
          definition = function(prm, x, distr) {

  sum(dgamma(x = x, shape = prm[1], scale = prm[2], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MLE                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Gammad"),
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
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  dn <- list(prm = c("shape", "scale"), sam = seqcol(x))
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
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::param(distr)@scale

  psi1a <- trigamma(a)
  d <- 1 / (a * psi1a - 1)

  v11 <- a
  v21 <- - b
  v22 <- b ^ 2 * psi1a

  D <- matrix(d * c(v11, v21, v21, v22), nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  me(x, distr)[ , 1]

})

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr) {

  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  s2 <- m2 - m ^ 2
  a  <- m ^ 2 / s2
  b  <- s2 / m

  dn <- list(prm = c("shape", "scale"), sam = seqcol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)

  s11 <- 2 * a * (a + 1)
  s22 <- b ^ 2 * (2 * a + 3) / a
  s12 <- - 2 * b * (a + 1)
  D <- matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  same(x, distr)[ , 1]

})

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr) {

  mx  <- colMeans(x)
  mlx <- colMeans(log(x))
  mxlx <- colMeans(x * log(x))
  b <- mxlx - mx * mlx
  a  <- mx / b

  dn <- list(prm = c("shape", "scale"), sam = seqcol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)

  c1 <- 1 + a * trigamma(a + 1)
  c2 <- 1 + a * trigamma(a)

  v11 <- a ^ 2 * c1
  v21 <- - a * b * c1
  v22 <- b ^ 2 * c2

  D <- matrix(c(v11, v21, v21, v22), 2, 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})
