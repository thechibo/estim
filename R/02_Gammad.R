# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gammad Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llgamma <- function(x, shape, scale) {
  ll(x, prm = c(shape, scale), distr = distr::Gammad())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Gammad"),
          definition = function(x, prm, distr) {

  sum(dgamma(x = x, shape = prm[1], scale = prm[2], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gammad"),
          definition = function(par, tx, distr) {

  par * log(par) - lgamma(par) - (tx[1] + 1) * par + (par - 1) * tx[2]

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gammad"),
          definition = function(par, tx, distr) {

  log(par) - digamma(par) - tx[1] + tx[2]

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
egamma <- function(x, type = "mle", ...) {

  estim(x, distr::Gammad(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx <- c(log(mean(x)), mean(log(x)))

  par <- optim(par = do.call(par0, list(x = x, distr = distr))[1],
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  par <- c(par, mean(x) / par)

  names(par) <- c("shape", "scale")
  par

})

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  s2 <- m2 - m ^ 2
  c(shape = m ^ 2 / s2, scale = s2 / m)

})

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  mx  <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  cxlx <- mxlx - mx * mlx

  c(shape = mx / cxlx, scale = cxlx)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vgamma <- function(shape, scale, type = "mle") {

  avar(distr::Gammad(shape = shape, scale = scale), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
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

#' @rdname avar_me
setMethod("avar_me",
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

#' @rdname avar_same
setMethod("avar_same",
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
