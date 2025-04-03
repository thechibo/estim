# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gam Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Gam",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Gamma Distribution
#' @name Gam
#'
#' @param x an object of class `Gam`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Gam`.
#' @param shape,scale numeric. The distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Gam <- function(shape = 1, scale = 1) {
  new("Gam", shape = shape, scale = scale)
}

setValidity("Gam", function(object) {
  if(length(object@shape) != 1) {
    stop("shape has to be a numeric of length 1")
  }
  if(object@shape <= 0) {
    stop("shape has to be positive")
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

#' @rdname Gam
setMethod("d", signature = c(distr = "Gam", x = "numeric"),
          function(distr, x) {
            dgamma(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Gam
setMethod("p", signature = c(distr = "Gam", x = "numeric"),
          function(distr, x) {
            pgamma(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Gam
setMethod("qn", signature = c(distr = "Gam", x = "numeric"),
          function(distr, x) {
            qgamma(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Gam
setMethod("r", signature = c(distr = "Gam", n = "numeric"),
          function(distr, n) {
            rgamma(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
setMethod("mean",
          signature  = c(x = "Gam"),
          definition = function(x) {

  x@shape * x@scale

})

#' @rdname Gam
#' @export
setMethod("median",
          signature  = c(x = "Gam"),
          definition = function(x) {

  qgamma(0.5, shape = x@shape, scale = x@scale)

})

#' @rdname Gam
setMethod("mode",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  b <- x@scale
  if(a >= 1) {
    return((a - 1) / b)
  } else {
    return(0)
  }

})

#' @rdname Gam
setMethod("var",
          signature  = c(x = "Gam"),
          definition = function(x) {

  x@shape * x@scale ^ 2

})

#' @rdname Gam
setMethod("sd",
          signature  = c(x = "Gam"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Gam
setMethod("skew",
          signature  = c(x = "Gam"),
          definition = function(x) {

  2 / sqrt(x@shape)

})

#' @rdname Gam
setMethod("kurt",
          signature  = c(x = "Gam"),
          definition = function(x) {

  6 / x@shape

})

#' @rdname Gam
setMethod("entro",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  a + log(x@scale) + lgamma(a) + (1 - a) * digamma(a)

})

#' @rdname Gam
setMethod("finf",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  b <- x@scale

  D <- matrix(c(trigamma(a), 1 / b, 1 / b, a / b ^ 2),
              nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
llgamma <- function(x, shape, scale) {
  ll(Gam(shape, scale), x)
}

#' @rdname Gam
setMethod("ll",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x) {

  a <- distr@shape
  b <- distr@scale
  - length(x) * (lgamma(a) + a * log(b)) + (a - 1) * sum(log(x)) - sum(x) / b

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gam"),
          definition = function(par, tx, distr) {

  par * log(par) - lgamma(par) - (tx[1] + 1) * par + (par - 1) * tx[2]

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gam"),
          definition = function(par, tx, distr) {

  log(par) - digamma(par) - tx[1] + tx[2]

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
egamma <- function(x, type = "mle", ...) {

  e(Gam(), x, type, ...)

}

#' @rdname Gam
setMethod("mle",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx <- c(log(mean(x)), mean(log(x)))

  par <- optim(par = do.call(par0, list(distr = distr, x = x))$shape,
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  list(shape = par, scale = mean(x) / par)

})

#' @rdname Gam
setMethod("me",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  s2 <- m2 - m ^ 2

  list(shape = m ^ 2 / s2, scale = s2 / m)

})

#' @rdname Gam
setMethod("same",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x) {

  mx  <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  cxlx <- mxlx - mx * mlx

  list(shape = mx / cxlx, scale = cxlx)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
vgamma <- function(shape, scale, type = "mle") {

  avar(Gam(shape = shape, scale = scale), type = type)

}

#' @rdname Gam
setMethod("avar_mle",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Gam
setMethod("avar_me",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  a <- distr@shape
  b <- distr@scale

  s11 <- 2 * a * (a + 1)
  s22 <- b ^ 2 * (2 * a + 3) / a
  s12 <- - 2 * b * (a + 1)
  D <- matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

#' @rdname Gam
setMethod("avar_same",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  a <- distr@shape
  b <- distr@scale

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
