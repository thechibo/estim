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
#' @param distr an object of class `Gam`.
#' @param shape,scale numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
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
setMethod("d", signature = c(x = "Gam"),
          function(x) {
            function(y, log = FALSE) {
              dgamma(y, shape = x@shape, scale = x@scale, log = log)
            }
          })

#' @rdname Gam
setMethod("p", signature = c(x = "Gam"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pgamma(q, shape = x@shape, scale = x@scale,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Gam
setMethod("qn", signature = c(x = "Gam"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qgamma(p, shape = x@shape, scale = x@scale,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Gam
setMethod("r", signature = c(x = "Gam"),
          function(x) {
            function(n) {
              rgamma(n, shape = x@shape, scale = x@scale)
            }
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

#' @rdname ll
#' @export
llgamma <- function(x, shape, scale) {
  ll(x, prm = c(shape, scale), distr = Gam())
}

#' @rdname Gam
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Gam"),
          definition = function(x, prm, distr) {

  - length(x) * (lgamma(prm[1]) + prm[1] * log(prm[2])) +
  (prm[1] - 1) * sum(log(x)) - sum(x) / prm[2]

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

#' @rdname estimation
#' @export
egamma <- function(x, type = "mle", ...) {

  estim(x, Gam(), type, ...)

}

#' @rdname Gam
setMethod("mle",
          signature  = c(x = "numeric", distr = "Gam"),
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

#' @rdname Gam
setMethod("me",
          signature  = c(x = "numeric", distr = "Gam"),
          definition = function(x, distr) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  s2 <- m2 - m ^ 2
  c(shape = m ^ 2 / s2, scale = s2 / m)

})

#' @rdname Gam
setMethod("same",
          signature  = c(x = "numeric", distr = "Gam"),
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
