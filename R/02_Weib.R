# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Weib Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Weib",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Weibull Distribution
#' @name Weib
#'
#' @param x an object of class `Weib`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Weib`.
#' @param shape,scale numeric. The distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Weib <- function(shape = 1, scale = 1) {
  new("Weib", shape = shape, scale = scale)
}

setValidity("Weib", function(object) {
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

#' @rdname Weib
setMethod("d", signature = c(distr = "Weib", x = "numeric"),
          function(distr, x) {
            dweibull(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Weib
setMethod("p", signature = c(distr = "Weib", x = "numeric"),
          function(distr, x) {
            pweibull(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Weib
setMethod("qn", signature = c(distr = "Weib", x = "numeric"),
          function(distr, x) {
            qweibull(x, shape = distr@shape, scale = distr@scale)
          })

#' @rdname Weib
setMethod("r", signature = c(distr = "Weib", n = "numeric"),
          function(distr, n) {
            rweibull(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
setMethod("mean",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * gamma(1 + 1 / x@shape)

})

#' @rdname Weib
setMethod("median",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * log(2) ^ (1 / x@shape)

})

#' @rdname Weib
setMethod("mode",
          signature  = c(x = "Weib"),
          definition = function(x) {

  if (x@shape > 1) {
    return(x@scale * (1 - 1 / x@shape) ^ (1 / x@shape))
  } else {
    return(0)
  }

})

#' @rdname Weib
setMethod("var",
          signature  = c(x = "Weib"),
          definition = function(x) {

  (gamma(1 + 2 / x@shape) - gamma(1 + 1 / x@shape) ^ 2) * x@scale ^ 2

})

#' @rdname Weib
setMethod("sd",
          signature  = c(x = "Weib"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Weib
setMethod("skew",
          signature  = c(x = "Weib"),
          definition = function(x) {

  m <- mean(x)
  s2 <- var(x)

  (gamma(1 + 3 / x@shape) * x@scale ^ 3 - 3 * m * s2 - m ^ 3) / s2 ^ 1.5

})

#' @rdname Weib
setMethod("kurt",
          signature  = c(x = "Weib"),
          definition = function(x) {

  g1 <- gamma(1 + 1 / x@shape)
  g2 <- gamma(1 + 2 / x@shape)
  g3 <- gamma(1 + 3 / x@shape)
  g4 <- gamma(1 + 4 / x@shape)

  (- 6 * g1 ^ 4 + 12 * g1 ^ 2 * g2 - 3 * g2 ^ 2 - 4 * g1 * g3 + g4) /
    (g2 - g1 ^ 2) ^ 2

})

#' @rdname Weib
setMethod("entro",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * (1 - 1 / x@shape) + log(x@scale / x@shape) + 1

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
#' @export
llweib <- function(x, shape, scale) {
  ll(distr = Weib(shape, scale), x)
}

#' @rdname Weib
setMethod("ll",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x) {


})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Weib"),
          definition = function(par, tx, distr) {


          })

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Weib"),
          definition = function(par, tx, distr) {


          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
#' @export
eweib <- function(x, type = "mle", ...) {

  e(Weib(), x, type, ...)

}

#' @rdname Weib
setMethod("mle",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

            tx <- c(log(mean(x)), mean(log(x)))

            par <- optim(par = do.call(par0, list(distr = distr, x = x))[1],
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

#' @rdname Weib
setMethod("me",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x) {


          })

#' @rdname Weib
setMethod("same",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x) {


          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
#' @export
vweib <- function(shape, scale, type = "mle") {

  avar(Weib(shape = shape, scale = scale), type = type)

}

#' @rdname Weib
setMethod("avar_mle",
          signature  = c(distr = "Weib"),
          definition = function(distr) {


          })

#' @rdname Weib
setMethod("avar_me",
          signature  = c(distr = "Weib"),
          definition = function(distr) {


          })

#' @rdname Weib
setMethod("avar_same",
          signature  = c(distr = "Weib"),
          definition = function(distr) {


          })
