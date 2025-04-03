# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pois Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Pois",
         contains = "Distribution",
         slots = c(lambda = "numeric"),
         prototype = list(lambda = 1))

#' @title Poisson Distribution
#' @name Pois
#'
#' @param x an object of class `Pois`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Pois`.
#' @param lambda numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Pois <- function(lambda = 1) {
  new("Pois", lambda = lambda)
}

setValidity("Pois", function(object) {
  if(length(object@lambda) != 1) {
    stop("lambda has to be a numeric of length 1")
  }
  if(object@lambda <= 0) {
    stop("lambda has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Pois
setMethod("d", signature = c(distr = "Pois", x = "numeric"),
          function(distr, x) {
            dpois(x, lambda = distr@lambda)
          })

#' @rdname Pois
setMethod("p", signature = c(distr = "Pois", x = "numeric"),
          function(distr, x) {
            ppois(x, lambda = distr@lambda)
          })

#' @rdname Pois
setMethod("qn", signature = c(distr = "Pois", x = "numeric"),
          function(distr, x) {
            qpois(x, lambda = distr@lambda)
          })

#' @rdname Pois
setMethod("r", signature = c(distr = "Pois", n = "numeric"),
          function(distr, n) {
            rpois(n, lambda = distr@lambda)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Pois
setMethod("mean",
          signature  = c(x = "Pois"),
          definition = function(x) {

  x@lambda

})

#' @rdname Pois
setMethod("median",
          signature  = c(x = "Pois"),
          definition = function(x) {

            warning("The median of a Pois(l) distribution is given by the
                    inequality: l - ln2 <= median < l + 1/3. The lower bound is
                    returned.")

            x@lambda - log(2)

          })


#' @rdname Pois
setMethod("mode",
          signature  = c(x = "Pois"),
          definition = function(x) {

            floor(x@lambda)

          })

#' @rdname Pois
setMethod("var",
          signature  = c(x = "Pois"),
          definition = function(x) {

  x@lambda

})

#' @rdname Pois
setMethod("sd",
          signature  = c(x = "Pois"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Pois
setMethod("skew",
          signature  = c(x = "Pois"),
          definition = function(x) {

  1 / sqrt(x@lambda)

})

#' @rdname Pois
setMethod("kurt",
          signature  = c(x = "Pois"),
          definition = function(x) {

  1 / x@lambda

})

#' @rdname Pois
setMethod("entro",
          signature  = c(x = "Pois"),
          definition = function(x) {

  warning("The entropy given is an approximation in the O(1 / l ^ 4) order.")
  l <- x@lambda
  0.5 * log(2 * pi  * exp(1) * l) - 1 / (12 * l) - 1 / (24 * l ^ 2) -
    19 / (360 * l ^ 3)

})

#' @rdname Pois
setMethod("finf",
          signature  = c(x = "Pois"),
          definition = function(x) {

  1 / x@lambda

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Pois
#' @export
llpois <- function(x, lambda) {
  ll(Pois(lambda), x)
}

#' @rdname Pois
setMethod("ll",
          signature  = c(distr = "Pois", x = "numeric"),
          definition = function(distr, x) {

  lam <- distr@lambda
  log(lam) * sum(x) - length(x) * lam - sum(log(factorial(x)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Pois
#' @export
epois <- function(x, type = "mle", ...) {

  e(Pois(), x, type, ...)

}

#' @rdname Pois
setMethod("mle",
          signature  = c(distr = "Pois", x = "numeric"),
          definition = function(distr, x) {

  list(lambda = mean(x))

})

#' @rdname Pois
setMethod("me",
          signature  = c(distr = "Pois", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Pois
#' @export
vpois <- function(lambda, type = "mle") {

  avar(Pois(lambda = lambda), type = type)

}

#' @rdname Pois
setMethod("avar_mle",
          signature  = c(distr = "Pois"),
          definition = function(distr) {

  c(lambda = distr@lambda)

})

#' @rdname Pois
setMethod("avar_me",
          signature  = c(distr = "Pois"),
          definition = function(distr) {

  avar_mle(distr)

})
