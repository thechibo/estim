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
#' @param distr an object of class `Pois`.
#' @param lambda numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Pois"),
          function(x) {
            function(y, log = FALSE) {
              dpois(y, lambda = x@lambda, log = log)
            }
          })

#' @rdname Pois
setMethod("p", signature = c(x = "Pois"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              ppois(q, lambda = x@lambda,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Pois
setMethod("qn", signature = c(x = "Pois"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qpois(p, lambda = x@lambda,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Pois
setMethod("r", signature = c(x = "Pois"),
          function(x) {
            function(n) {
              rpois(n, lambda = x@lambda)
            }
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
setMethod("finf",
          signature  = c(x = "Pois"),
          definition = function(x) {

  1 / x@lambda

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llpois <- function(x, lambda) {
  ll(x, prm = lambda, distr = Pois())
}

#' @rdname Pois
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Pois"),
          definition = function(x, prm, distr) {

  log(prm) * sum(x) - length(x) * prm - sum(log(factorial(x)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
epois <- function(x, type = "mle", ...) {

  estim(x, Pois(), type, ...)

}

#' @rdname Pois
setMethod("mle",
          signature  = c(x = "numeric", distr = "Pois"),
          definition = function(x, distr) {

  c(lambda = mean(x))

})

#' @rdname Pois
setMethod("me",
          signature  = c(x = "numeric", distr = "Pois"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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
