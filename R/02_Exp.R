# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exp Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Exp",
         contains = "Distribution",
         slots = c(rate = "numeric"),
         prototype = list(rate = 1))

#' @title Exponential Distribution
#' @name Exp
#'
#' @param x an object of class `Exp`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Exp`.
#' @param rate numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Exp <- function(rate = 1) {
  new("Exp", rate = rate)
}

setValidity("Exp", function(object) {
  if(length(object@rate) != 1) {
    stop("rate has to be a numeric of length 1")
  }
  if(object@rate <= 0) {
    stop("rate has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
setMethod("d", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            dexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("p", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            pexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("qn", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            qexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("r", signature = c(distr = "Exp", n = "numeric"),
          function(distr, n) {
            rexp(n, rate = distr@rate)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
setMethod("mean",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate

})

#' @rdname Exp
setMethod("median",
          signature  = c(x = "Exp"),
          definition = function(x) {

  log(2) / x@rate

})

#' @rdname Exp
setMethod("mode",
          signature  = c(x = "Exp"),
          definition = function(x) {

  0

})

#' @rdname Exp
setMethod("var",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate ^ 2

})

#' @rdname Exp
setMethod("sd",
          signature  = c(x = "Exp"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Exp
setMethod("skew",
          signature  = c(x = "Exp"),
          definition = function(x) {

  2

})

#' @rdname Exp
setMethod("kurt",
          signature  = c(x = "Exp"),
          definition = function(x) {

  6

})

#' @rdname Exp
setMethod("entro",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 - log(x@rate)

})

#' @rdname Exp
setMethod("finf",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
llexp <- function(x, rate) {
  ll(Exp(rate), x)
}

#' @rdname Exp
setMethod("ll",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  rate <- distr@rate
  length(x) * log(rate) - rate * sum(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
eexp <- function(x, type = "mle", ...) {

  e(Exp(), x, type, ...)

}

#' @rdname Exp
setMethod("mle",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  list(rate = 1 / mean(x))

})

#' @rdname Exp
setMethod("me",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
vexp <- function(rate, type = "mle") {

  avar(Exp(rate = rate), type = type)

}

#' @rdname Exp
setMethod("avar_mle",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  rate <- distr@rate
  c(rate = rate ^ 2)

})

#' @rdname Exp
setMethod("avar_me",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  avar_mle(distr)

})
