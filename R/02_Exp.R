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
#' @param distr an object of class `Exp`.
#' @param rate numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Exp"),
          function(x) {
            function(y, log = FALSE) {
              dexp(y, rate = x@rate, log = log)
            }
          })

#' @rdname Exp
setMethod("p", signature = c(x = "Exp"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pexp(q, rate = x@rate, lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Exp
setMethod("qn", signature = c(x = "Exp"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qexp(p, rate = x@rate, lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Exp
setMethod("r", signature = c(x = "Exp"),
          function(x) {
            function(n) {
              rexp(n, rate = x@rate)
            }
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

#' @rdname ll
#' @export
llexp <- function(x, rate) {
  ll(x, prm = rate, distr = Exp())
}

#' @rdname Exp
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Exp"),
          definition = function(x, prm, distr) {

  length(x) * log(prm) - prm * sum(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
eexp <- function(x, type = "mle", ...) {

  estim(x, Exp(), type, ...)

}

#' @rdname Exp
setMethod("mle",
          signature  = c(x = "numeric", distr = "Exp"),
          definition = function(x, distr) {

  c(rate = 1 / mean(x))

})

#' @rdname Exp
setMethod("me",
          signature  = c(x = "numeric", distr = "Exp"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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
