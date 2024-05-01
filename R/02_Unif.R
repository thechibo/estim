# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unif Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Unif",
         contains = "Distribution",
         slots = c(min = "numeric", max = "numeric"),
         prototype = list(min = 0, max = 1))

#' @title Uniform Distribution
#' @name Unif
#'
#' @param x an object of class `Unif`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Unif`.
#' @param min,max numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Unif <- function(min = 0, max = 1) {
  new("Unif", min = min, max = max)
}

setValidity("Unif", function(object) {
  if(length(object@min) != 1) {
    stop("min has to be a numeric of length 1")
  }
  if(length(object@max) != 1) {
    stop("max has to be a numeric of length 1")
  }
  if(object@min < object@max) {
    stop("min must be less than max")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
setMethod("d", signature = c(x = "Unif"),
          function(x) {
            function(y, log = FALSE) {
              dunif(y, min = x@min, max = x@max, log = log)
            }
          })

#' @rdname Unif
setMethod("p", signature = c(x = "Unif"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              punif(q, min = x@min, max = x@max,
                     lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Unif
setMethod("q2", signature = c(x = "Unif"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qunif(p, min = x@min, max = x@max,
                     lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Unif
setMethod("r", signature = c(x = "Unif"),
          function(x) {
            function(n) {
              runif(n, min = x@min, max = x@max)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
setMethod("mean",
          signature  = c(x = "Unif"),
          definition = function(x) {

  (x@max + x@min) / 2

})

#' @rdname Unif
setMethod("var",
          signature  = c(x = "Unif"),
          definition = function(x) {

  (x@max - x@min) ^ 2 / 12

})

#' @rdname Unif
setMethod("sd",
          signature  = c(x = "Unif"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Unif
setMethod("skew",
          signature  = c(x = "Unif"),
          definition = function(x) {

  0

})

#' @rdname Unif
setMethod("kurt",
          signature  = c(x = "Unif"),
          definition = function(x) {

  - 1.2

})

#' @rdname Unif
setMethod("entro",
          signature  = c(x = "Unif"),
          definition = function(x) {

  log(x@max - x@min)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llunif <- function(x, min, max) {
  ll(x, prm = c(min, max), distr = Unif())
}

#' @rdname Unif
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Unif"),
          definition = function(x, prm, distr) {

  if (max(x) > prm[2] || min(x) < prm[1]) {
    return(0)
  } else {
    return(- length(x) * log(prm[2] - prm[1]))
  }

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
eunif <- function(x, type = "mle", ...) {

  estim(x, Unif(), type, ...)

}

#' @rdname Unif
setMethod("mle",
          signature  = c(x = "numeric", distr = "Unif"),
          definition = function(x, distr) {

  c(min = min(x), max = max(x))

})

#' @rdname Unif
setMethod("me",
          signature  = c(x = "numeric", distr = "Unif"),
          definition = function(x, distr) {

  m <- mean(x)
  s <- sqrt(3) * bsd(x)
  c(min = m - s, max = m + s)

})
