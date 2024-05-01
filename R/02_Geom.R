# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Geom Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Geom",
         contains = "Distribution",
         slots = c(prob = "numeric"),
         prototype = list(prob = 0.5))

#' @title Geometric Distribution
#' @name Geom
#'
#' @param x an object of class `Geom`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Geom`.
#' @param prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Geom <- function(prob = 0.5) {
  new("Geom", prob = prob)
}

setValidity("Geom", function(object) {
  if(length(object@prob) != 1) {
    stop("prob has to be a numeric of length 1")
  }
  if(object@prob <= 0 || object@prob >= 1) {
    stop("prob has to be between 0 and 1")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Geom
setMethod("d", signature = c(x = "Geom"),
          function(x) {
            function(y, log = FALSE) {
              dgeom(y, prob = x@prob, log = log)
            }
          })

#' @rdname Geom
setMethod("p", signature = c(x = "Geom"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pgeom(q, prob = x@prob,
                      lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Geom
setMethod("q2", signature = c(x = "Geom"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qgeom(p, prob = x@prob,
                      lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Geom
setMethod("r", signature = c(x = "Geom"),
          function(x) {
            function(n) {
              rgeom(n, prob = x@prob)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Geom
setMethod("mean",
          signature  = c(x = "Geom"),
          definition = function(x) {

  1 / x@prob - 1

})

#' @rdname Geom
setMethod("mode",
          signature  = c(x = "Geom"),
          definition = function(x) {

  0

})

#' @rdname Geom
setMethod("var",
          signature  = c(x = "Geom"),
          definition = function(x) {

  (1 - x@prob) / x@prob ^ 2

})

#' @rdname Geom
setMethod("sd",
          signature  = c(x = "Geom"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Geom
setMethod("skew",
          signature  = c(x = "Geom"),
          definition = function(x) {

  (2 - x@prob) / sqrt(1 - x@prob)

})

#' @rdname Geom
setMethod("kurt",
          signature  = c(x = "Geom"),
          definition = function(x) {

  6 + x@prob ^ 2 / (1 - x@prob)

})

#' @rdname Geom
setMethod("entro",
          signature  = c(x = "Geom"),
          definition = function(x) {

  p <- x@prob
  (- (1 - p) * log(1 - p) - p * log(p)) / p

})

#' @rdname Geom
setMethod("finf",
          signature  = c(x = "Geom"),
          definition = function(x) {

  1 / (x@prob ^ 2 (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llgeom <- function(x, prob) {
  ll(x, prm = prob, distr = Geom())
}

#' @rdname Geom
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Geom"),
          definition = function(x, prm, distr) {

  log(1 - prm) * sum(x) + log(prm) * length(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
egeom <- function(x, type = "mle", ...) {

  estim(x, Geom(), type, ...)

}

#' @rdname Geom
setMethod("mle",
          signature  = c(x = "numeric", distr = "Geom"),
          definition = function(x, distr) {

  c(prob = 1 / (1 + mean(x)))

})

#' @rdname Geom
setMethod("me",
          signature  = c(x = "numeric", distr = "Geom"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vgeom <- function(prob, type = "mle") {

  avar(Geom(prob = prob), type = type)

}

#' @rdname Geom
setMethod("avar_mle",
          signature  = c(distr = "Geom"),
          definition = function(distr) {

  prob <- distr@prob
  c(prob = prob ^ 2 (1 - prob))

})

#' @rdname Geom
setMethod("avar_me",
          signature  = c(distr = "Geom"),
          definition = function(distr) {

  avar_mle(distr)

})
