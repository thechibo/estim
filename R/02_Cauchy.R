# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cauchy Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cauchy",
         contains = "Distribution",
         slots = c(location = "numeric", scale = "numeric"),
         prototype = list(location = 1, scale = 1))

#' @title Cauchy Distribution
#' @name Cauchy
#'
#' @param x an object of class `Cauchy`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Cauchy`.
#' @param location,scale numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Cauchy <- function(location = 1, scale = 1) {
  new("Cauchy", location = location, scale = scale)
}

setValidity("Cauchy", function(object) {
  if(length(object@location) != 1) {
    stop("location has to be a numeric of length 1")
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

#' @rdname Cauchy
setMethod("d", signature = c(x = "Cauchy"),
          function(x) {
            function(y, log = FALSE) {
              dcauchy(y, location = x@location, scale = x@scale, log = log)
            }
          })

#' @rdname Cauchy
setMethod("p", signature = c(x = "Cauchy"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pcauchy(q, location = x@location, scale = x@scale,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Cauchy
setMethod("qn", signature = c(x = "Cauchy"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qcauchy(p, location = x@location, scale = x@scale,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Cauchy
setMethod("r", signature = c(x = "Cauchy"),
          function(x) {
            function(n) {
              rcauchy(n, location = x@location, scale = x@scale)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
setMethod("median",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

            x@location

          })

#' @rdname Cauchy
setMethod("mode",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

            x@location

          })

#' @rdname Cauchy
setMethod("entro",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

            log(4 * pi * x@scale)

          })

#' @rdname Cauchy
setMethod("finf",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

            1 / (2 * x@scale ^ 2)

          })


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llcauchy <- function(x, location, scale) {
  ll(x, prm = c(location, scale), distr = Cauchy())
}

#' @rdname Cauchy
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Cauchy"),
          definition = function(x, prm, distr) {

  - length(x) * log(pi * prm[2]) -
      sum(log (1 + ((x - prm[1]) / prm[2]) ^ 2))

})
