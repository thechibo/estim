# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chisq Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Chisq",
         contains = "Distribution",
         slots = c(df = "numeric", ncp = "numeric"),
         prototype = list(df = 1, ncp = 0))

#' @title Chi-Square Distribution
#' @name Chisq
#'
#' @param x an object of class `Chisq`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param df,ncp numeric. The distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Chisq <- function(df = 1, ncp = 0) {
  new("Chisq", df = df, ncp = ncp)
}

setValidity("Chisq", function(object) {
  if(length(object@df) != 1) {
    stop("df has to be a numeric of length 1")
  }
  if(length(object@ncp) != 1) {
    stop("ncp has to be a numeric of length 1")
  }
  if(object@df < 0) {
    stop("df has to be non-negative")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
setMethod("d", signature = c(x = "Chisq"),
          function(x) {
            function(y, log = FALSE) {
              dchisq(y, df = x@df, ncp = x@ncp, log = log)
            }
          })

#' @rdname Chisq
setMethod("p", signature = c(x = "Chisq"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pchisq(q, df = x@df, ncp = x@ncp,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Chisq
setMethod("qn", signature = c(x = "Chisq"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qchisq(p, df = x@df, ncp = x@ncp,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Chisq
setMethod("r", signature = c(x = "Chisq"),
          function(x) {
            function(n) {
              rchisq(n, df = x@df, ncp = x@ncp)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
setMethod("mean",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  x@df + x@ncp

})

#' @rdname Chisq
setMethod("var",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  2 * x@df + 4 * x@ncp

})

#' @rdname Chisq
setMethod("sd",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Chisq
setMethod("skew",
          signature  = c(x = "Chisq"),
          definition = function(x) {

 2 ^ 1.5 * (x@df + 3 * x@ncp) / (x@df + 2 * x@ncp) ^ 1.5

})

#' @rdname Chisq
setMethod("kurt",
          signature  = c(x = "Chisq"),
          definition = function(x) {

12 * (x@df + 4 * x@ncp) / (x@df + 2 * x@ncp) ^ 2

})
