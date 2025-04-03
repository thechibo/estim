# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chisq Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Chisq",
         contains = "Distribution",
         slots = c(df = "numeric"),
         prototype = list(df = 1))

#' @title Chi-Square Distribution
#' @name Chisq
#'
#' @param x an object of class `Chisq`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Chisq`.
#' @param df numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Chisq <- function(df = 1) {
  new("Chisq", df = df)
}

setValidity("Chisq", function(object) {
  if(length(object@df) != 1) {
    stop("df has to be a numeric of length 1")
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
setMethod("d", signature = c(distr = "Chisq", x = "numeric"),
          function(distr, x) {
            dchisq(x, df = distr@df, ncp = 0)
          })

#' @rdname Chisq
setMethod("p", signature = c(distr = "Chisq", x = "numeric"),
          function(distr, x) {
            pchisq(x, df = distr@df, ncp = 0)
          })

#' @rdname Chisq
setMethod("qn", signature = c(distr = "Chisq", x = "numeric"),
          function(distr, x) {
            qchisq(x, df = distr@df, ncp = 0)
          })

#' @rdname Chisq
setMethod("r", signature = c(distr = "Chisq", n = "numeric"),
          function(distr, n) {
            rchisq(n, df = distr@df, ncp = 0)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
setMethod("mean",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  x@df

})

#' @rdname Chisq
setMethod("median",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  warning("The median of the Chi-Squared Distribution is not
          available in closed-form. An approximation is provided.")

  k <- x@df
  k * (1 - 2 / (9 * k)) ^ 3

})

#' @rdname Chisq
setMethod("mode",
          signature  = c(x = "Chisq"),
          definition = function(x) {

            max(x@df - 2, 0)

          })

#' @rdname Chisq
setMethod("var",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  2 * x@df

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

 2 ^ 1.5 / sqrt(x@df)

})

#' @rdname Chisq
setMethod("kurt",
          signature  = c(x = "Chisq"),
          definition = function(x) {

12 / x@df

})

#' @rdname Chisq
setMethod("entro",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  k2 <- x@df / 2

  k2 + log(2 * gamma(k2)) + (1 - k2) * digamma(k2)

})

#' @rdname Chisq
setMethod("finf",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  c(df = 0.25 * trigamma(x@df / 2))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
llchisq <- function(x, df) {
  ll(Chisq(df), x)
}

#' @rdname Chisq
setMethod("ll",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x) {

    df <- distr@df
  - length(x) * (lgamma(df / 2) + df * log(2) / 2) +
    (df / 2 - 1) * sum(log(x)) - sum(x) / 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
echisq <- function(x, type = "mle", ...) {

  e(Chisq(), x, type, ...)

}

#' @rdname Chisq
setMethod("mle",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x) {

  list(df = 2 * idigamma(mean(log(x)) - log(2)))

})

#' @rdname Chisq
setMethod("me",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x) {

  list(df = mean(x))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
vchisq <- function(df, type = "mle") {

  avar(Chisq(df), type = type)

}

#' @rdname Chisq
setMethod("avar_mle",
          signature  = c(distr = "Chisq"),
          definition = function(distr) {

  1 / finf(distr)

})

#' @rdname Chisq
setMethod("avar_me",
          signature  = c(distr = "Chisq"),
          definition = function(distr) {

  c(df = 2 * distr@df)

})
