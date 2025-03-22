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
#' @param distr an object of class `Chisq`.
#' @param df numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Chisq"),
          function(x) {
            function(y, log = FALSE) {
              dchisq(y, df = x@df, ncp = 0, log = log)
            }
          })

#' @rdname Chisq
setMethod("p", signature = c(x = "Chisq"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pchisq(q, df = x@df, ncp = 0,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Chisq
setMethod("qn", signature = c(x = "Chisq"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qchisq(p, df = x@df, ncp = 0,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Chisq
setMethod("r", signature = c(x = "Chisq"),
          function(x) {
            function(n) {
              rchisq(n, df = x@df, ncp = 0)
            }
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

#' @rdname ll
#' @export
llchisq <- function(x, df) {
  ll(x, prm = df, distr = Chisq())
}

#' @rdname Chisq
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Chisq"),
          definition = function(x, prm, distr) {

  - length(x) * (lgamma(prm / 2) + prm * log(2) / 2) +
    (prm / 2 - 1) * sum(log(x)) - sum(x) / 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
echisq <- function(x, type = "mle", ...) {

  estim(x, Chisq(), type, ...)

}

#' @rdname Chisq
setMethod("mle",
          signature  = c(x = "numeric", distr = "Chisq"),
          definition = function(x, distr) {

  c(df = 2 * idigamma(mean(log(x)) - log(2)))

})

#' @rdname Chisq
setMethod("me",
          signature  = c(x = "numeric", distr = "Chisq"),
          definition = function(x, distr) {

  c(df = mean(x))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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
