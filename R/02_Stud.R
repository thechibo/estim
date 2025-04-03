# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stud Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Stud",
         contains = "Distribution",
         slots = c(df = "numeric"),
         prototype = list(df = 1))

#' @title Student Distribution
#' @name Stud
#'
#' @param x an object of class `Stud`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Stud`.
#' @param df numeric. The distribution parameter.
#'
#' @inherit Distributions return
#'
#' @export
Stud <- function(df = 1) {
  new("Stud", df = df)
}

setValidity("Stud", function(object) {
  if(length(object@df) != 1) {
    stop("df has to be a numeric of length 1")
  }
  if(object@df <= 0) {
    stop("df has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Stud
setMethod("d", signature = c(distr = "Stud", x = "numeric"),
          function(distr, x) {
            dt(x, df = distr@df)
          })

#' @rdname Stud
setMethod("p", signature = c(distr = "Stud", x = "numeric"),
          function(distr, x) {
            pt(x, df = distr@df)
          })

#' @rdname Stud
setMethod("qn", signature = c(distr = "Stud", x = "numeric"),
          function(distr, x) {
            qt(x, df = distr@df)
          })

#' @rdname Stud
setMethod("r", signature = c(distr = "Stud", n = "numeric"),
          function(distr, n) {
            rt(n, df = distr@df)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Stud
setMethod("mean",
          signature  = c(x = "Stud"),
          definition = function(x) {

  df <- x@df

  if (df > 1) {
    return(0)
  } else {
    stop("Expectation is undefined for Student's t distribution
          no more than 1 df.")
  }

})

#' @rdname Stud
setMethod("median",
          signature  = c(x = "Stud"),
          definition = function(x) {

  0

})

#' @rdname Stud
setMethod("mode",
          signature  = c(x = "Stud"),
          definition = function(x) {

  0

})

#' @rdname Stud
setMethod("var",
          signature  = c(x = "Stud"),
          definition = function(x) {

  df <- x@df

  if (df > 2) {
    return(df / (df - 2))
  } else {
    stop("Variance is undefined for Student's t distribution with
          no more than 2 df.")
  }

})

#' @rdname Stud
setMethod("sd",
          signature  = c(x = "Stud"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Stud
setMethod("skew",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@df > 3) {
    0
  }  else {
    stop("Skewness is undefined for Student's t distribution
    with no more than 3 df.")
  }

})

#' @rdname Stud
setMethod("kurt",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@df > 4) {
    6 / (x@df - 4)
  } else if (x@df > 2) {
    Inf
  } else {
    stop("Kurtosis is undefined for Student's t distribution
  with no more than 2 df.")
  }

})

#' @rdname Stud
setMethod("entro",
          signature  = c(x = "Stud"),
          definition = function(x) {

  df <- x@df
  ((df + 1) / 2) * (digamma((df + 1) / 2) - digamma(df / 2)) +
              log(df) / 2 + lbeta(df / 2, 1 / 2)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Stud
#' @export
llt <- function(x, df) {
  ll(Stud(df), x)
}

#' @rdname Stud
setMethod("ll",
          signature  = c(distr = "Stud", x = "numeric"),
          definition = function(distr, x) {

  df <- distr@df
  n <- length(x)
  s <- sum(log(1 + x ^ 2 / df))

  n * lgamma((df + 1) / 2) - n * lgamma(df / 2) - n * log(pi * df) / 2 -
    (df + 1) * s / 2

})
