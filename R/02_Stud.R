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
#' @description
#' The Student's t-distribution is a continuous probability distribution used
#' primarily in hypothesis testing and in constructing confidence intervals for
#' small sample sizes. It is defined by one parameter: the degrees of freedom
#' \eqn{\nu > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Stud` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Stud` instead.
#' @param df numeric. The distribution parameter.
#'
#' @details
#' The probability density function (PDF) of the Student's t-distribution is:
#' \deqn{ f(x; \nu) = \frac{\Gamma\left(\frac{\nu + 1}{2}\right)}{\sqrt{\nu\pi}\ \Gamma\left(\frac{\nu}{2}\right)}\left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu + 1}{2}} .}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dt()], [pt()], [qt()], [rt()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Student Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' df <- 12
#' D <- Stud(df)
#' x <- c(-3, 0, 3)
#' n <- 100
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, x) # density function
#' p(D, x) # distribution function
#' qn(D, 0.8) # inverse distribution function
#' x <- r(D, n) # random generator function
#'
#' # alternative way to use the function
#' d1 <- d(D) ; d1(x) # d1 is a function itself
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
#' median(D) # Median
#' mode(D) # Mode
#' var(D) # Variance
#' sd(D) # Standard Deviation
#' skew(D) # Skewness
#' kurt(D) # Excess Kurtosis
#' entro(D) # Entropy
#'
#' # List of all available moments
#' mom <- moments(D)
#' mom$mean # expectation
#'
#' # ------------------
#' # Point Estimation
#' # ------------------
#'
#' ll(D, x)
#' llt(x, df)
#'
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
