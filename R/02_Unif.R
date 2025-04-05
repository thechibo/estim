# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unif Distribution                                                         ----
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
#' @description
#' The Uniform distribution is an absolute continuous probability distribution
#' where all intervals of the same length within the distribution's support are
#' equally probable. It is defined by two parameters: the lower bound \eqn{a}
#' and the upper bound \eqn{b}, with \eqn{a < b}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Unif` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Unif` instead.
#' @param min,max numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Uniform distribution is:
#' \deqn{ f(x; a, b) = \frac{1}{b - a}, \quad a \le x \le b .}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dunif()], [punif()], [qunif()], [runif()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Uniform Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- 3 ; b <- 5
#' D <- Unif(a, b)
#' x <- c(0.3, 0.8, 0.5)
#' n <- 100
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, x) # density function
#' p(D, x) # distribution function
#' qn(D, x) # inverse distribution function
#' x <- r(D, n) # random generator function
#'
#' # alternative way to use the function
#' df <- d(D) ; df(x) # df is a function itself
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
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
#' llunif(x, a, b)
#'
#' eunif(x, type = "mle")
#' eunif(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("unif", x) # the distr argument can be a character
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
  if(object@min >= object@max) {
    stop("min must be less than max")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
setMethod("d", signature = c(distr = "Unif", x = "numeric"),
          function(distr, x) {
            dunif(x, min = distr@min, max = distr@max)
          })

#' @rdname Unif
setMethod("p", signature = c(distr = "Unif", x = "numeric"),
          function(distr, x) {
            punif(x, min = distr@min, max = distr@max)
          })

#' @rdname Unif
setMethod("qn", signature = c(distr = "Unif", x = "numeric"),
          function(distr, x) {
            qunif(x, min = distr@min, max = distr@max)
          })

#' @rdname Unif
setMethod("r", signature = c(distr = "Unif", n = "numeric"),
          function(distr, n) {
            runif(n, min = distr@min, max = distr@max)
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
setMethod("median",
          signature  = c(x = "Unif"),
          definition = function(x) {

            (x@max + x@min) / 2

          })

#' @rdname Unif
setMethod("mode",
          signature  = c(x = "Unif"),
          definition = function(x) {

            warning("The mode is any element in the support (or its interior) of
            a Uniform distribution. The mean is returned by default.")
            return((x@max + x@min) / 2)

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

#' @rdname Unif
#' @export
llunif <- function(x, min, max) {
  ll(Unif(min, max), x)
}

#' @rdname Unif
setMethod("ll",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x) {

  m <- distr@min
  M <- distr@max
  if (max(x) > M || min(x) < m) {
    return(0)
  } else {
    return(- length(x) * log(M - m))
  }

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
#' @export
eunif <- function(x, type = "mle", ...) {

  e(Unif(), x, type, ...)

}

#' @rdname Unif
setMethod("mle",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x) {

  list(min = min(x), max = max(x))

})

#' @rdname Unif
setMethod("me",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x) {

  m <- mean(x)
  s <- sqrt(3) * bsd(x)

  list(min = m - s, max = m + s)

})
