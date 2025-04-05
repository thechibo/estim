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
#' @description
#' The Exponential distribution is a continuous probability distribution often
#' used to model the time between independent events that occur at a constant
#' average rate. It is defined by the rate parameter \eqn{\lambda > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Exp` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Exp` instead.
#' @param rate numeric. The distribution rate parameter.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Exponential distribution is
#' given by: \deqn{ f(x; \lambda) = \lambda e^{-\lambda x}, \quad x \geq 0 .}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dexp()], [pexp()], [qexp()], [rexp()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Exp Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' rate <- 5
#' D <- Exp(rate)
#' x <- c(0.3, 2, 10)
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
#' df <- d(D) ; df(x) # df is a function itself
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
#' finf(D) # Fisher Information Matrix
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
#' llexp(x, rate)
#'
#' eexp(x, type = "mle")
#' eexp(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("exp", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vexp(rate, type = "mle")
#' vexp(rate, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
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
setMethod("d", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            dexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("p", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            pexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("qn", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x) {
            qexp(x, rate = distr@rate)
          })

#' @rdname Exp
setMethod("r", signature = c(distr = "Exp", n = "numeric"),
          function(distr, n) {
            rexp(n, rate = distr@rate)
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

#' @rdname Exp
#' @export
llexp <- function(x, rate) {
  ll(Exp(rate), x)
}

#' @rdname Exp
setMethod("ll",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  rate <- distr@rate
  length(x) * log(rate) - rate * sum(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
eexp <- function(x, type = "mle", ...) {

  e(Exp(), x, type, ...)

}

#' @rdname Exp
setMethod("mle",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  list(rate = 1 / mean(x))

})

#' @rdname Exp
setMethod("me",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
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
