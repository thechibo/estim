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
#' @description
#' The Geometric distribution is a discrete probability distribution that models
#' the number of failures before the first success in a sequence of independent
#' Bernoulli trials, each with success probability \eqn{0 < p \leq 1}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Geom` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Geom` instead.
#' @param prob numeric. The distribution parameter.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the Geometric distribution is:
#' \deqn{ P(X = k) = (1 - p)^k p, \quad k \in \mathbb{N}_0.}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dgeom()], [pgeom()], [qgeom()], [rgeom()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Geom Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' p <- 0.4
#' D <- Geom(p)
#' x <- 0:4
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
#' llgeom(x, p)
#'
#' egeom(x, type = "mle")
#' egeom(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("geom", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vgeom(p, type = "mle")
#' vgeom(p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
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
setMethod("d", signature = c(distr = "Geom", x = "numeric"),
          function(distr, x) {
            dgeom(x, prob = distr@prob)
          })

#' @rdname Geom
setMethod("p", signature = c(distr = "Geom", x = "numeric"),
          function(distr, x) {
            pgeom(x, prob = distr@prob)
          })

#' @rdname Geom
setMethod("qn", signature = c(distr = "Geom", x = "numeric"),
          function(distr, x) {
            qgeom(x, prob = distr@prob)
          })

#' @rdname Geom
setMethod("r", signature = c(distr = "Geom", n = "numeric"),
          function(distr, n) {
            rgeom(n, prob = distr@prob)
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
setMethod("median",
          signature  = c(x = "Geom"),
          definition = function(x) {

  y <- - 1 / log(1 - x@prob, base = 2)
  if ((y %% 1) == 0 & x@prob != 0.5) {
    warning("The median of the Geom distribution is not uniquely defined in this case.")
  }

  ceiling(y) - 1

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

  1 / (x@prob ^ 2 * (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Geom
#' @export
llgeom <- function(x, prob) {
  ll(Geom(prob), x)
}

#' @rdname Geom
setMethod("ll",
          signature  = c(distr = "Geom", x = "numeric"),
          definition = function(distr, x) {

  p <- distr@prob
  log(1 - p) * sum(x) + log(p) * length(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Geom
#' @export
egeom <- function(x, type = "mle", ...) {

  e(Geom(), x, type, ...)

}

#' @rdname Geom
setMethod("mle",
          signature  = c(distr = "Geom", x = "numeric"),
          definition = function(distr, x) {

  list(prob = 1 / (1 + mean(x)))

})

#' @rdname Geom
setMethod("me",
          signature  = c(distr = "Geom", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Geom
#' @export
vgeom <- function(prob, type = "mle") {

  avar(Geom(prob = prob), type = type)

}

#' @rdname Geom
setMethod("avar_mle",
          signature  = c(distr = "Geom"),
          definition = function(distr) {

  prob <- distr@prob
  c(prob = prob ^ 2 * (1 - prob))

})

#' @rdname Geom
setMethod("avar_me",
          signature  = c(distr = "Geom"),
          definition = function(distr) {

  avar_mle(distr)

})
