# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bern Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Bern",
         contains = "Distribution",
         slots = c(prob = "numeric"),
         prototype = list(prob = 0.5))

#' @title Bern Distribution
#' @name Bern
#'
#' @description
#' The Bernoulli distribution is a discrete probability distribution which takes
#' the value 1 with probability \eqn{p} and the value 0 with probability
#' \eqn{1 - p}, where \eqn{0 \leq p \leq 1}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Bern` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Bern` instead.
#' @param prob numeric. The distribution parameter, within the (0, 1) interval.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param log,log.p logical. Should the logarithm of the probability be returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the Bernoulli distribution is given
#' by: \deqn{ f(x; p) = p^x (1 - p)^{1 - x}, \quad p \in (0, 1), \quad x \in \{0, 1\}.}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dbinom()], [pbinom()], [qbinom()], [rbinom()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Bernoulli Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' p <- 0.7
#' D <- Bern(p)
#' x <- c(0, 1)
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
#' llbern(x, p)
#'
#' ebern(x, type = "mle")
#' ebern(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("bern", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vbern(p, type = "mle")
#' vbern(p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Bern <- function(prob = 0.5) {
  new("Bern", prob = prob)
}

setValidity("Bern", function(object) {
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

#' @rdname Bern
#' @export
dbern <- function(x, prob, log = FALSE) {
  dbinom(x, size = 1, prob, log)
}

#' @rdname Bern
#' @export
pbern <- function(x, prob, lower.tail = TRUE, log.p = FALSE) {
  pbinom(x, size = 1, prob, lower.tail, log.p)
}

#' @rdname Bern
#' @export
qbern <- function(x, prob, lower.tail = TRUE, log.p = FALSE) {
  qbinom(x, size = 1, prob, lower.tail, log.p)
}

#' @rdname Bern
#' @export
rbern <- function(n, prob) {
  rbinom(n, size = 1, prob)
}

#' @rdname Bern
setMethod("d", signature = c(distr = "Bern", x = "numeric"),
          function(distr, x) {
            dbinom(x, size = 1, prob = distr@prob)
          })

#' @rdname Bern
setMethod("p", signature = c(distr = "Bern", x = "numeric"),
          function(distr, x) {
            pbinom(x, size = 1, prob = distr@prob)
          })

#' @rdname Bern
setMethod("qn", signature = c(distr = "Bern", x = "numeric"),
          function(distr, x) {
            qbinom(x, size = 1, prob = distr@prob)
          })

#' @rdname Bern
setMethod("r", signature = c(distr = "Bern", n = "numeric"),
          function(distr, n) {
            rbinom(n, size = 1, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Bern
setMethod("mean",
          signature  = c(x = "Bern"),
          definition = function(x) {

  x@prob

})

#' @rdname Bern
setMethod("median",
          signature  = c(x = "Bern"),
          definition = function(x) {

  if (x@prob < 0.5) {
    return(0)
  } else if (x@prob > 0.5) {
    return(1)
  } else {
    warning("Bernoulli prob is equal to 0.5, therefore the median is any element
            in the [0, 1] interval. 0.5 is returned by default.")
    return(0.5)
  }

})

#' @rdname Bern
setMethod("mode",
          signature  = c(x = "Bern"),
          definition = function(x) {

  if (x@prob < 0.5) {
    return(0)
  } else if (x@prob > 0.5) {
    return(1)
  } else {
    warning("Bernoulli prob is equal to 0.5, therefore the mode is both 0 and 1.
            1 is returned by default.")
    return(1)
  }

})

#' @rdname Bern
setMethod("var",
          signature  = c(x = "Bern"),
          definition = function(x) {

  x@prob * (1 - x@prob)

})

#' @rdname Bern
setMethod("sd",
          signature  = c(x = "Bern"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Bern
setMethod("skew",
          signature  = c(x = "Bern"),
          definition = function(x) {

  p <- x@prob
  (1 - 2 * p) / sqrt(p * (1 - p))

})

#' @rdname Bern
setMethod("kurt",
          signature  = c(x = "Bern"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p
  (1 - 6 * p * q) / (p * q)

})

#' @rdname Bern
setMethod("entro",
          signature  = c(x = "Bern"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p
  - (q * log(q) + p * log(p))

})

#' @rdname Bern
setMethod("finf",
          signature  = c(x = "Bern"),
          definition = function(x) {

  1 / (x@prob * (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Bern
#' @export
llbern <- function(x, prob) {
  ll(distr = Bern(prob), x = x)
}

#' @rdname Bern
setMethod("ll",
          signature  = c(distr = "Bern", x = "numeric"),
          definition = function(distr, x) {

  p <- distr@prob
  n <- length(x)
  s <- sum(x)

  log(p) * s + log(1 - p) * (n - s)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Bern
#' @export
ebern <- function(x, type = "mle", ...) {

  e(Bern(), x = x, type = type, ...)

}

#' @rdname Bern
setMethod("mle",
          signature  = c(distr = "Bern", x = "numeric"),
          definition = function(distr, x) {

  list(prob = mean(x))

})

#' @rdname Bern
setMethod("me",
          signature  = c(distr = "Bern", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Bern
#' @export
vbern <- function(prob, type = "mle") {

  avar(Bern(prob), type = type)

}

#' @rdname Bern
setMethod("avar_mle",
          signature  = c(distr = "Bern"),
          definition = function(distr) {

  p <- distr@prob
  c(prob = p * (1 - p))

})

#' @rdname Bern
setMethod("avar_me",
          signature  = c(distr = "Bern"),
          definition = function(distr) {

  avar_mle(distr)

})
