# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Laplace Distribution                                                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Laplace",
         contains = "Distribution",
         slots = c(mu = "numeric", sigma = "numeric"),
         prototype = list(mu = 0, sigma = 1))

#' @title Laplace Distribution
#' @name Laplace
#'
#' @description
#' The Laplace distribution, also known as the double exponential distribution,
#' is a continuous probability distribution that is often used to model data
#' with sharp peaks and heavy tails. It is parameterized by a location parameter
#' \eqn{\mu} and a scale parameter \eqn{b > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Laplace` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Laplace` instead.
#' @param mu,sigma numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Laplace distribution is:
#' \deqn{ f(x; \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right) .}
#'
#' @inherit Distributions return
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Laplace Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' m <- 3 ; s <- 5
#' D <- Laplace(m, s)
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
#' elaplace(x, type = "mle")
#' elaplace(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("laplace", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vlaplace(m, s, type = "mle")
#' vlaplace(m, s, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Laplace <- function(mu = 0, sigma = 1) {
  new("Laplace", mu = mu, sigma = sigma)
}

setValidity("Laplace", function(object) {
  if(length(object@mu) != 1) {
    stop("mu has to be a numeric of length 1")
  }
  if(length(object@mu) != 1) {
    stop("mu has to be a numeric of length 1")
  }
  if(object@sigma <= 0) {
    stop("sigma has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
dlaplace <- function(x, mu, sigma) {
  0.5 * dexp(abs(x - mu), rate = 1 / sigma)
}

plaplace <- function(x, mu, sigma) {
  unlist(lapply(x, function(x) {
    if (x >= mu) {
      return(1 - 0.5 * exp((mu - x) / sigma))
    } else {
      return(0.5 * exp((x - mu) / sigma))
    }
  }))
}

qlaplace <- function(x, mu, sigma) {
  unlist(lapply(x, function(x) {
    if (x >= 0.5) {
      return(mu + qexp(2 * (x - 0.5), rate = 1 / sigma))
    } else {
      return(mu - qexp(2 * x, rate = 1 / sigma, lower.tail = FALSE))
    }
  }))
}

rlaplace <- function(n, mu, sigma) {
  (2 * rbern(n, 0.5) - 1) * rexp(n, rate = 1 / sigma) + mu
}

#' @rdname Laplace
setMethod("d", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            dlaplace(x, distr@mu, distr@sigma)
          })

#' @rdname Laplace
setMethod("p", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            plaplace(x, distr@mu, distr@sigma)
          })

#' @rdname Laplace
setMethod("qn", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            qlaplace(x, distr@mu, distr@sigma)
          })

#' @rdname Laplace
setMethod("r", signature = c(distr = "Laplace", n = "numeric"),
          function(distr, n) {
            rlaplace(n, distr@mu, distr@sigma)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
setMethod("mean",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("median",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("mode",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("var",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  2 * x@sigma ^ 2

})

#' @rdname Laplace
setMethod("sd",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Laplace
setMethod("skew",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  0

})

#' @rdname Laplace
setMethod("kurt",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  6

})

#' @rdname Laplace
setMethod("entro",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  log(2 * x@sigma * exp(1))

})

#' @rdname Laplace
setMethod("finf",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  mat <- matrix(c(1, 0, 0, 1 / x@sigma), 2, 2)
  prm_names <- c("mu", "sigma")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
lllaplace <- function(x, mu, sigma) {
  ll(distr = Laplace(mu, sigma), x)
}

#' @rdname Laplace
setMethod("ll",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x) {

  - length(x) * log(2 * distr@sigma) - sum(abs(x - distr@mu)) / distr@sigma

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
elaplace <- function(x, type = "mle", ...) {

  e(Laplace(), x, type, ...)

}

#' @rdname Laplace
setMethod("mle",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x) {

  m <- median(x)

  list(mu = m, sigma = mean(abs(x - m)))

})

#' @rdname Laplace
setMethod("me",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
vlaplace <- function(mu, sigma, type = "mle") {

  avar(Laplace(mu = mu, sigma = sigma), type = type)

}

#' @rdname Laplace
setMethod("avar_mle",
          signature  = c(distr = "Laplace"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Laplace
setMethod("avar_me",
          signature  = c(distr = "Laplace"),
          definition = function(distr) {

  avar_mle(distr)

})
