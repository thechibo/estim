# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lnorm Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Lnorm",
         contains = "Distribution",
         slots = c(meanlog = "numeric", sdlog = "numeric"),
         prototype = list(meanlog = 0, sdlog = 1))

#' @title Log-Normal Distribution
#' @name Lnorm
#'
#' @description
#' The Lognormal distribution is an absolute continuous probability distribution
#' of a random variable whose logarithm is normally distributed. It is defined
#' by parameters \eqn{\mu} and \eqn{\sigma > 0}, which are the mean and standard
#' deviation of the underlying normal distribution.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Lnorm` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Lnorm` instead.
#' @param meanlog,sdlog numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Lognormal distribution is:
#' \deqn{ f(x; \mu, \sigma) = \frac{1}{x \sigma \sqrt{2\pi}} e^{-\frac{(\log x - \mu)^2}{2 \sigma^2}}, \quad x > 0 .}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dlnorm()], [plnorm()], [qlnorm()], [rlnorm()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Lnorm Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' m <- 3 ; s <- 5
#' D <- Lnorm(m, s)
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
#' elnorm(x, type = "mle")
#' elnorm(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("lnorm", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vlnorm(m, s, type = "mle")
#' vlnorm(m, s, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Lnorm <- function(meanlog = 0, sdlog = 1) {
  new("Lnorm", meanlog = meanlog, sdlog = sdlog)
}

setValidity("Lnorm", function(object) {
  if(length(object@meanlog) != 1) {
    stop("meanlog has to be a numeric of length 1")
  }
  if(length(object@sdlog) != 1) {
    stop("sdlog has to be a numeric of length 1")
  }
  if(object@sdlog <= 0) {
    stop("sdlog has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
setMethod("d", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            dlnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("p", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            plnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("qn", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            qlnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("r", signature = c(distr = "Lnorm", n = "numeric"),
          function(distr, n) {
            rlnorm(n, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
setMethod("mean",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog + x@sdlog ^ 2 / 2)

})

#' @rdname Lnorm
setMethod("median",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog)

})

#' @rdname Lnorm
setMethod("mode",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog - x@sdlog ^ 2)

})

#' @rdname Lnorm
setMethod("var",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  (exp(x@sdlog ^ 2) - 1) * exp(2 * x@meanlog + x@sdlog ^ 2)

})

#' @rdname Lnorm
setMethod("sd",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Lnorm
setMethod("skew",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  s <- x@sdlog
  (exp(s ^ 2) + 2) * sqrt(exp(s ^ 2) - 1)

})

#' @rdname Lnorm
setMethod("kurt",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  s <- x@sdlog
  exp(4 * s ^ 2) + 2 * exp(3 * s ^ 2) + 3 * exp(2 * s ^ 2) - 6

})

#' @rdname Lnorm
setMethod("entro",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  m <- x@meanlog
  s <- x@sdlog

  log(sqrt(2 * pi) * s * exp(m + 0.5), base = 2)

})

#' @rdname Lnorm
setMethod("finf",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  mat <- matrix(c(1, 0, 0, 2) / x@sdlog, 2, 2)
  prm_names <- c("meanlog", "sdlog")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
lllnorm <- function(x, meanlog, sdlog) {
  ll(Lnorm(meanlog, sdlog), x)
}

#' @rdname Lnorm
setMethod("ll",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

    m <- distr@meanlog
    s <- distr@sdlog
  - 0.5 * length(x) * log(2 * pi * s ^ 2) - sum(log(x)) -
    0.5 * sum((log(x) - m) ^ 2) / s ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
elnorm <- function(x, type = "mle", ...) {

  e(Lnorm(), x, type, ...)

}

#' @rdname Lnorm
setMethod("mle",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

  list(meanlog = mean(log(x)), sdlog = bsd(log(x)))

})

#' @rdname Lnorm
setMethod("me",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
vlnorm <- function(meanlog, sdlog, type = "mle") {

  avar(Lnorm(meanlog = meanlog, sdlog = sdlog), type = type)

}

#' @rdname Lnorm
setMethod("avar_mle",
          signature  = c(distr = "Lnorm"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Lnorm
setMethod("avar_me",
          signature  = c(distr = "Lnorm"),
          definition = function(distr) {

  avar_mle(distr)

})
