# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Norm Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Norm",
         contains = "Distribution",
         slots = c(mean = "numeric", sd = "numeric"),
         prototype = list(mean = 0, sd = 1))

#' @title Normal Distribution
#' @name Norm
#'
#' @description
#' The Normal or Gaussian distribution, is an absolute continuous probability
#' distribution characterized by two parameters: the mean \eqn{\mu} and the
#' standard deviation \eqn{\sigma > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Norm` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Norm` instead.
#' @param mean,sd numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Normal distribution is:
#' \deqn{ f(x; \mu, \sigma) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{1}{2} \left(\frac{x - \mu}{\sigma}\right)^2} .}
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dnorm()], [pnorm()], [qnorm()], [rnorm()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Normal Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' m <- 3 ; s <- 5
#' D <- Norm(m, s)
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
#' enorm(x, type = "mle")
#' enorm(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("norm", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vnorm(m, s, type = "mle")
#' vnorm(m, s, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Norm <- function(mean = 0, sd = 1) {
  new("Norm", mean = mean, sd = sd)
}

setValidity("Norm", function(object) {
  if(length(object@mean) != 1) {
    stop("mean has to be a numeric of length 1")
  }
  if(length(object@sd) != 1) {
    stop("sd has to be a numeric of length 1")
  }
  if(object@sd <= 0) {
    stop("sd has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("d", signature = c(distr = "Norm", x = "numeric"),
          function(distr, x) {
            dnorm(x, mean = distr@mean, sd = distr@sd)
          })

#' @rdname Norm
setMethod("p", signature = c(distr = "Norm", x = "numeric"),
          function(distr, x) {
            pnorm(x, mean = distr@mean, sd = distr@sd)
          })

#' @rdname Norm
setMethod("qn", signature = c(distr = "Norm", x = "numeric"),
          function(distr, x) {
            qnorm(x, mean = distr@mean, sd = distr@sd)
          })

#' @rdname Norm
setMethod("r", signature = c(distr = "Norm", n = "numeric"),
          function(distr, n) {
            rnorm(n, mean = distr@mean, sd = distr@sd)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("mean",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("median",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("mode",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("var",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd ^ 2

})

#' @rdname Norm
setMethod("sd",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd

})

#' @rdname Norm
setMethod("skew",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("kurt",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("entro",
          signature  = c(x = "Norm"),
          definition = function(x) {

  log(2 * pi * exp(1) * x@sd ^ 2) / 2

})

#' @rdname Norm
setMethod("finf",
          signature  = c(x = "Norm"),
          definition = function(x) {

  sd <- x@sd

  mat <- matrix(c(1 / sd ^ 2, 0, 0, 2 / sd ^ 2), 2, 2)
  prm_names <- c("mean", "sd")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
llnorm <- function(x, mean, sd) {
  ll(Norm(mean, sd), x)
}

#' @rdname Norm
setMethod("ll",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x) {

  m <- distr@mean
  s <- distr@sd

  - 0.5 * length(x) * log(2 * pi * s ^ 2) -
              0.5 * sum((x - m) ^ 2) / s ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
enorm <- function(x, type = "mle", ...) {

  e(Norm(), x, type, ...)

}

#' @rdname Norm
setMethod("mle",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x) {

  list(mean = mean(x), sd = bsd(x))

})

#' @rdname Norm
setMethod("me",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
vnorm <- function(mean, sd, type = "mle") {

  avar(Norm(mean = mean, sd = sd), type = type)

}

#' @rdname Norm
setMethod("avar_mle",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Norm
setMethod("avar_me",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  avar_mle(distr)

})
