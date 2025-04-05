# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Binom Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Binom",
         contains = "Distribution",
         slots = c(size = "numeric", prob = "numeric"),
         prototype = list(size = 1, prob = 0.5))

#' @title Binom Distribution
#' @name Binom
#'
#' @description
#' The binomial distribution is a discrete probability distribution which models
#' the probability of having x successes in n independent Binomoulli trials with
#' success probability p.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Binom` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Binom` instead.
#' @param size,prob numeric. The distribution parameters, `size` must be a
#' positive integer and `prob` must be within the (0, 1) interval.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the binomial distribution is given
#' by: \deqn{ f(x; n, p) = \binom{n}{x} p^x (1 - p)^{n - x}, \quad N \in
#' \mathbb{N}, \quad p \in (0, 1),} with \eqn{x \in \{0, 1, \dots, N\}}.
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
#' # Binomial Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' N <- 10 ; p <- 0.7
#' D <- Binom(N, p)
#' x <- 0:N
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
#' llbinom(x, N, p)
#'
#' ebinom(x, size = N, type = "mle")
#' ebinom(x, size = N, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vbinom(N, p, type = "mle")
#' vbinom(N, p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Binom <- function(size = 1, prob = 0.5) {
  new("Binom", size = size, prob = prob)
}

setValidity("Binom", function(object) {
  if(length(object@size) != 1) {
    stop("size has to be a numeric of length 1")
  }
  if(!is_natural(object@size)) {
    stop("size has to be a natural number")
  }
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

#' @rdname Binom
setMethod("d", signature = c(distr = "Binom", x = "numeric"),
          function(distr, x) {
            dbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Binom
setMethod("p", signature = c(distr = "Binom", x = "numeric"),
          function(distr, x) {
            pbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Binom
setMethod("qn", signature = c(distr = "Binom", x = "numeric"),
          function(distr, x) {
            qbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Binom
setMethod("r", signature = c(distr = "Binom", n = "numeric"),
          function(distr, n) {
            rbinom(n, size = distr@size, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
setMethod("mean",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size * x@prob

})

#' @rdname Binom
setMethod("var",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size * x@prob * (1 - x@prob)

})

#' @rdname Binom
setMethod("sd",
          signature  = c(x = "Binom"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Binom
setMethod("skew",
          signature  = c(x = "Binom"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p
  (q - p) / sqrt(x@size * p * q)

})

#' @rdname Binom
setMethod("kurt",
          signature  = c(x = "Binom"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p

  (1 - 6 * p * q) / (x@size * p * q)

})

#' @rdname Binom
setMethod("entro",
          signature  = c(x = "Binom"),
          definition = function(x) {

  warning("The entropy given is an approximation in the O(1 / n) order.")
  p <- x@prob
  0.5 * log(2 * pi  * exp(1) * x@size * p * (1 - p), base = 2)

})

#' @rdname Binom
setMethod("finf",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size / (x@prob * (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
llbinom <- function(x, size, prob) {
  ll(distr = Binom(size, prob), x)
}

#' @rdname Binom
setMethod("ll",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x) {

  N <- distr@size
  p <- distr@prob
  n <- length(x)
  s <- sum(x)
  y <- sum(unlist(lapply(x, FUN = function(x) { lchoose(N, x) })))

  log(p) * s + log(1 - p) * (n * N - s) + y

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
ebinom <- function(x, size, type = "mle", ...) {

  e(Binom(size = size), x, type, ...)

}

#' @rdname Binom
setMethod("mle",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x) {

  p <- mean(x) / distr@size

  if (p > 1) {
    stop("Success probability ", p, ", greater than 1.
          Did you forget to specify the size of the Binomial?")
  }

  list(prob = p)

})

#' @rdname Binom
setMethod("me",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
vbinom <- function(size, prob, type = "mle") {

  avar(Binom(size = size, prob = prob), type = type)

}

#' @rdname Binom
setMethod("avar_mle",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  prob <- distr@prob
  size <- distr@size
  c(prob = prob * (1 - prob) / size)

})

#' @rdname Binom
setMethod("avar_me",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  avar_mle(distr)

})
