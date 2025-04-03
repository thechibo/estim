# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nbinom Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Nbinom",
         contains = "Distribution",
         slots = c(size = "numeric", prob = "numeric"),
         prototype = list(size = 1, prob = 0.5))

#' @title Negative Binomial Distribution
#' @name Nbinom
#'
#' @param x an object of class `Nbinom`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Nbinom`.
#' @param size,prob numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @importFrom stats integrate
#'
#' @inherit Distributions return
#'
#' @export
Nbinom <- function(size = 1, prob = 0.5) {
  new("Nbinom", size = size, prob = prob)
}

setValidity("Nbinom", function(object) {
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

#' @rdname Nbinom
setMethod("d", signature = c(distr = "Nbinom", x = "numeric"),
          function(distr, x) {
            dnbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Nbinom
setMethod("p", signature = c(distr = "Nbinom", x = "numeric"),
          function(distr, x) {
            pnbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Nbinom
setMethod("qn", signature = c(distr = "Nbinom", x = "numeric"),
          function(distr, x) {
            qnbinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Nbinom
setMethod("r", signature = c(distr = "Nbinom", n = "numeric"),
          function(distr, n) {
            rnbinom(n, size = distr@size, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Nbinom
setMethod("mean",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  x@size * (1 - x@prob) / x@prob

})

#' @rdname Nbinom
setMethod("median",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  qnbinom(0.5, size = x@size, prob = x@prob)

})

#' @rdname Nbinom
setMethod("mode",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  floor((x@size - 1) * (1 - x@prob) / x@prob)

})

#' @rdname Nbinom
setMethod("var",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  x@size * (1 - x@prob) / x@prob ^ 2

})

#' @rdname Nbinom
setMethod("sd",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Nbinom
setMethod("skew",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  (2 - x@prob) / sqrt((1 - x@prob) * x@size)

})

#' @rdname Nbinom
setMethod("kurt",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  k <- x@size
  p <- x@prob

  6 / k + p ^ 2 / ((1 - p) * k)

})

#' @rdname Nbinom
setMethod("entro",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  # https://arxiv.org/pdf/1708.06394
  # Expressions for the Entropy of Binomial-Type Distributions
  k <- x@size
  p <- x@prob
  h <- - p * log(p) - (1 - p) * log(1 - p)

  f <- function(z) {
    ((1 - z)^(k - 1) - 1) * ((1 + p*z / (1 - p)) ^ (- k) + p*k*z / (1 - p) - 1) /
      z * log(1 - z)
  }

  c <- integrate(f, lower = 0, upper = 1)$value

  k * (h - p * log(k)) / (1 - p) + c

})

#' @rdname Nbinom
setMethod("finf",
          signature  = c(x = "Nbinom"),
          definition = function(x) {

  size <- x@size
  prob <- x@prob
  c(prob = size / (prob ^ 2 * (1 - prob)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Nbinom
#' @export
llnbinom <- function(x, size, prob) {
  ll(Nbinom(size, prob), x)
}

#' @rdname Nbinom
setMethod("ll",
          signature  = c(distr = "Nbinom", x = "numeric"),
          definition = function(distr, x) {

  N <- distr@size
  p <- distr@prob

  n <- length(x)
  s <- sum(x)
  y <- sum(unlist(lapply(x, FUN = function(x) { lchoose(x + N - 1, x) })))

  log(1 - p) * s + n * N * log(p) + y

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Nbinom
#' @export
enbinom <- function(x, size, type = "mle", ...) {

  e(Nbinom(size = size), x, type, ...)

}

#' @rdname Nbinom
setMethod("mle",
          signature  = c(distr = "Nbinom", x = "numeric"),
          definition = function(distr, x) {

  size <- distr@size

  list(prob = size / (size + mean(x)))

})

#' @rdname Nbinom
setMethod("me",
          signature  = c(distr = "Nbinom", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Nbinom
#' @export
vnbinom <- function(size, prob, type = "mle") {

  avar(Nbinom(size = size, prob = prob), type = type)

}

#' @rdname Nbinom
setMethod("avar_mle",
          signature  = c(distr = "Nbinom"),
          definition = function(distr) {

  size <- distr@size
  prob <- distr@prob
  c(prob = prob ^ 2 * (1 - prob) / size)

})

#' @rdname Nbinom
setMethod("avar_me",
          signature  = c(distr = "Nbinom"),
          definition = function(distr) {

  avar_mle(distr)

})
