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

#' @title Binomial Distribution
#' @name Binom
#'
#' @param x an object of class `Binom`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Binom`.
#' @param size,prob numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
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
