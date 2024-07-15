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
#' @param distr an object of class `Nbinom`.
#' @param size,prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Nbinom"),
          function(x) {
            function(y, log = FALSE) {
              dnbinom(y, size = x@size, prob = x@prob, log = log)
            }
          })

#' @rdname Nbinom
setMethod("p", signature = c(x = "Nbinom"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pnbinom(q, size = x@size, prob = x@prob,
                     lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Nbinom
setMethod("qn", signature = c(x = "Nbinom"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qnbinom(p, size = x@size, prob = x@prob,
                     lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Nbinom
setMethod("r", signature = c(x = "Nbinom"),
          function(x) {
            function(n) {
              rnbinom(n, size = x@size, prob = x@prob)
            }
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

#' @rdname ll
#' @export
llnbinom <- function(x, size, prob) {
  ll(x, prm = c(size, prob), distr = Nbinom())
}

#' @rdname Nbinom
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Nbinom"),
          definition = function(x, prm, distr) {

  n <- length(x)
  s <- sum(x)
  y <- sum(unlist(lapply(x, FUN = function(x) { lchoose(x + prm[1] - 1, x) })))

  log(1 - prm[2]) * s + n * prm[1] * log(prm[2]) + y

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
enbinom <- function(x, type = "mle", ...) {

  estim(x, Nbinom(), type, ...)

}

#' @rdname Nbinom
setMethod("mle",
          signature  = c(x = "numeric", distr = "Nbinom"),
          definition = function(x, distr) {

  size <- distr@size
  c(prob = size / (size + mean(x)))

})

#' @rdname Nbinom
setMethod("me",
          signature  = c(x = "numeric", distr = "Nbinom"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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
