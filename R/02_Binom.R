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
#' @param distr an object of class `Binom`.
#' @param size,prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Binom"),
          function(x) {
            function(y, log = FALSE) {
              dbinom(y, size = x@size, prob = x@prob, log = log)
            }
          })

#' @rdname Binom
setMethod("p", signature = c(x = "Binom"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pbinom(q, size = x@size, prob = x@prob,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Binom
setMethod("q2", signature = c(x = "Binom"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qbinom(p, size = x@size, prob = x@prob,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Binom
setMethod("r", signature = c(x = "Binom"),
          function(x) {
            function(n) {
              rbinom(n, size = x@size, prob = x@prob)
            }
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
setMethod("finf",
          signature  = c(x = "Binom"),
          definition = function(x) {

  1 / (x@prob * (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llbinom <- function(x, size, prob) {
  ll(x, prm = c(size, prob), distr = Binom())
}

#' @rdname Binom
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Binom"),
          definition = function(x, prm, distr) {

  n <- length(x)
  s <- sum(x)
  y <- sum(unlist(lapply(x, FUN = function(x) { lchoose(prm[1], x) })))

  log(prm[2]) * s + log(1 - prm[2]) * (n * prm[1] - s) + y

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
ebinom <- function(x, type = "mle", ...) {

  estim(x, Binom(), type, ...)

}

#' @rdname Binom
setMethod("mle",
          signature  = c(x = "numeric", distr = "Binom"),
          definition = function(x, distr) {

  c(prob = mean(x))

})

#' @rdname Binom
setMethod("me",
          signature  = c(x = "numeric", distr = "Binom"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vbinom <- function(size, prob, type = "mle") {

  avar(Binom(size = size, prob = prob), type = type)

}

#' @rdname Binom
setMethod("avar_mle",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  prob <- distr@prob
  c(prob = prob * (1 - prob))

})

#' @rdname Binom
setMethod("avar_me",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  avar_mle(distr)

})
