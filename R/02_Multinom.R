# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multinom Distribution                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Multinom",
         contains = "Distribution",
         slots = c(size = "numeric", prob = "numeric"),
         prototype = list(size = 1, prob = c(0.5, 0.5)))

#' @title Multinomial Distribution
#' @name Multinom
#'
#' @param x an object of class `Multinom`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Multinom`.
#' @param size,prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Multinom <- function(size = 1, prob = c(0.5, 0.5)) {
  new("Multinom", size = size, prob = prob)
}

setValidity("Multinom", function(object) {
  if(length(object@size) != 1) {
    stop("size has to be a numeric of length 1")
  }
  if(!is_natural(object@size)) {
    stop("size has to be a natural number")
  }
  if(length(object@prob) > 1) {
    stop("prob has to be a numeric of length at least 2")
  }
  if(any(object@prob <= 0 || object@prob >= 1)) {
    stop("prob has to be between 0 and 1")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
setMethod("d", signature = c(x = "Multinom"),
          function(x) {
            function(y, log = FALSE) {
              dmultinom(y, size = x@size, prob = x@prob, log = log)
            }
          })

#' @rdname Multinom
setMethod("r", signature = c(x = "Multinom"),
          function(x) {
            function(n) {
              rmultinom(n, size = x@size, prob = x@prob)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
setMethod("mean",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  x@prob

})

#' @rdname Multinom
setMethod("var",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  k <- length(x@prob)

  x@size * (diag(x@prob) - matrix(x@prob, k, 1) %*% matrix(x@prob, 1, k))

})

#' @rdname Multinom
setMethod("finf",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  k <- length(x@prob)

  x@size * (diag(1 / x@prob) - matrix(1, k, 1) %*% matrix(1, 1, k))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llMultinom <- function(x, size, prob) {
  ll(x, prm = c(size, prob), distr = Multinom())
}

#' @rdname Multinom
setMethod("ll",
          signature  = c(x = "matrix", prm = "numeric", distr = "Multinom"),
          definition = function(x, prm, distr) {

  ncol(x) * lfactorial(prm[1]) - sum(lfactorial(x)) +
  sum(t(x) %*% diag(log(prm[-1])))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
emultinom <- function(x, type = "mle", ...) {

  estim(x, Multinom(), type, ...)

}

#' @rdname Multinom
setMethod("mle",
          signature  = c(x = "matrix", distr = "Multinom"),
          definition = function(x, distr) {

  c(prob = colMeans(x))

})

#' @rdname Multinom
setMethod("me",
          signature  = c(x = "matrix", distr = "Multinom"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vmultinom <- function(size, prob, type = "mle") {

  avar(Multinom(size = size, prob = prob), type = type)

}

#' @rdname Multinom
setMethod("avar_mle",
          signature  = c(distr = "Multinom"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Multinom
setMethod("avar_me",
          signature  = c(distr = "Multinom"),
          definition = function(distr) {

  avar_mle(distr)

})
