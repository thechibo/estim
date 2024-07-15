# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cat Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cat",
         contains = "Distribution",
         slots = c(prob = "numeric"),
         prototype = list(prob = c(0.5, 0.5)))

#' @title Categorical Distribution
#' @name Cat
#'
#' @param x an object of class `Cat`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Cat`.
#' @param prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @importFrom extraDistr dcat rcat
#' @export
Cat <- function(prob = c(0.5, 0.5)) {
  new("Cat", prob = prob)
}

setValidity("Cat", function(object) {
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

#' @rdname Cat
setMethod("d", signature = c(x = "Cat"),
          function(x) {
            function(y, log = FALSE) {
              dcat(y, prob = x@prob, log = log)
            }
          })

#' @rdname Cat
setMethod("r", signature = c(x = "Cat"),
          function(x) {
            function(n) {
              rcat(n, prob = x@prob)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
setMethod("mean",
          signature  = c(x = "Cat"),
          definition = function(x) {

  x@prob

})

#' @rdname Cat
setMethod("var",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  diag(x@prob) - matrix(x@prob, k, 1) %*% matrix(x@prob, 1, k)

})

#' @rdname Cat
setMethod("finf",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  diag(1 / x@prob) - matrix(1, k, 1) %*% matrix(1, 1, k)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llcat <- function(x, prob) {
  ll(x, prm = prob, distr = Cat())
}

#' @rdname Cat
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Cat"),
          definition = function(x, prm, distr) {

  sum(log(prm[x]))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
ecat <- function(x, type = "mle", ...) {

  estim(x, Cat(), type, ...)

}

#' @rdname Cat
setMethod("mle",
          signature  = c(x = "numeric", distr = "Cat"),
          definition = function(x, distr) {

  c(prob = unname(table(x) / length(x)))

})

#' @rdname Cat
setMethod("me",
          signature  = c(x = "numeric", distr = "Cat"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vcat <- function(prob, type = "mle") {

  avar(Cat(prob = prob), type = type)

}

#' @rdname Cat
setMethod("avar_mle",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Cat
setMethod("avar_me",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  avar_mle(distr)

})
