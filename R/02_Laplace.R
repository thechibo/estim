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
#' @param x an object of class `Laplace`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Laplace`.
#' @param mu,sigma numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @importFrom extraDistr dlaplace plaplace qlaplace rlaplace
#' @export
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
setMethod("d", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            extraDistr::dlaplace(x, mu = distr@mu, sigma = distr@sigma)
          })

#' @rdname Laplace
setMethod("p", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            extraDistr::plaplace(x, mu = distr@mu, sigma = distr@sigma)
          })


#' @rdname Laplace
setMethod("qn", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x) {
            extraDistr::qlaplace(x, mu = distr@mu, sigma = distr@sigma)
          })


#' @rdname Laplace
setMethod("r", signature = c(distr = "Laplace", n = "numeric"),
          function(distr, n) {
            extraDistr::rlaplace(n, mu = distr@mu, sigma = distr@sigma)
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
