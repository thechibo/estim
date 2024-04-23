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
#' @param distr an object of class `Laplace`.
#' @param mu,sigma numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Laplace"),
          function(x) {
            function(y, log = FALSE) {
              extraDistr::dlaplace(y, mu = x@mu, sigma = x@sigma, log = log)
            }
          })

#' @rdname Laplace
setMethod("p", signature = c(x = "Laplace"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              extraDistr::plaplace(q, mu = x@mu, sigma = x@sigma,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Laplace
setMethod("q2", signature = c(x = "Laplace"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              extraDistr::qlaplace(p, mu = x@mu, sigma = x@sigma,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Laplace
setMethod("r", signature = c(x = "Laplace"),
          function(x) {
            function(n) {
              extraDistr::rlaplace(n, mu = x@mu, sigma = x@sigma)
            }
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

#' @rdname ll
#' @export
lllaplace <- function(x, mu, sigma) {
  ll(x, prm = c(mu, sigma), distr = Laplace())
}

#' @rdname Laplace
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Laplace"),
          definition = function(x, prm, distr) {

  - length(x) * log(2 * prm[2]) - sum(abs(x - prm[1])) / prm[2]

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
elaplace <- function(x, type = "mle", ...) {

  estim(x, Laplace(), type, ...)

}

#' @rdname Laplace
setMethod("mle",
          signature  = c(x = "numeric", distr = "Laplace"),
          definition = function(x, distr) {

  m <- median(x)
  c(mu = m, sigma = mean(abs(x - m)))

})

#' @rdname Laplace
setMethod("me",
          signature  = c(x = "numeric", distr = "Laplace"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
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
