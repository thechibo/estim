# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lnorm Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Lnorm",
         contains = "Distribution",
         slots = c(meanlog = "numeric", sdlog = "numeric"),
         prototype = list(meanlog = 0, sdlog = 1))

#' @title Lnorm Distribution
#' @name Lnorm
#'
#' @param x an object of class `Lnorm`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Lnorm`.
#' @param meanlog,sdlog numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Lnorm <- function(meanlog = 0, sdlog = 1) {
  new("Lnorm", meanlog = meanlog, sdlog = sdlog)
}

setValidity("Lnorm", function(object) {
  if(length(object@meanlog) != 1) {
    stop("meanlog has to be a numeric of length 1")
  }
  if(length(object@sdlog) != 1) {
    stop("sdlog has to be a numeric of length 1")
  }
  if(object@sdlog <= 0) {
    stop("sdlog has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
setMethod("d", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            dlnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("p", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            plnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("qn", signature = c(distr = "Lnorm", x = "numeric"),
          function(distr, x) {
            qlnorm(x, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

#' @rdname Lnorm
setMethod("r", signature = c(distr = "Lnorm", n = "numeric"),
          function(distr, n) {
            rlnorm(n, meanlog = distr@meanlog, sdlog = distr@sdlog)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
setMethod("mean",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog + x@sdlog ^ 2 / 2)

})

#' @rdname Lnorm
setMethod("median",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog)

})

#' @rdname Lnorm
setMethod("mode",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@meanlog - x@sdlog ^ 2)

})

#' @rdname Lnorm
setMethod("var",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  (exp(x@sdlog ^ 2) - 1) * exp(2 * x@meanlog + x@sdlog ^ 2)

})

#' @rdname Lnorm
setMethod("sd",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Lnorm
setMethod("skew",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  s <- x@sdlog
  (exp(s ^ 2) + 2) * sqrt(exp(s ^ 2) - 1)

})

#' @rdname Lnorm
setMethod("kurt",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  s <- x@sdlog
  exp(4 * s ^ 2) + 2 * exp(3 * s ^ 2) + 3 * exp(2 * s ^ 2) - 6

})

#' @rdname Lnorm
setMethod("entro",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  m <- x@meanlog
  s <- x@sdlog

  log(sqrt(2 * pi) * s * exp(m + 0.5), base = 2)

})

#' @rdname Lnorm
setMethod("finf",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  mat <- matrix(c(1, 0, 0, 2) / x@sdlog, 2, 2)
  prm_names <- c("meanlog", "sdlog")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
lllnorm <- function(x, meanlog, sdlog) {
  ll(Lnorm(meanlog, sdlog), x)
}

#' @rdname Lnorm
setMethod("ll",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

    m <- distr@meanlog
    s <- distr@sdlog
  - 0.5 * length(x) * log(2 * pi * s ^ 2) - sum(log(x)) -
    0.5 * sum((log(x) - m) ^ 2) / s ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
elnorm <- function(x, type = "mle", ...) {

  e(Lnorm(), x, type, ...)

}

#' @rdname Lnorm
setMethod("mle",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

  list(meanlog = mean(log(x)), sdlog = bsd(log(x)))

})

#' @rdname Lnorm
setMethod("me",
          signature  = c(distr = "Lnorm", x = "numeric"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
#' @export
vlnorm <- function(meanlog, sdlog, type = "mle") {

  avar(Lnorm(meanlog = meanlog, sdlog = sdlog), type = type)

}

#' @rdname Lnorm
setMethod("avar_mle",
          signature  = c(distr = "Lnorm"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Lnorm
setMethod("avar_me",
          signature  = c(distr = "Lnorm"),
          definition = function(distr) {

  avar_mle(distr)

})
