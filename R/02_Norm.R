# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Norm Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Norm",
         contains = "Distribution",
         slots = c(mean = "numeric", sd = "numeric"),
         prototype = list(mean = 0, sd = 1))

#' @title Normal Distribution
#' @name Norm
#'
#' @param x an object of class `Norm`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Norm`.
#' @param mean,sd numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Norm <- function(mean = 0, sd = 1) {
  new("Norm", mean = mean, sd = sd)
}

setValidity("Norm", function(object) {
  if(length(object@mean) != 1) {
    stop("mean has to be a numeric of length 1")
  }
  if(length(object@sd) != 1) {
    stop("sd has to be a numeric of length 1")
  }
  if(object@sd <= 0) {
    stop("sd has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("d", signature = c(x = "Norm"),
          function(x) {
            function(y, log = FALSE) {
              dnorm(y, mean = x@mean, sd = x@sd, log = log)
            }
          })

#' @rdname Norm
setMethod("p", signature = c(x = "Norm"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pnorm(q, mean = x@mean, sd = x@sd,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Norm
setMethod("q2", signature = c(x = "Norm"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qnorm(p, mean = x@mean, sd = x@sd,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Norm
setMethod("r", signature = c(x = "Norm"),
          function(x) {
            function(n) {
              rnorm(n, mean = x@mean, sd = x@sd)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("mean",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("median",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("mode",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("var",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd ^ 2

})

#' @rdname Norm
setMethod("sd",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd

})

#' @rdname Norm
setMethod("skew",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("kurt",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("entro",
          signature  = c(x = "Norm"),
          definition = function(x) {

  log(2 * pi * exp(1) * x@sd ^ 2) / 2

})

#' @rdname Norm
setMethod("finf",
          signature  = c(x = "Norm"),
          definition = function(x) {

  sd <- x@sd

  mat <- matrix(c(1 / sd ^ 2, 0, 0, 2 / sd ^ 2), 2, 2)
  prm_names <- c("mean", "sd")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llnorm <- function(x, mean, sd) {
  ll(x, prm = c(mean, sd), distr = Norm())
}

#' @rdname Norm
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Norm"),
          definition = function(x, prm, distr) {

  - 0.5 * length(x) * log(2 * pi * prm[2] ^ 2) -
              0.5 * sum((x - prm[1]) ^ 2) / prm[2] ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
enorm <- function(x, type = "mle", ...) {

  estim(x, Norm(), type, ...)

}

#' @rdname Norm
setMethod("mle",
          signature  = c(x = "numeric", distr = "Norm"),
          definition = function(x, distr) {

  c(mean = mean(x), sd = bsd(x))

})

#' @rdname Norm
setMethod("me",
          signature  = c(x = "numeric", distr = "Norm"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vnorm <- function(mean, sd, type = "mle") {

  avar(Norm(mean = mean, sd = sd), type = type)

}

#' @rdname Norm
setMethod("avar_mle",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Norm
setMethod("avar_me",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  avar_mle(distr)

})
