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
#' @param distr an object of class `Lnorm`.
#' @param meanlog,sdlog numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
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
setMethod("d", signature = c(x = "Lnorm"),
          function(x) {
            function(y, log = FALSE) {
              dlnorm(y, meanlog = x@meanlog, sdlog = x@sdlog, log = log)
            }
          })

#' @rdname Lnorm
setMethod("p", signature = c(x = "Lnorm"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              plnorm(q, meanlog = x@meanlog, sdlog = x@sdlog,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Lnorm
setMethod("qn", signature = c(x = "Lnorm"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qlnorm(p, meanlog = x@meanlog, sdlog = x@sdlog,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Lnorm
setMethod("r", signature = c(x = "Lnorm"),
          function(x) {
            function(n) {
              rlnorm(n, meanlog = x@meanlog, sdlog = x@sdlog)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Lnorm
setMethod("mean",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@logmean + x@logsd ^ 2 / 2)

})

#' @rdname Lnorm
setMethod("median",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@logmean)

})

#' @rdname Lnorm
setMethod("mode",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  exp(x@logmean - x@logsd ^ 2)

})

#' @rdname Lnorm
setMethod("var",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  (exp(x@logsd ^ 2) - 1) * exp(2 * x@logmean + x@logsd ^ 2)

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

  s <- x@logsd
  (exp(s ^ 2) + 2) * sqrt(exp(s ^ 2) - 1)

})

#' @rdname Lnorm
setMethod("kurt",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  s <- x@logsd
  exp(4 * s ^ 2) + 2 * exp(3 * s ^ 2) + 3 * exp(2 * s ^ 2) - 6

})

#' @rdname Lnorm
setMethod("entro",
          signature  = c(x = "Lnorm"),
          definition = function(x) {

  m <- x@logmean
  s <- x@logsd

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

#' @rdname ll
#' @export
lllnorm <- function(x, meanlog, sdlog) {
  ll(x, prm = c(meanlog, sdlog), distr = Lnorm())
}

#' @rdname Lnorm
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Lnorm"),
          definition = function(x, prm, distr) {

  - 0.5 * length(x) * log(2 * pi * prm[2] ^ 2) - sum(log(x)) -
    0.5 * sum((log(x) - prm[1]) ^ 2) / prm[2] ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
elnorm <- function(x, type = "mle", ...) {

  estim(x, Lnorm(), type, ...)

}

#' @rdname Lnorm
setMethod("mle",
          signature  = c(x = "numeric", distr = "Lnorm"),
          definition = function(x, distr) {

  c(meanlog = mean(log(x)), sdlog = bsd(log(x)))

})

#' @rdname Lnorm
setMethod("me",
          signature  = c(x = "numeric", distr = "Lnorm"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
