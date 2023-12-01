# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exp Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llexp <- function(x, rate) {
  ll(x, prm = rate, distr = distr::Exp())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Exp"),
          definition = function(x, prm, distr) {

  sum(dexp(x = x, rate = prm, log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
eexp <- function(x, type = "mle", ...) {

  estim(x, distr::Exp(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Exp"),
          definition = function(x, distr) {

  c(rate = 1 / mean(x))

})

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Exp"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vexp <- function(rate, type = "mle") {

  avar(distr::Exp(rate = rate), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  rate <- distr::rate(distr)
  c(rate = rate ^ 2)

})

#' @rdname avar_me
setMethod("avar_me",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  avar_mle(distr)

})
