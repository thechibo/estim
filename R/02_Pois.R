# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pois Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llpois <- function(x, lambda) {
  ll(x, prm = lambda, distr = distr::Pois())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Pois"),
          definition = function(x, prm, distr) {

  sum(dpois(x = x, lambda = prm, log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
epois <- function(x, type = "mle", ...) {

  estim(x, distr::Pois(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Pois"),
          definition = function(x, distr) {

  c(lambda = mean(x))

})

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Pois"),
          definition = function(x, distr) {

  mle(x, distr)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vpois <- function(lambda, type = "mle") {

  avar(distr::Pois(lambda = lambda), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
          signature  = c(distr = "Pois"),
          definition = function(distr) {

  c(lambda = distr::lambda(distr))

})

#' @rdname avar_me
setMethod("avar_me",
          signature  = c(distr = "Pois"),
          definition = function(distr) {

  avar_mle(distr)

})
