# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Binom Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llbinom <- function(x, size, prob) {
  ll(x, prm = c(size, prob), distr = distr::Binom())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Binom"),
          definition = function(x, prm, distr) {

  sum(dbinom(x = x, size = prm[1], prob = prm[2], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
ebinom <- function(x, type = "mle", ...) {

  estim(x, distr::Binom(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Binom"),
          definition = function(x, distr) {

  c(prob = mean(x))

})

#' @rdname me
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

  avar(distr::Binom(size = size, prob = prob), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  prob <- distr::prob(distr)
  c(prob = prob * (1 - prob))

})

#' @rdname avar_me
setMethod("avar_me",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  avar_mle(distr)

})
