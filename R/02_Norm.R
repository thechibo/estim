# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Norm Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llnorm <- function(x, mean, sd) {
  ll(x, prm = c(mean, sd), distr = distr::Norm())
}

#' @rdname ll
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Norm"),
          definition = function(x, prm, distr) {

  sum(dnorm(x = x, mean = prm[1], sd = prm[2], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estim
#' @export
enorm <- function(x, type = "mle", ...) {

  estim(x, distr::Norm(), type, ...)

}

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Norm"),
          definition = function(x, distr) {

  c(mean = mean(x), sd = bsd(x))

})

#' @rdname me
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

  avar(distr::Norm(mean = mean, sd = sd), type = type)

}

#' @rdname avar_mle
setMethod("avar_mle",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  mu <- distr::mean(distr)
  sd <- distr::sd(distr)

  mat <- matrix(c(sd ^ 2, 0, 0, sd ^ 2 / 2), 2, 2)
  prm_names <- c("mean", "sd")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})

#' @rdname avar_me
setMethod("avar_me",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  avar_mle(distr)

})
