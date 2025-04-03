# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cauchy Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cauchy",
         contains = "Distribution",
         slots = c(location = "numeric", scale = "numeric"),
         prototype = list(location = 0, scale = 1))

#' @title Cauchy Distribution
#' @name Cauchy
#'
#' @param x an object of class `Cauchy`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param n numeric. The sample size.
#' @param distr an object of class `Cauchy`.
#' @param location,scale numeric. The distribution parameters.
#' @param par0,method,lower,upper arguments passed to optim.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @inherit Distributions return
#'
#' @export
Cauchy <- function(location = 0, scale = 1) {
  new("Cauchy", location = location, scale = scale)
}

setValidity("Cauchy", function(object) {
  if(length(object@location) != 1) {
    stop("location has to be a numeric of length 1")
  }
  if(length(object@scale) != 1) {
    stop("scale has to be a numeric of length 1")
  }
  if(object@scale <= 0) {
    stop("scale has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
setMethod("d", signature = c(distr = "Cauchy", x = "numeric"),
          function(distr, x) {
          dcauchy(x, location = distr@location, scale = distr@scale)
          })

#' @rdname Cauchy
setMethod("p", signature = c(distr = "Cauchy", x = "numeric"),
          function(distr, x) {
            pcauchy(x, location = distr@location, scale = distr@scale)
          })

#' @rdname Cauchy
setMethod("qn", signature = c(distr = "Cauchy", x = "numeric"),
          function(distr, x) {
            qcauchy(x, location = distr@location, scale = distr@scale)
          })

#' @rdname Cauchy
setMethod("r", signature = c(distr = "Cauchy", n = "numeric"),
          function(distr, n) {
            rcauchy(n, location = distr@location, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
setMethod("mean",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The mean of the Cauchy distribution is not defined. NaN is returned")
  NaN

})


#' @rdname Cauchy
setMethod("median",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  x@location

})

#' @rdname Cauchy
setMethod("mode",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  x@location

})

#' @rdname Cauchy
setMethod("var",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The variance of the Cauchy distribution is not defined.
          NaN is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("sd",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The standard deviation of the Cauchy distribution is not
          defined. NaN is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("skew",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The skewness of the Cauchy distribution is not defined. NaN
          is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("kurt",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The kurtosis of the Cauchy distribution is not defined. NaN
          is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("entro",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  log(4 * pi * x@scale)

})

#' @rdname Cauchy
setMethod("finf",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  D <- diag(2) / (2 * x@scale ^ 2)

  rownames(D) <- c("location", "scale")
  colnames(D) <- c("location", "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
llcauchy <- function(x, location, scale) {
  ll(distr = Cauchy(location, scale), x)
}

#' @rdname Cauchy
setMethod("ll",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x) {


  - length(x) * log(pi * distr@scale) -
      sum(log (1 + ((x - distr@location) / distr@scale) ^ 2))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Cauchy"),
          definition = function(par, tx, distr) {

  log(par[2]) - mean(log((tx - par[1]) ^ 2 + par[2] ^ 2))

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Cauchy"),
          definition = function(par, tx, distr) {

  c(2 * mean((tx - par[1]) / ((tx - par[1]) ^ 2 + par[2] ^ 2)),
    1 / par[2] - 2 * par[2] * mean(1 / ((tx - par[1]) ^ 2 + par[2] ^ 2)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
ecauchy <- function(x, type = "mle", ...) {

  e(Cauchy(), x, type, ...)

}

#' @rdname Cauchy
setMethod("mle",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x,
                                par0 = "me",
                                method = "L-BFGS-B",
                                lower = c(-Inf, 1e-5),
                                upper = c(Inf, Inf)) {

  par <- optim(par = unlist(do.call(par0, list(distr = distr, x = x))),
               fn = lloptim,
               gr = dlloptim,
               tx = x,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  par <- unname(par)
  list(location = par[1], scale = par[2])

})

#' @rdname Cauchy
#' @export
setMethod("me",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x) {

  x0 <- median(x)

  list(location = x0, scale = median(abs(x - x0)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
vcauchy <- function(location, scale, type = "mle") {

  avar(Cauchy(location = location, scale = scale), type = type)

}

#' @rdname Cauchy
setMethod("avar_mle",
          signature  = c(distr = "Cauchy"),
          definition = function(distr) {

  inv2x2(finf(distr))

})
