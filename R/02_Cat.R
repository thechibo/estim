# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cat Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cat",
         contains = "Distribution",
         slots = c(prob = "numeric"),
         prototype = list(prob = c(0.5, 0.5)))

#' @title Categorical Distribution
#' @name Cat
#'
#' @param x an object of class `Cat`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Cat`.
#' @param prob numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#' @param dim numeric. The parameter dimension. See details.
#' @inherit Distributions return
#'
#' @details
#' The estimation of `prob` from a sample would by default return a vector of
#' probabilities corresponding to the categories that appeared in the sample and
#' 0 for the rest. However, the parameter dimension cannot be uncovered by the
#' sample, it has to be provided separately. This can be done with the argument
#' `dim`. If `dim` is not supplied, the dimension will be retrieved from the
#' `distr` argument. Categories that did not appear in the sample will have 0
#' probabilities appended to the end of the prob vector.
#'
#' Note that the actual dimension of the probability parameter vector is `k-1`,
#' therefore the Fisher information matrix and the asymptotic variance -
#' covariance matrix of the estimators is of dimension `(k-1)x(k-1)`.
#'
#' @importFrom extraDistr dcat rcat
#' @export
Cat <- function(prob = c(0.5, 0.5)) {
  new("Cat", prob = prob)
}

setValidity("Cat", function(object) {
  if(length(object@prob) <= 1) {
    stop("prob has to be a numeric of length at least 2")
  }
  if(any(object@prob <= 0) || any(object@prob >= 1)) {
    stop("prob has to be between 0 and 1")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
setMethod("d", signature = c(x = "Cat"),
          function(x) {
            function(y, log = FALSE) {
              dcat(y, prob = x@prob, log = log)
            }
          })

#' @rdname Cat
setMethod("r", signature = c(x = "Cat"),
          function(x) {
            function(n) {
              rcat(n, prob = x@prob)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
setMethod("mean",
          signature  = c(x = "Cat"),
          definition = function(x) {

  x@prob

})

#' @rdname Cat
setMethod("mode",
          signature  = c(x = "Cat"),
          definition = function(x) {

  which(x@prob == max(x@prob))

})

#' @rdname Cat
setMethod("var",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  diag(x@prob) - matrix(x@prob, k, 1) %*% matrix(x@prob, 1, k)

})

#' @rdname Cat
setMethod("entro",
          signature  = c(x = "Cat"),
          definition = function(x) {

  p <- x@prob

  - p * log(p) - (1 - p) * log(1 - p)

})

#' @rdname Cat
setMethod("finf",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  D <- diag(1 / x@prob[-k]) + matrix(1, k - 1, 1) %*% matrix(1, 1, k - 1) /
    x@prob[k]

  rownames(D) <- paste0("prob", seq_along(x@prob[-k]))
  colnames(D) <- paste0("prob", seq_along(x@prob[-k]))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llcat <- function(x, prob) {
  ll(x, prm = prob, distr = Cat())
}

#' @rdname Cat
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Cat"),
          definition = function(x, prm, distr) {

  sum(log(prm[x]))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
ecat <- function(x, type = "mle", ...) {

  estim(x, Cat(), type, ...)

}

#' @rdname Cat
setMethod("mle",
          signature  = c(x = "numeric", distr = "Cat"),
          definition = function(x, distr, dim = NULL) {

  if (is.null(dim)) {
    dim <- length(distr@prob)
  }

  p <- unname(table(x) / length(x))

  if (dim < length(p)) {
    stop("Dimension of Cat distribution supplied was ", dim, ", but ",
         length(p), " categories found in the sample.")
  }

  p <- c(p, rep(0, length = dim - length(p)))

  c(prob = p)

})

#' @rdname Cat
setMethod("me",
          signature  = c(x = "numeric", distr = "Cat"),
          definition = function(x, distr, dim = NULL) {

  mle(x, distr, dim)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname avar
#' @export
vcat <- function(prob, type = "mle") {

  avar(Cat(prob = prob), type = type)

}

#' @rdname Cat
setMethod("avar_mle",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  as.matrix(nearPD(solve(finf(distr))))

})

#' @rdname Cat
setMethod("avar_me",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  avar_mle(distr)

})
