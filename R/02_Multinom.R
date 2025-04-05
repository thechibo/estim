# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multinom Distribution                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Multinom",
         contains = "Distribution",
         slots = c(size = "numeric", prob = "numeric"),
         prototype = list(size = 1, prob = c(0.5, 0.5)))

#' @title Multinomial Distribution
#' @name Multinom
#'
#' @description
#' The multinomial distribution is a discrete probability distribution which models
#' the probability of having x successes in n independent categorical trials
#' with success probability vector p.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Multinom` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Multinom` instead.
#' @param size,prob numeric. The distribution parameters, `size` must be a
#' positive integer and `prob` must be a probability vector.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the Multinomial distribution is:
#' \deqn{ P(X_1 = x_1, ..., X_k = x_k) = \frac{n!}{x_1! x_2! ... x_k!} \prod_{i=1}^k p_i^{x_i}, }
#' subject to \eqn{ \sum_{i=1}^{k} x_i = n }.
#'
#' @inherit Distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dmultinom()], [rmultinom()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Multinomial Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' N <- 10 ; p <- c(0.1, 0.2, 0.7)
#' D <- Multinom(N, p)
#' x <- c(2, 3, 5)
#' n <- 100
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, x) # density function
#'
#' # alternative way to use the function
#' df <- d(D) ; df(x) # df is a function itself
#'
#' x <- r(D, n) # random generator function
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
#' mode(D) # Mode
#' var(D) # Variance
#' entro(D) # Entropy
#' finf(D) # Fisher Information Matrix
#'
#' # List of all available moments
#' mom <- moments(D)
#' mom$mean # expectation
#'
#' # ------------------
#' # Point Estimation
#' # ------------------
#'
#' ll(D, x)
#' llmultinom(x, N, p)
#'
#' emultinom(x, type = "mle")
#' emultinom(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("multinom", x) # the distr argument can be a character
#'
#' # ------------------
#' # As. Variance
#' # ------------------
#'
#' vmultinom(N, p, type = "mle")
#' vmultinom(N, p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' avar(D, type = "mle")
Multinom <- function(size = 1, prob = c(0.5, 0.5)) {
  new("Multinom", size = size, prob = prob)
}

setValidity("Multinom", function(object) {
  if(length(object@size) != 1) {
    stop("size has to be a numeric of length 1")
  }
  if(!is_natural(object@size)) {
    stop("size has to be a natural number")
  }
  if(length(object@prob) < 2) {
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

#' @rdname Multinom
setMethod("d", signature = c(distr = "Multinom", x = "numeric"),
          function(distr, x) {
            dmultinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Multinom
setMethod("d", signature = c(distr = "Multinom", x = "matrix"),
          function(distr, x) {
            dmultinom(x, size = distr@size, prob = distr@prob)
          })

#' @rdname Multinom
setMethod("r", signature = c(distr = "Multinom", n = "numeric"),
          function(distr, n) {
            rmultinom(n, size = distr@size, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
setMethod("mean",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  x@size * x@prob

})

#' @rdname Multinom
setMethod("mode",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  which(x@prob == max(x@prob))

})

#' @rdname Multinom
setMethod("var",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  k <- length(x@prob)

  x@size * (diag(x@prob) - matrix(x@prob, k, 1) %*% matrix(x@prob, 1, k))

})

#' @rdname Multinom
setMethod("entro",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  N <- x@size
  p <- x@prob

  z <- 0
  for (x in 0:N) {
    z <- z + choose(N, x) * p ^ x * (1 - p) ^ (N - x) * log(factorial(x))
  }

  - log(factorial(N)) - N * sum(p * log(p)) + sum(z)

})

#' @rdname Multinom
setMethod("finf",
          signature  = c(x = "Multinom"),
          definition = function(x) {

  k <- length(x@prob)

  if (k == 1) {
    y <- 1 / x@prob[-k]
  } else {
    y <- diag(1 / x@prob[-k])
  }

  D <- x@size * (y - matrix(1, k - 1, 1) %*% matrix(1, 1, k - 1) /
              x@prob[k])

  rownames(D) <- paste0("prob", seq_along(x@prob[-k]))
  colnames(D) <- paste0("prob", seq_along(x@prob[-k]))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
#' @export
llmultinom <- function(x, size, prob) {
  ll(distr = Multinom(size, prob), x)
}

#' @rdname Multinom
setMethod("ll",
          signature  = c(distr = "Multinom", x = "matrix"),
          definition = function(distr, x) {

  N <- unique(colSums(x))
  if (length(N) != 1) {
    stop("ColSums of x need to be equal. Found multiple values: ", paste(N, " "))
  }

  ncol(x) * lfactorial(distr@size) - sum(lfactorial(x)) +
  sum(t(x) %*% diag(log(distr@prob)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
#' @export
emultinom <- function(x, type = "mle", ...) {

  e(Multinom(), x, type, ...)

}

#' @rdname Multinom
setMethod("mle",
          signature  = c(distr = "Multinom", x = "matrix"),
          definition = function(distr, x) {

  N <- unique(colSums(x))
  if (length(N) != 1) {
    stop("ColSums of x need to be equal. Found multiple values: ", paste(N, " "))
  }

  list(prob = rowMeans(x) / N)

})

#' @rdname Multinom
setMethod("me",
          signature  = c(distr = "Multinom", x = "matrix"),
          definition = function(distr, x) {

  mle(distr, x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multinom
#' @export
vmultinom <- function(size, prob, type = "mle") {

  avar(Multinom(size = size, prob = prob), type = type)

}

#' @rdname Multinom
setMethod("avar_mle",
          signature  = c(distr = "Multinom"),
          definition = function(distr) {

  as.matrix(nearPD(solve(finf(distr))))

})

#' @rdname Multinom
setMethod("avar_me",
          signature  = c(distr = "Multinom"),
          definition = function(distr) {

  avar_mle(distr)

})
