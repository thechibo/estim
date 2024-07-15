# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Beta",
  contains = "Distribution",
  slots = c(shape1 = "numeric", shape2 = "numeric"),
  prototype = list(shape1 = 1, shape2 = 1))

#' @title Beta Distribution
#' @name Beta
#'
#' @param x an object of class `Beta`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Beta`.
#' @param shape1,shape2 numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param par0,method,lower,upper arguments passed to optim.
#'
#' @inherit Distributions return
#'
#' @export
Beta <- function(shape1 = 1, shape2 = 1) {
  new("Beta", shape1 = shape1, shape2 = shape2)
}

setValidity("Beta", function(object) {
  if(length(object@shape1) != 1) {
    stop("shape1 has to be a numeric of length 1")
  }
  if(object@shape1 <= 0) {
    stop("shape1 has to be positive")
  }
  if(length(object@shape2) != 1) {
    stop("shape2 has to be a numeric of length 1")
  }
  if(object@shape2 <= 0) {
    stop("shape2 has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Beta
#' @export
setMethod("d", signature = c(x = "Beta"),
          function(x) {
            function(y, log = FALSE) {
              dbeta(y, shape1 = x@shape1, shape2 = x@shape2, ncp = 0,
                    log = log)
            }
          })

#' @rdname Beta
#' @export
setMethod("p", signature = c(x = "Beta"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pbeta(q, shape1 = x@shape1, shape2 = x@shape2, ncp = 0,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Beta
#' @export
setMethod("qn", signature = c(x = "Beta"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qbeta(p, shape1 = x@shape1, shape2 = x@shape2, ncp = 0,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Beta
#' @export
setMethod("r", signature = c(x = "Beta"),
          function(x) {
            function(n) {
              rbeta(n, shape1 = x@shape1, shape2 = x@shape2, ncp = 0)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Beta
#' @export
setMethod("mean",
          signature  = c(x = "Beta"),
          definition = function(x) {

  x@shape1 / (x@shape1 + x@shape2)

})

#' @rdname Beta
#' @export
setMethod("median",
          signature  = c(x = "Beta"),
          definition = function(x) {

  qbeta(0.5, shape1 = x@shape1, shape2 = x@shape2)

})

#' @rdname Beta
#' @export
setMethod("mode",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  if (a > 1 && b > 1) {
    return((a - 1) / (a + b - 2))
  } else if (a == 1 && b == 1) {
    warning("In Beta(1, 1), all elements in the [0, 1] interval are modes.
             0.5 is returned by default.")
    return(0.5)
  } else if (a < 1 && b < 1) {
    warning("In Beta(a, b) with a < 1 and b < 1, both 0 and 1 are modes.
             1 is returned by default.")
    return(1)
  } else if (a <= 1) {
    return(0)
  } else {
    return(1)
  }

})

#' @rdname Beta
#' @export
setMethod("var",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  (a * b) / ((a + b) ^ 2 * (a + b + 1))

})

#' @rdname Beta
#' @export
setMethod("sd",
          signature  = c(x = "Beta"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Beta
#' @export
setMethod("skew",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  (2 * (b - a) * sqrt(a + b + 1)) / ((a + b + 2) * sqrt(a * b))

})

#' @rdname Beta
#' @export
setMethod("kurt",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  (6 * (a - b) ^ 2 * (a + b + 1) - a * b * (a + b + 2)) /
    (a * b * (a + b + 2) * (a + b + 3))

})

#' @rdname Beta
#' @export
setMethod("entro",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  lbeta(a, b) - (a - 1) * digamma(a) - (b - 1) * digamma(b) +
  (a + b - 2) * digamma(a + b)

})

#' @rdname Beta
#' @export
setMethod("finf",
          signature  = c(x = "Beta"),
          definition = function(x) {

  a <- x@shape1
  b <- x@shape2

  p1a  <- trigamma(a)
  p1b  <- trigamma(b)
  p1   <- trigamma(a + b)

  D <- matrix(c(p1a - p1, - p1, - p1, p1b - p1), nrow = 2, ncol = 2)

  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Beta
#' @export
llbeta <- function(x, shape1, shape2) {
  ll(x, prm = c(shape1, shape2), distr = Beta())
}

#' @rdname Beta
#' @export
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Beta"),
          definition = function(x, prm, distr) {

  length(x) * (lgamma(sum(prm)) - lgamma(prm[1]) - lgamma(prm[2])) +
   (prm[1] - 1) * sum(log(x)) + (prm[2] - 1) * sum(log(1 - x))

})

# Bias Corrected log-likelihood
# (Firth, 1993, Cribari-Neto and Vasconcellos, 2010)
#ll = function(prm, x) {
#  p1a = trigamma(prm[1])
#  p1b = trigamma(prm[2])
#  p1  = trigamma(sum(prm))
#  d   = p1a * p1b - p1 * (p1a + p1b)
#  ld = log((length(x) ^ 2) * d) / 2
#  sum(do.call(dbeta, c(list(x = x, log = TRUE), prm))) + ld
#}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Beta"),
          definition = function(par, tx, distr) {

  a <- idigamma(digamma(par) + tx)
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * tx)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Beta"),
          definition = function(par, tx, distr) {

  # Shape parameters (a, b) as a function of a0
  a <- idigamma(digamma(par) + tx)

  # a_i derivative wrt a0
  da <- trigamma(par) / trigamma(a)

  # lloptim derivative wrt a0 (par)
  digamma(sum(a)) * sum(da) - sum(digamma(a) * da) + sum(tx * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname estimation
#' @export
ebeta <- function(x, type = "mle", ...) {

  estim(x, Beta(), type, ...)

}

#' @rdname Beta
#' @export
setMethod("mle",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- c(mean(log(x)), mean(log(1 - x)))

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  shape <- idigamma(digamma(par) + tx)

  names(shape) <- paste0("shape", seq_along(shape))
  shape

})

#' @rdname Beta
#' @export
setMethod("me",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  d  <- (m - m2) / (m2 - m ^ 2)

  c(shape1 = d * m, shape2 = d * (1 - m))

})

#' @rdname Beta
#' @export
setMethod("same",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  mx <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  my <- 1 - mx
  mly <- mean(log(1 - x))
  myly <- mean((1 - x) * log(1 - x))
  s <- mxlx - mx * mlx + myly - my * mly

  c(shape1 = mx / s, shape2 = my / s)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Beta
#' @export
vbeta <- function(shape1, shape2, type = "mle") {

  avar(Beta(shape1 = shape1, shape2 = shape2), type = type)

}

#' @rdname Beta
#' @export
setMethod("avar_mle",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Beta
#' @export
setMethod("avar_me",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr@shape1
  b <- distr@shape2

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  s4  <- s2 ^ 2
  m3  <- 2 * (b - a) * s2 / (th * (th + 2))
  m4  <- 3 * prd * (prd * (th + 2) + 2 * (b - a) ^2) /
    ((th ^ 4) * (th + 1) * (th + 2) * (th + 3))
  d   <- (th + 1) ^ 2 * (th + 2) ^ 2 * s2
  e   <- (th + 1) ^ 3 * (m4 - s4 - m3 ^ 2 / s2) / s2

  s11 <- (a * (a + 1)) ^ 2 / d + a * e / b
  s22 <- (b * (b + 1)) ^ 2 / d + b * e / a
  s12 <- - a * (a + 1) * b * (b + 1) / d + e

  D <- matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)
  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})

#' @rdname Beta
#' @export
setMethod("avar_same",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr@shape1
  b <- distr@shape2

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  p1a <- trigamma(a)
  p1b <- trigamma(b)
  m1  <- matrix(c(a ^ 2, prd, prd, b ^ 2), nrow = 2, ncol = 2)
  m2  <- matrix(c(prd, th2 - prd, th2 - prd, prd), nrow = 2, ncol = 2)

  D <- (s2 * th2 * (p1a + p1b) + 1) * m1 - m2 / (th + 1)
  rownames(D) <- c("shape1", "shape2")
  colnames(D) <- c("shape1", "shape2")

  D

})
