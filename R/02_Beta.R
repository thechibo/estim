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
#' @description
#' The Beta distribution is a absolute continuous probability distribution with
#' support \eqn{S = [0,1]}, parameterized by two shape parameters,
#' \eqn{\alpha > 0} and \eqn{\beta > 0}.
#'
#' @param n numeric. The sample size.
#' @param distr,x If both arguments coexist, `distr` is an object of class
#' `Beta` and `x` is a numeric vector, the sample of observations. For the
#' moment functions that only take an `x` argument, `x` is an object of class
#' `Beta` instead.
#' @param shape1,shape2 numeric. The distribution parameters (positive real
#' numbers).
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization.
#'
#' @details
#' The probability density function (PDF) of the Beta distribution is given by:
#' \deqn{ f(x; \alpha, \beta) = \frac{x^{\alpha - 1} (1 - x)^{\beta - 1}}{B(\alpha, \beta)},
#' \quad \alpha\in\mathbb{R}_+, \, \beta\in\mathbb{R}_+}
#' for \eqn{x \in S = [0, 1]}, where \eqn{B(\alpha, \beta)} is the Beta function:
#' \deqn{ B(\alpha, \beta) = \int_0^1 t^{\alpha - 1} (1 - t)^{\beta - 1} dt.}
#'
#' @section Moments:
#' The Beta distribution has the following known moments:
#' - **Mean**: \eqn{E[X] = \frac{\alpha}{\alpha + \beta}.}
#' - **Median**: Approximated by \eqn{\frac{\alpha - 1/3}{\alpha + \beta - 2/3}}, for \eqn{\alpha, \beta > 1}
#' - **Mode**: \eqn{\frac{\alpha - 1}{\alpha + \beta - 2}}, for \eqn{\alpha, \beta > 1.}
#' - **Variance**: \eqn{Var(X) = \frac{\alpha\beta}{(\alpha + \beta)^2 (\alpha + \beta + 1)}.}
#' - **Skewness**: \eqn{\frac{2(\beta - \alpha)\sqrt{\alpha + \beta + 1}}{(\alpha + \beta + 2)\sqrt{\alpha\beta}}}
#' - **Kurtosis**: \eqn{\frac{6[(\alpha - \beta)^2 (\alpha + \beta + 1) - \alpha\beta (\alpha + \beta + 2)]}{\alpha\beta (\alpha + \beta + 2) (\alpha + \beta + 3)}}
#' - **Entropy**: \eqn{H(X) = \log B(\alpha, \beta) - (\alpha - 1)\psi(\alpha) - (\beta - 1)\psi(\beta) + (\alpha + \beta - 2)\psi(\alpha + \beta)}
#' - **Fisher Information**: The Fisher information matrix for \eqn{(\alpha, \beta)} is given by:
#'   \deqn{ I(\alpha, \beta) = \begin{bmatrix} \psi' (\alpha) - \psi' (\alpha + \beta) & -\psi' (\alpha + \beta) \\ -\psi' (\alpha + \beta) & \psi' (\beta) - \psi' (\alpha + \beta) \end{bmatrix} }
#'
#' @section Estimation:
#' The parameters \eqn{(\alpha, \beta)} can be estimated using:
#' - **Method of Moments**: Solving for \eqn{\alpha} and \eqn{\beta} using sample mean and variance.
#' - **Maximum Likelihood Estimation (MLE)**: Numerical methods are typically used to maximize the log-likelihood function.
#'
#' @inherit Distributions return
#'
#' @references
#'
#' - Tamae, H., Irie, K. & Kubokawa, T. (2020), A score-adjusted approach to
#' closed-form estimators for the gamma and beta distributions, Japanese Journal
#' of Statistics and Data Science 3, 543–561.
#'
#' - Papadatos, N. (2022), On point estimators for gamma and beta distributions,
#' arXiv preprint arXiv:2205.10799.
#'
#' - Oikonomidis, I. & Trevezas, S. (2025), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @seealso
#' Functions from the `stats` package: [dbeta()], [pbeta()], [qbeta()], [rbeta()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- 3 ; b <- 5
#' D <- Beta(a, b)
#' x <- c(0.3, 0.8, 0.5)
#' n <- 5
#'
#' # Density function
#' df <- d(D) ; df(x) # df is a function itself
#' d(D, x) # alternative way to use the function
#'
#' # Distribution function
#' pf <- p(D) ; pf(x)
#' p(D, x)
#'
#' # Inverse distribution function
#' qf <- qn(D) ; qf(x)
#' qn(D, x)
#'
#' # Random Generator function
#' rf <- r(D) ; rf(n)
#' r(D, n)
#'
#' # List of all available moments
#' mom <- moments(D)
#'
#' # Expectation
#' mean(D)
#' mom$mean
#'
#' # Variance and Standard Deviation
#' var(D)
#' sd(D)
#'
#' # Skewness and Excess Kurtosis
#' skew(D)
#' kurt(D)
#'
#' # Entropy
#' entro(D)
#'
#' # Fisher Information Matrix
#' finf(D)
#'
#' # Point Estimation - The e Functions
#'
#' ebeta(x, type = "mle")
#' ebeta(x, type = "me")
#' ebeta(x, type = "same")
#'
#' mle(D, x)
#' me(D, x)
#' same(D, x)
#' e(D, x, type = "mle")
#'
#' mle("beta", x) # the distr argument can be a character
#'
#' # Asymptotic Variance - The v Functions
#'
#' vbeta(a, b, type = "mle")
#' vbeta(a, b, type = "me")
#' vbeta(a, b, type = "same")
#'
#' avar_mle(D)
#' avar_me(D)
#' avar_same(D)
#'
#' avar(D, type = "mle")
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
setMethod("d", signature = c(distr = "Beta", x = "numeric"),
          function(distr, x) {
            dbeta(x, shape1 = distr@shape1, shape2 = distr@shape2)
          })

#' @rdname Beta
#' @export
setMethod("p", signature = c(distr = "Beta", x = "numeric"),
          function(distr, x) {
            pbeta(x, shape1 = distr@shape1, shape2 = distr@shape2)
          })

#' @rdname Beta
#' @export
setMethod("qn", signature = c(distr = "Beta", x = "numeric"),
          function(distr, x) {
            qbeta(x, shape1 = distr@shape1, shape2 = distr@shape2)
          })

#' @rdname Beta
#' @export
setMethod("r", signature = c(distr = "Beta", n = "numeric"),
          function(distr, n) {
            rbeta(n, shape1 = distr@shape1, shape2 = distr@shape2)
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
    warning("In Beta(1, 1), all elements in the [0, 1] interval are modes. 0.5 is returned by default.")
    return(0.5)
  } else if (a < 1 && b < 1) {
    warning("In Beta(a, b) with a < 1 and b < 1, both 0 and 1 are modes. 1 is returned by default.")
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
  ll(distr = Beta(shape1, shape2), x = x)
}

#' @rdname Beta
#' @export
setMethod("ll",
          signature  = c(distr = "Beta", x = "numeric"),
          definition = function(distr, x) {

  a <- distr@shape1
  b <- distr@shape2
  a0 <- a + b

  length(x) * (lgamma(a0) - lgamma(a) - lgamma(b)) +
   (a - 1) * sum(log(x)) + (b - 1) * sum(log(1 - x))

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

  # Shape parameters (a, b) as a function of a0
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

#' @rdname Beta
#' @export
ebeta <- function(x, type = "mle", ...) {

  e(Beta(), x = x, type = type, ...)

}

#' @rdname Beta
#' @export
setMethod("mle",
          signature  = c(distr = "Beta", x = "numeric"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- c(mean(log(x)), mean(log(1 - x)))

  par <- optim(par = sum(unlist(do.call(par0, list(distr = distr, x = x)))),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  shape <- idigamma(digamma(par) + tx)

  list(shape1 = shape[1], shape2 = shape[2])

})

#' @rdname Beta
#' @export
setMethod("me",
          signature  = c(distr = "Beta", x = "numeric"),
          definition = function(distr, x) {

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  d  <- (m - m2) / (m2 - m ^ 2)

  list(shape1 = d * m, shape2 = d * (1 - m))

})

#' @rdname Beta
#' @export
setMethod("same",
          signature  = c(distr = "Beta", x = "numeric"),
          definition = function(distr, x) {

  mx <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  my <- 1 - mx
  mly <- mean(log(1 - x))
  myly <- mean((1 - x) * log(1 - x))
  s <- mxlx - mx * mlx + myly - my * mly

  list(shape1 = mx / s, shape2 = my / s)

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
