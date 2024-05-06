# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generics & Classes                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Distribution S4 Classes
#' @name Distributions
#'
#' @description
#' A collection of classes that provide a flexible and structured way to work
#' with probability distributions.
#'
#' @return
#' The dpqr family of functions return the evaluated density, cumulative
#' probability, quantile, and random sample, respectively.
#' The moments family of functions return the appropriate theoretical moment,
#' as calculated by the distribution true parameters.
#' The ll function returns the evaluated log-likelihood, given a sample and the
#' theoretical parameters.
#' The estim family of functions return the estimated parameters of the
#' distribution, given a sample.
#' The avar family of functions return the asymptotic variance or variance -
#' covariance matrix (if there are two or more parameters) of the corresponding
#' estimation method.
#' Calculus performed on Distribution objects returns a Distribution object of
#' the appropriate class and with the appropriate parameters.
#'
NULL

setClass("Distribution")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title The d p q r Functions
#' @name dpqr
#' @aliases d p q r
#'
#' @description
#' Four generic functions that take a distribution object (e.g. `Bern`) and
#' return the density, cumulative probability, quantile, and random generator
#' functions, respectively.
#'
#' @param x an object of subclass `Distribution`.
#' @param ... extra arguments.
#'
#' @return The d p q r functions return the density, cumulative probability,
#' quantile, and random generator functions, respectively.
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' library(estimators)
#'
#' # Create the distribution
#' x <- Beta(3, 5)
#'
#' # Density function
#' df <- d(x)
#' df(c(0.3, 0.8, 0.5))
#'
#' # Probability function
#' pf <- p(x)
#' pf(c(0.3, 0.8, 0.5))
#'
#' # Density function
#' qf <- qn(x)
#' qf(c(0.3, 0.8, 0.5))
#'
#' # Random Generator function
#' rf <- r(x)
#' rf(5)
setGeneric("d", function(x, ...) {
  standardGeneric("d")
})

#' @rdname dpqr
#' @export
setGeneric("p", function(x, ...) {
  standardGeneric("p")
})

#' @rdname dpqr
#' @export
setGeneric("qn", function(x, ...)
  standardGeneric("qn"))

#' @rdname dpqr
#' @export
setGeneric("r", function(x, ...) {
  standardGeneric("r")
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Moments - Parametric Quantities of Interest
#' @name moments
#'
#' @description A set of functions that calculate the theoretical moments
#' (expectation, variance, skewness, excess kurtosis) and other important
#' parametric functions (median, mode, entropy, Fisher information) of a
#' distribution.
#'
#' @param x an object of a `Distribution` subclass.
#' @param y,use,na.rm arguments in `mean` and `var` standard methods from the
#' `stats` package not used here.
#' @param ... extra arguments.
#'
#' @details
#' The `moments()` function automatically finds the available methods for a
#' given distribution and results all of the results in a list.
#'
#' Not all functions are available for distributions; for example, the `sd()`
#' is available only for univariate distributions.
#'
#' @return Numeric, either vector or matrix depending on the moment and the
#' distribution. Function `moments()` returns a list of all available methods.
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' library(estimators)
#'
#' # Create the distribution
#' x <- Beta(3, 5)
#'
#' # List of all available moments
#' mom <- moments(x)
#'
#' # Expectation
#' mean(x)
#' mom$mean
#'
#' # Variance and Standard Deviation
#' var(x)
#' sd(x)
#'
#' # Skewness and Excess Kurtosis
#' skew(x)
#' kurt(x)
#'
#' # Entropy
#' entro(x)
#'
#' # Fisher Information Matrix
#' finf(x)
moments <- function(x) {
  mom <- get_moment_methods(x)
  y <- lapply(mom, FUN = function(m) { do.call(m, list(x = x)) })
  names(y) <- mom
  y
}

#' @rdname moments
#' @name mean
#' @usage mean(x, ...)
NULL

#' @rdname moments
setGeneric("median")

#' @rdname moments
setGeneric("mode")

#' @rdname moments
setGeneric("var")

#' @rdname moments
setGeneric("sd")

#' @rdname moments
setGeneric("skew", function(x, ...) {
  standardGeneric("skew")
})

#' @rdname moments
setGeneric("kurt", function(x, ...) {
  standardGeneric("kurt")
})

#' @rdname moments
setGeneric("entro", function(x, ...) {
  standardGeneric("entro")
})

#' @rdname moments
setGeneric("finf", function(x, ...) {
  standardGeneric("finf")
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Log-Likelihood
#'
#' @description
#' These functions calculate the log-likelihood of an IID sample for specific
#' values of the distribution parameters. See Details.
#'
#' @param x numeric. A sample under estimation.
#' @param prm numeric. A vector of the distribution parameters.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param location,alpha,mu,sigma,meanlog,sdlog,min,max,size,prob,shape,rate,scale,mean,sd,lambda numeric.
#' Distribution parameters.
#' @param ... extra arguments.
#'
#' @details
#' The log-likelihood functions are provided in two forms: the `ll<name>`
#' distribution-specific version that follows the base R conventions, and the
#' S4 generic `ll`.
#'
#' @return Numeric. The value of the log-likelihood function.
#' @export
setGeneric("ll", signature = c("x", "prm", "distr"),
           function(x, prm, distr, ...) { standardGeneric("ll") })

#' @rdname ll
setMethod("ll",
          signature  = c(x = "ANY", "prm" = "ANY", distr = "character"),
          definition = function(x, prm, distr, ...) {

            ll(x, prm, get_distr_class(distr), ...)

          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setGeneric("lloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("lloptim") })

setGeneric("dlloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("dlloptim") })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Parameter Estimation
#'
#' @description
#' Estimates the parameters of a random sample according to a
#' specified family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#'
#' @importFrom Matrix Matrix nearPD Cholesky
#' @export
#'
#' @references
#' Ye, Z.-S. & Chen, N. (2017), Closed-form estimators for the gamma
#' distribution derived from likelihood equations, The American Statistician
#' 71(2), 177–181.
#'
#' Van der Vaart, A. W. (2000), Asymptotic statistics, Vol. 3,
#' Cambridge university press.
#'
#' Tamae, H., Irie, K. & Kubokawa, T. (2020), A score-adjusted approach to
#' closed-form estimators for the gamma and beta distributions, Japanese Journal
#' of Statistics and Data Science 3, 543–561.
#'
#' Mathal, A. & Moschopoulos, P. (1992), A form of multivariate gamma
#' distribution, Annals of the Institute of Statistical Mathematics 44, 97–106.
#'
#' Oikonomidis, I. & Trevezas, S. (2023), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @seealso [mle], [me], [same]
#'
#' @examples
#' # -----------------------------------------------------
#' # Beta Distribution Example
#' # -----------------------------------------------------
#'
#' # Simulation
#' set.seed(1)
#' shape1 <- 1
#' shape2 <- 2
#' D <- Beta(shape1, shape2)
#' x <- r(D)(100)
#'
#' # Likelihood - The ll Functions
#'
#' llbeta(x, shape1, shape2)
#' ll(x, c(shape1, shape2), D)
#' ll(x, c(shape1, shape2), "beta")
#'
#' # Point Estimation - The e Functions
#'
#' ebeta(x, type = "mle")
#' ebeta(x, type = "me")
#' ebeta(x, type = "same")
#'
#' mle(x, D)
#' me(x, D)
#' same(x, D)
#'
#' estim(x, D, type = "mle")
#'
#' # Asymptotic Variance - The v Functions
#'
#' vbeta(shape1, shape2, type = "mle")
#' vbeta(shape1, shape2, type = "me")
#' vbeta(shape1, shape2, type = "same")
#'
#' avar_mle(D)
#' avar_me(D)
#' avar_same(D)
#'
#' avar(D, type = "mle")
estim <- function(x, distr, type = "mle", ...) {
  type <- tolower(type)
  if (type %in% c("mle", "me", "same")) {
    return(do.call(type, list(x = x, distr = distr, ...)))
  } else {
    stop("Type must be one of mle, me, or same, case ignored. Instead got",
         type)
  }
}

#' @title Maximum Likelihood Estimation
#'
#' @description
#' Calculates the MLE under the assumption the sample observations are
#' independent and identically distributed (iid) according to a
#' specified family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [me], [same]
#'
#' @inherit estim references examples
setGeneric("mle", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("mle") })

#' @rdname mle
setMethod("mle",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            mle(x, get_distr_class(distr), ...)

          })

#' @title Moment Estimation
#'
#' @description
#' Calculates the ME under the assumption the sample observations are
#' independent and identically distributed (iid) according to a
#' specified family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [mle], [same]
#'
#' @inherit estim references examples
setGeneric("me", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("me") })

#' @rdname me
setMethod("me",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            me(x, get_distr_class(distr), ...)

          })

#' @title Score - Adjusted Moment Estimation
#'
#' @description
#' Calculates the SAME under the assumption the sample observations are
#' independent and identically distributed (iid) according to a
#' specified family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [mle], [me]
#'
#' @inherit estim references examples
setGeneric("same", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("same") })

#' @rdname same
setMethod("same",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            same(x, get_distr_class(distr), ...)

          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Avar                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Asymptotic Variance
#'
#' @description
#' Calculates the asymptotic variance (or variance - covariance matrix in the
#' multidimensional case) of an estimator, given a specified family of
#' distributions and the true parameter values.
#'
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param alpha,mu,sigma,size,prob,shape,rate,scale,mean,sd,lambda numeric.
#' Distribution parameters.
#' @param ... extra arguments.
#'
#' @return A named matrix. The asymptotic covariance matrix of the estimator.
#'
#' @export
#'
#' @seealso [avar_mle], [avar_me], [avar_same]
#'
#' @inherit estim references examples
avar <- function(distr, type, ...) {
  type <- tolower(type)
  if (type %in% c("mle", "me", "same")) {
    return(do.call(paste0("avar_", type), list(distr = distr, ...)))
  } else {
    stop("Method must be one of mle, me, or same, case ignored. Instead got",
         type)
  }
}

#' @title MLE Asymptotic Variance
#'
#' @description
#' Calculates the asymptotic variance (or variance - covariance matrix in the
#' multidimensional case) of the MLE, given a specified family of
#' distributions and the true parameter values.
#'
#' @export
#'
#' @seealso [avar], [avar_me], [avar_same]
#'
#' @inherit avar params return examples
setGeneric("avar_mle", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_mle") })

#' @title ME Asymptotic Variance
#'
#' @description
#' Calculates the asymptotic variance (or variance - covariance matrix in the
#' multidimensional case) of the ME, given a specified family of
#' distributions and the true parameter values.
#'
#' @export
#'
#' @seealso [avar], [avar_mle], [avar_same]
#'
#' @inherit avar params return references examples
setGeneric("avar_me", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_me") })

#' @title SAME Asymptotic Variance
#'
#' @description
#' Calculates the asymptotic variance (or variance - covariance matrix in the
#' multidimensional case) of the SAME, given a specified family of
#' distributions and the true parameter values.
#'
#' @export
#'
#' @seealso [avar], [avar_mle], [avar_me]
#'
#' @inherit avar params return references examples
setGeneric("avar_same", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_same") })
