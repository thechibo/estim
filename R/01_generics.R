# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generics                                                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#' @param shape1,shape2,size,prob,shape,rate,scale,mean,sd,lambda numeric.
#' Distribution parameters.
#' @param mar numeric. In univariate distributions, a matrix can be given
#' instead of a vector. In this case, the apply function is utilized and `mar`
#' is passed to the `MARGIN` argument. If each sample is a row, set `mar = 1`,
#' if it is a column, set `mar = 2` (the default).
#' @param ... extra arguments.
#'
#' @details
#' The log-likelihood functions are provided in two forms: the `ll<name>`
#' distribution-specific version that follows the base R conventions, and the
#' S4 generic `ll` that complements the `distr` package functionalities.
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
#' @importClassesFrom distr Beta Binom Exp Gammad Norm Pois
#' @importFrom Matrix Matrix nearPD Cholesky
#' @export
#'
#' @seealso [mle], [me], [same]
#'
#' @examples \dontrun{
#' # -------------------------------------------
#' # Beta Distribution Example
#' # -------------------------------------------
#'
#' # Simulation
#'
#' set.seed(1)
#' x <- rbeta(100, shape1, shape2)
#' D <- Beta(shape1, shape2)
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
#' }
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
#' @param par0 function. The estimator to use for initialization of the
#' likelihood maximization algorithm.
#' @param method,lower,upper arguments passed to optim.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [me], [same]
#'
#' @inherit estim examples
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
#' @param dirich logical. Should the Dirichlet-based estimator be calculated
#' instead? Applies only to the Multivariate Gamma distribution.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [mle], [same]
#'
#' @inherit estim examples
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
#' @param dirich logical. Should the Dirichlet-based estimator be calculated
#' instead? Applies only to the Multivariate Gamma distribution.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#' @export
#'
#' @seealso [estim], [mle], [me]
#'
#' @inherit estim examples
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
#' @param shape1,shape2,size,prob,shape,rate,scale,mean,sd,lambda numeric.
#' Distribution parameters.
#' @param ... extra arguments.
#'
#' @return A named matrix. The asymptotic covariance matrix of the estimator.
#'
#' @export
#'
#' @seealso [avar_mle], [avar_me], [avar_same]
#'
#' @inherit estim examples
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
#' @param dirich logical. Should the Dirichlet-based estimator be calculated
#' instead? Applies only to the Multivariate Gamma distribution.
#' @param comp logical. Should the component matrices A, B be returned instead?
#' This argument is used internally by the `avar_me` method for the
#' Multivariate Gamma Dirichlet-based ME that needs these matrices.
#'
#' @export
#'
#' @seealso [avar], [avar_mle], [avar_same]
#'
#' @inherit avar params return examples
setGeneric("avar_me", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_me") })

#' @title SAME Asymptotic Variance
#'
#' @description
#' Calculates the asymptotic variance (or variance - covariance matrix in the
#' multidimensional case) of the SAME, given a specified family of
#' distributions and the true parameter values.
#'
#' @param dirich logical. Should the Dirichlet-based estimator be calculated
#' instead? Applies only to the Multivariate Gamma distribution.
#' @param comp logical. Should the component matrices A, B be returned instead?
#' This argument is used internally by the `avar_same` method for the
#' Multivariate Gamma Dirichlet-based SAME that needs these matrices.
#'
#' @export
#'
#' @seealso [avar], [avar_mle], [avar_me]
#'
#' @inherit avar params return examples
setGeneric("avar_same", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_same") })
