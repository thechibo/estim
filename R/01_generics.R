# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generics                                                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Log-Likelihood
#'
#' @param prm numeric. A vector of the distribution parameters.
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution of interest.
#' @param ... extra arguments.
#'
#' @return Numeric. The value of the log-likelihood function.
#' @export
setGeneric("ll", signature = c("prm", "x", "distr"),
           function(prm, x, distr, ...) { standardGeneric("ll") })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setGeneric("lloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("lloptim") })

setGeneric("dlloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("dlloptim") })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MLE                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Maximum Likelihood Estimator
#'
#' @description
#' Calculates the MLE of a sample under the assumption the observations are
#' independent and identically distributed (iid) according to a specified
#' family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param par0 function. The estimator to use for initialization of the
#' likelihood maximization algorithm.
#' @param method,lower,upper arguments passed to optim.
#' @param ... extra arguments.
#'
#' @return numeric. The estimator produced by the sample.
#'
#' @importClassesFrom distr Beta Gammad
#' @importFrom Matrix Matrix nearPD Cholesky
#' @export
#'
#' @examples \dontrun{
#' # Distribution
#' D <- Beta(shape1 = 1, shape2 = 1.5)
#'
#' # Simulation
#' set.seed(2)
#' x <- r(D)(50)
#'
#' # Estimation
#' me(x, "Beta")
#' same(x, "Beta")
#' mle(x, "Beta")
#'
#' # Asymptotic Covariance Matrix
#' acov_me(D)
#' acov_same(D)
#' acov_mle(D)
#' }
setGeneric("mle", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("mle") })

#' Asymptotic covariance matrix of estimator
#'
#' @param distr A subclass of `Distribution`. The distribution of interest.
#' @param ... extra arguments.
#'
#' @return A named matrix. The asymptotic covariance matrix of the estimator.
#'
#' @export
#'
#'@inherit mle examples
setGeneric("acov_mle", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_mle") })

#' @rdname mle
setMethod("mle",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            mle(x, get_distr_class(distr), ...)

          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Moment Estimator
#'
#' @description
#' Calculates the MΕ of a sample under the assumption the observations are
#' independent and identically distributed (iid) according to a specified
#' family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param ... extra arguments.
#'
#' @inherit mle return examples
#' @export
setGeneric("me", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("me") })

#' @inherit acov_mle title params return
#' @inherit mle examples
#'
#' @param comp logical. Should the components A, B of the delta method be
#' returned instead of the final variance-covariance matrix?
#'
#' @export
setGeneric("acov_me", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_me") })

#' @rdname me
setMethod("me",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            me(x, get_distr_class(distr), ...)

          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Score-Adjusted Moment Estimator
#'
#' @description
#' Calculates the SAMΕ of a sample under the assumption the observations are
#' independent and identically distributed (iid) according to a specified
#' family of distributions.
#'
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution family assumed.
#' @param ... extra arguments.
#'
#' @inherit mle return examples
#' @export
setGeneric("same", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("same") })

#' @inherit acov_mle title params return
#' @inherit mle examples
#'
#' @param comp logical. Should the components A, B of the delta method be
#' returned instead of the final variance-covariance matrix?
#'
#' @export
setGeneric("acov_same", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_same") })

#' @rdname same
setMethod("same",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

            same(x, get_distr_class(distr), ...)

          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME2                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setGeneric("me2", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("me2") })

#' @rdname acov_me
setGeneric("acov_me2", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_me2") })

#' @rdname me
setMethod("me2",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

  me2(x, get_distr_class(distr), ...)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME2                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setGeneric("same2", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("same2") })

#' @rdname acov_same
setGeneric("acov_same2", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_same2") })

#' @rdname same
setMethod("same2",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

  same2(x, get_distr_class(distr), ...)

})
