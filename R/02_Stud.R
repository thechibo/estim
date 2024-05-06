# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stud Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Stud",
         contains = "Distribution",
         slots = c(df = "numeric", ncp = "numeric"),
         prototype = list(df = 1, ncp = 0))

#' @title Student Distribution
#' @name Stud
#'
#' @param x an object of class `Stud`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param df,ncp numeric. The distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Stud <- function(df = 1, ncp = 0) {
  new("Stud", df = df, ncp = ncp)
}

setValidity("Stud", function(object) {
  if(length(object@df) != 1) {
    stop("df has to be a numeric of length 1")
  }
  if(length(object@ncp) != 1) {
    stop("ncp has to be a numeric of length 1")
  }
  if(object@df <= 0) {
    stop("df has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Stud
setMethod("d", signature = c(x = "Stud"),
          function(x) {
            function(y, log = FALSE) {
              dt(y, df = x@df, ncp = x@ncp, log = log)
            }
          })

#' @rdname Stud
setMethod("p", signature = c(x = "Stud"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pt(q, df = x@df, ncp = x@ncp,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Stud
setMethod("qn", signature = c(x = "Stud"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qt(p, df = x@df, ncp = x@ncp,
                    lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Stud
setMethod("r", signature = c(x = "Stud"),
          function(x) {
            function(n) {
              rt(n, df = x@df, ncp = x@ncp)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Stud
setMethod("mean",
          signature  = c(x = "Stud"),
          definition = function(x) {

  df <- x@df

  if (df > 1) {
    return(x@ncp * sqrt(df / 2) * gamma((df - 1) / 2) / gamma(df / 2))
  } else {
    stop("Expectation is undefined for Student's t distribution
          no more than 1 df.")
  }

})

#' @rdname Stud
setMethod("median",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@ncp == 0) {
    0
  } else {
    stop("Median not available for non-central Student's t.")
  }

})

#' @rdname Stud
setMethod("mode",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@ncp == 0) {
    0
  } else {
    stop("Mode not available for non-central Student's t.")
  }

})

#' @rdname Stud
setMethod("var",
          signature  = c(x = "Stud"),
          definition = function(x) {

  df <- x@df

  if (df > 2) {
    return(df * (1 + x@ncp ^ 2) / (df - 2) -
             0.5 * x@ncp ^ 2 * df * (gamma((df - 1) / 2) / gamma(df / 2)) ^ 2)
  } else {
    stop("Variance is undefined for Student's t distribution
          no more than 2 df.")
  }

})

#' @rdname Stud
setMethod("sd",
          signature  = c(x = "Stud"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Stud
setMethod("skew",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@ncp == 0) {
    if (x@df > 3) {
      0
    }  else {
        stop("Skewness is undefined for Student's t distribution
              with no more than 3 df.")
      }
  } else {
    stop("Skewness not available for non-central Student's t.")
  }

})

#' @rdname Stud
setMethod("kurt",
          signature  = c(x = "Stud"),
          definition = function(x) {

  if (x@ncp == 0) {
    if (x@df > 4) {
      6 / (x@df - 4)
    } else if (x@df > 2) {
      Inf
    } else {
      stop("Kurtosis is undefined for Student's t distribution
            with no more than 2 df.")
    }
  } else {
    stop("Skewness not available for non-central Student's t.")
  }

})
