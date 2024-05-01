# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fisher Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Fisher",
         contains = "Distribution",
         slots = c(df = "numeric", ncp = "numeric"),
         prototype = list(df = 1, ncp = 0))

#' @title Fisher Distribution
#' @name Fisher
#'
#' @param x an object of class `Fisher`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param df1,df2,ncp numeric. The distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Fisher <- function(df1 = 1, df2 = 1, ncp = 0) {
  new("Fisher", df1 = df1, df2 = df2, ncp = ncp)
}

setValidity("Fisher", function(object) {
  if(length(object@df1) != 1) {
    stop("df1 has to be a numeric of length 1")
  }
  if(length(object@df2) != 1) {
    stop("df2 has to be a numeric of length 1")
  }
  if(length(object@ncp) != 1) {
    stop("ncp has to be a numeric of length 1")
  }
  if(object@df1 <= 0) {
    stop("df1 has to be positive")
  }
  if(object@df2 <= 0) {
    stop("df2 has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Fisher
setMethod("d", signature = c(x = "Fisher"),
          function(x) {
            function(y, log = FALSE) {
              df(y, df1 = x@df1, df2 = x@df2, ncp = x@ncp, log = log)
            }
          })

#' @rdname Fisher
setMethod("p", signature = c(x = "Fisher"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pf(q, df1 = x@df1, df2 = x@df2, ncp = x@ncp,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Fisher
setMethod("q2", signature = c(x = "Fisher"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qf(p, df1 = x@df1, df2 = x@df2, ncp = x@ncp,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Fisher
setMethod("r", signature = c(x = "Fisher"),
          function(x) {
            function(n) {
              rf(n, df1 = x@df1, df2 = x@df2, ncp = x@ncp)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Fisher
setMethod("mean",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  if (x@df2 > 2) {
    return(x@df2 * (x@df1 + x@ncp) / (x@df1 * (x@df2 - 2)))
  } else {
    stop("Expectation is undefined for F distribution with df2 <= 2.")
  }

})

#' @rdname Fisher
setMethod("mode",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  if (x@ncp == 0) {
    if (x@df1 > 2) {
      return(x@df1 - 2) * x@df2 / (x@df1 * (x@df2 + 2))
    } else {
      stop("Expectation is undefined for F distribution with df1 <= 2.")
    }

  } else {
    stop("Skewness not available for non-central F.")
  }

})

#' @rdname Fisher
setMethod("var",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (x@df2 > 4) {
    return(2 * ((n1 + x@ncp) ^ 2 + (n1 + 2 * x@ncp) * (n2 - 2)) *
             (n2 / n1) ^ 2 / ((n2 - 2) ^ 2 * (n2 - 4)))
  } else {
    stop("Variance is undefined for Fdistribution with df2 <= 4.")
  }

})

#' @rdname Fisher
setMethod("sd",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Fisher
setMethod("skew",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (x@ncp == 0) {
    if (n2 > 6) {
      ((2 * n1 + n2 - 2) * sqrt(8 * (n2 - 4))) /
        ((n2 - 6) * sqrt(n1 * (n1 + n2 - 2)))
    }  else {
      stop("Skewness is undefined for F distribution with df2 <= 6.")
    }
  } else {
    stop("Skewness not available for non-central F.")
  }

})
