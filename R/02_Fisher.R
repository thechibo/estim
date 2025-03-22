# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fisher Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Fisher",
         contains = "Distribution",
         slots = c(df1 = "numeric", df2 = "numeric"),
         prototype = list(df1 = 1, df2 = 1))

#' @title Fisher Distribution
#' @name Fisher
#'
#' @param x an object of class `Fisher`. If the function also has a `distr`
#' argument, `x` is a numeric vector, a sample of observations.
#' @param distr an object of class `Fisher`.
#' @param df1,df2 numeric. The distribution parameters.
#' @param prm numeric. A vector including the distribution parameters.
#'
#' @inherit Distributions return
#'
#' @export
Fisher <- function(df1 = 1, df2 = 1) {
  new("Fisher", df1 = df1, df2 = df2)
}

setValidity("Fisher", function(object) {
  if(length(object@df1) != 1) {
    stop("df1 has to be a numeric of length 1")
  }
  if(length(object@df2) != 1) {
    stop("df2 has to be a numeric of length 1")
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
              df(y, df1 = x@df1, df2 = x@df2, ncp = 0, log = log)
            }
          })

#' @rdname Fisher
setMethod("p", signature = c(x = "Fisher"),
          function(x) {
            function(q, lower.tail = TRUE, log.p = FALSE) {
              pf(q, df1 = x@df1, df2 = x@df2, ncp = 0,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Fisher
setMethod("qn", signature = c(x = "Fisher"),
          function(x) {
            function(p, lower.tail = TRUE, log.p = FALSE) {
              qf(p, df1 = x@df1, df2 = x@df2, ncp = 0,
                 lower.tail = lower.tail, log.p = log.p)
            }
          })

#' @rdname Fisher
setMethod("r", signature = c(x = "Fisher"),
          function(x) {
            function(n) {
              rf(n, df1 = x@df1, df2 = x@df2, ncp = 0)
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
    return(x@df2 * x@df1 / (x@df1 * (x@df2 - 2)))
  } else {
    stop("Expectation is undefined for F distribution with df2 <= 2.")
  }

})

#' @rdname Fisher
setMethod("median",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  qf(0.5, df1 = x@df1, df2 = x@df2)

})

#' @rdname Fisher
setMethod("mode",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  if (x@df1 > 2) {
    return(x@df1 - 2) * x@df2 / (x@df1 * (x@df2 + 2))
  } else {
    stop("Expectation is undefined for F distribution with df1 <= 2.")
  }

})

#' @rdname Fisher
setMethod("var",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (x@df2 > 4) {
    return(2 * (n1 ^ 2 + n1 * (n2 - 2)) *
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

  if (n2 > 6) {
    ((2 * n1 + n2 - 2) * sqrt(8 * (n2 - 4))) /
      ((n2 - 6) * sqrt(n1 * (n1 + n2 - 2)))
  }  else {
    stop("Skewness is undefined for F distribution with df2 <= 6.")
  }

})

#' @rdname Fisher
setMethod("kurt",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (n2 > 8) {
    12 * (n1 * (5 * n2 - 22) * (n1 + n2 -2) + (n2 - 4) * (n2 - 2) ^ 2) /
      n1 * (n2 - 6) * (n2 - 8) * (n1 + n2 - 2)
  }  else {
    stop("Kurtosis is undefined for F distribution with df2 <= 8.")
  }

})

#' @rdname Fisher
setMethod("entro",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  lgamma(n1 / 2) + lgamma(n2 / 2) - lgamma((n1 + n2) / 2) +
    (1 - n1 / 2) * digamma(1 + n1 / 2) - (1 + n2 / 2) * digamma(1 + n2 / 2) +
    ((n1 + n2) / 2) * digamma((n1 + n2) / 2) + log(n2) - log(n1)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
#' @export
llf <- function(x, df1, df2) {
  ll(x, prm = c(df1, df2), distr = Fisher())
}

#' @rdname Fisher
setMethod("ll",
          signature  = c(x = "numeric", prm = "numeric", distr = "Fisher"),
          definition = function(x, prm, distr) {

  d1 <- prm[1]
  d2 <- prm[2]
  n <- length(x)
  s <- sum(log(x))
  t <- sum(log(d1 * x + d2))

  (n * d1 * log(d1) + n * d2 * log(d2) + d1 * s - (d1 + d2) * t) / 2 -
    s - n * lbeta(d1 / 2, d2 / 2)

})
