# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dirichlet Distribution                                                    ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title The Dirichlet Distribution
#'
#' @description
#' Density function and random generation for the Dirichlet distribution with
#' parameter vector `shape`.
#'
#' @param x numeric. The quantile vector.
#' @param shape numeric. The parameter vector.
#' @param log logical. If `TRUE`, probabilities `p` are given as `log(p)`.
#' @param n numeric. The number of observations.
#'
#' @return `dDirichlet` returns a numeric vector (the evaluated density
#' function). `rDirichlet` returns a matrix with `length(shape)` rows and `n`
#' columns.
#'
#' @export
#'
#' @examples \dontrun{
#' # Classic R Stats Format
#' dDirichlet(c(0.3, 0.7), shape = c(2, 3))
#' set.seed(1)
#' rDirichlet(10, shape = c(2, 3))
#'
#' # S4 Distribution Class
#' D <- Dirichlet(shape = c(2, 3))
#' d(D)(c(0.3, 0.7))
#' set.seed(1)
#' r(D)(10)
#' }
dDirichlet <- function(x, shape, log = FALSE) {

  if (length(x) != length(shape)) {
    stop("The lengths of x (", length(x), ") and shape (",
         length(shape), ") must be equal.")
  }

  ld <- lgamma(sum(shape)) - sum(lgamma(shape)) + sum((shape - 1) * log(x))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname dDirichlet
rDirichlet <- function(n, shape) {

  k <- length(shape)
  x <- matrix(nrow = n, ncol = k)
  for (j in 1:k) {
    x[, j] <- stats::rgamma(n, shape[j], 1)
  }

  t(x / rowSums(x))

}

setClass("DirichletParameter",
         representation = representation(shape = "numeric"),
         prototype = prototype(shape = c(1, 1),
                               name = gettext("Parameter of a Dirichlet distribution")),
         contains = "Parameter"
)

#' @title Dirichlet Distribution S4 Class
#'
#' @slot shape numeric. The parameter vector.
#'
#' @return An object of class `Dirichlet`.
#' @export
#'
#' @inherit dDirichlet examples
Dirichlet <- setClass("Dirichlet",
                      slots = list(shape = "numeric"),
                      prototype = prototype(
                        r = function(n) {
                          rDirichlet(n, shape = c(1, 1))
                        },
                        d = function(x, log = FALSE) {
                          dDirichlet(x, shape = c(1, 1), log = log)
                        },
                        param = new("DirichletParameter"),
                        .logExact = TRUE,
                        .lowerExact = TRUE
                      ),
                      contains = "AbscontDistribution"
)

# Access methods
setMethod("shape", "DirichletParameter", function(object) object@shape)

# Replace methods
setReplaceMethod("shape", "DirichletParameter",
                 function(object, value){ object@shape <- value; object})

setValidity("DirichletParameter", function(object){
  if (any(shape(object) <= 0)) {
    stop("shape has to be positive")
  } else {
    return(TRUE)
  }
})

# wrapped access methods
setMethod("shape", "Dirichlet",
          function(object) { distr::shape(distr::param(object)) })

# wrapped replace methods
setMethod("shape<-", "Dirichlet",
          function(object, value) { new("Dirichlet", shape = value(object)) })

setMethod("initialize", "Dirichlet",
          function(.Object, shape = c(1, 1)) {
            .Object@img <- new("Reals")
            .Object@param <- new("DirichletParameter", shape = shape)
            .Object@r <- function(n) {}
            .Object@d <- function(x, log = FALSE) {}
            body(.Object@r) <- substitute(
              { rDirichlet(n, shape = shapeSub) },
              list(shapeSub = shape)
            )
            body(.Object@d) <- substitute(
              { dDirichlet(x, shape = shapeSub, log = log) },
              list(shapeSub = shape)
            )
            .Object@.withSim   <- FALSE
            .Object@.withArith <- FALSE
            .Object
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "matrix", distr = "Dirichlet"),
          definition = function(prm, x, distr) {

  sum(apply(x, MARGIN = 2, FUN = dDirichlet, shape = prm, log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dirichlet"),
          definition = function(par, tx, distr) {

  a <- idigamma(digamma(par) + tx)
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * tx)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dirichlet"),
          definition = function(par, tx, distr) {

  # Shape parameters (a_i) as a function of a0
  a <- idigamma(digamma(par) + tx)

  # a_i derivative wrt a0
  da <- trigamma(par) / trigamma(a)

  # lloptim derivative wrt a0 (par)
  digamma(sum(a)) * sum(da) - sum(digamma(a) * da) + sum(tx * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MLE                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Dirichlet"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- rowMeans(log(x))

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

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "Dirichlet"),
          definition = function(distr) {

  a <- distr::shape(distr)
  k <- length(a)

  D <- solve(diag(trigamma(a)) - matrix(trigamma(sum(a)), k, k))
  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("shape", seq_along(a))
  colnames(D) <- paste0("shape", seq_along(a))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "Dirichlet"),
          definition = function(x, distr) {

  m  <- rowMeans(x)
  m2 <- rowMeans(x ^ 2)
  shape  <- m * (m - m2) / (m2 - m ^ 2)

  names(shape) <- paste0("shape", seq_along(shape))
  shape

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "Dirichlet"),
          definition = function(distr, comp = FALSE) {

  a <- distr::shape(distr)
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)

  A1 <- diag(a0 * (2 * a0 + 1) * (a + 1) / b)
  A2 <- diag(- a0 * (a0 + 1) ^ 2 / b)
  A <- Matrix(rbind(A1, A2))

  B11 <- - Matrix(a, k, 1) %*% Matrix(a, 1, k)
  diag(B11) <- a * b
  B11 <- B11 / (a0 ^ 2 * (a0 + 1))
  B11 <- nearPD(B11)

  B12 <- - Matrix(a, k, 1) %*% Matrix(a * (a + 1), 1, k)
  diag(B12) <- a * (a + 1) * b
  B12 <- B12 * 2 / (a0 ^ 2 * (a0 + 1) * (a0 + 2))

  c22 <- - 2 * (2 * a0 + 3) / (a0 ^ 2 * (a0 + 1) ^ 2 * (a0 + 2) * (a0 + 3))
  B22 <- c22 * Matrix(a * (a + 1), k, 1) %*% Matrix(a * (a + 1), 1, k)
  diag(B22) <- (a * (a + 1) * (a + 2) * (a + 3)) /
    (a0 * (a0 + 1) * (a0 + 2) * (a0 + 3)) -
    (a * (a + 1) / (a0 * (a0 + 1))) ^ 2
  B22 <- nearPD(B22)

  B <- rbind(cbind(B11, B12),
             cbind(Matrix::t(B12), B22))
  B <- nearPD(B)

  if (!comp) {
    D <- nearPD(Matrix::t(A) %*% B %*% A)
    D <- as.matrix(D)
    rownames(D) <- paste0("shape", seq_along(a))
    colnames(D) <- paste0("shape", seq_along(a))
    return(D)
  } else {
    return(list(A = as.matrix(A), B = as.matrix(B)))
  }

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "Dirichlet"),
          definition = function(x, distr) {

  m  <- rowMeans(x)
  logm  <- rowMeans(log(x))
  mlogm <- rowMeans(x * log(x))

  shape  <- (length(m) - 1) * m / sum(mlogm - m * logm)

  names(shape) <- paste0("shape", seq_along(shape))
  shape

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "Dirichlet"),
          definition = function(distr, comp = FALSE) {

  # Required variables
  a <- distr::shape(distr)
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)

  # Matrix A

  A1 <- Matrix(Ddigamma(a, a0), k, 1) %*% Matrix(a, 1, k)
  A2 <- Matrix(a, k, 1) %*% Matrix(a, 1, k)
  A3 <- - Matrix(1, k, 1) %*% Matrix(a, 1, k)
  Ik <- diag(k)
  Ok <- 0 * Ik

  A <- a0 * (rbind(A1, A2 / a0, A3) / (k - 1) + rbind(Ik, Ok, Ok))

  # Matrix B

  B11 <- - A2
  diag(B11) <- a * b
  B11 <- B11 / (a0 ^ 2 * (a0 + 1))
  B11 <- nearPD(B11)

  B22 <- trigamma(a) * Ik - Matrix(trigamma(a0), k, k)
  B22 <- nearPD(B22)

  c331 <- Ddigamma(a + 1, a0 + 2)
  B331 <- Matrix(c331, k, 1) %*% Matrix(c331, 1, k)
  c332 <- Ddigamma(a + 1, a0 + 1)
  B332 <- Matrix(c332, k, 1) %*% Matrix(c332, 1, k)

  B33 <- A2 * (B331 - trigamma(a0 + 2)) / (a0 * (a0 + 1)) - A2 * B332 / (a0 ^ 2)
  diag(B33) <- (Ddigamma(a + 2, a0 + 2) ^ 2 + Dtrigamma(a + 2, a0 + 2))  *
    a * (a + 1) / (a0 * (a0 + 1)) - (Ddigamma(a + 1, a0 + 1) * a / a0) ^ 2
  B33 <- nearPD(B33)

  B12 <- Ik / a0 - Matrix(a, k, 1) %*% Matrix(1, 1, k) / a0 ^ 2

  c13 <- a * (Ddigamma(a + 1, a0 + 2) + 1) / (a0 ^ 2 * (a0 + 1))
  B13 <- - Matrix(a, k, 1) %*% Matrix(c13, 1, k)
  diag(B13) <- - Matrix::diag(B13) * b / a

  c231 <- (Ddigamma(a + 1, a0 + 1) + a * trigamma(a + 1)) / a0
  c232 <- (Ddigamma(a + 1, a0 + 1) / a0 + trigamma(a0 + 1)) * a / a0
  B23 <- c231 * Ik - Matrix(1, k, 1) %*% Matrix(c232, 1, k)

  B <- rbind(cbind(B11, B12, B13),
             cbind(Matrix::t(B12), B22, B23),
             cbind(Matrix::t(B13), Matrix::t(B23), B33))
  B <- nearPD(B)

  if (!comp) {
    D <- nearPD(Matrix::t(A) %*% B %*% A)
    D <- as.matrix(D)
    rownames(D) <- paste0("shape", seq_along(a))
    colnames(D) <- paste0("shape", seq_along(a))
    return(D)
  } else {
    return(list(A = as.matrix(A), B = as.matrix(B)))
  }

})
