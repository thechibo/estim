# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MGamma Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title The Multivariate Gamma Distribution
#'
#' @description
#' Density function and random generation for the Multivariate Gamma
#' distribution with shape parameter vector `shape` and scale parameter `scale`.
#'
#' @param x numeric. The quantile vector.
#' @param shape numeric. The shape parameter vector.
#' @param scale numeric. The scale parameter vector.
#' @param log logical. If `TRUE`, probabilities `p` are given as `log(p)`.
#' @param n numeric. The number of observations.
#'
#' @return `dMGamma` returns a numeric vector (the evaluated density function).
#' `rMGamma` returns a matrix with `length(shape)` rows and `n` columns.
#'
#' @export
#'
#' @examples \dontrun{
#' # Classic R Stats Format
#' dMGamma(c(4, 6), shape = c(2, 3), scale = 2)
#' set.seed(1)
#' rMGamma(10, shape = c(2, 3), scale = 2)
#'
#' # S4 Distribution Class
#' D <- MGamma(shape = c(2, 3), scale = 2)
#' d(D)(c(4, 6))
#' set.seed(1)
#' r(D)(10)
#' }
dMGamma <- function(x, shape, scale, log = FALSE) {

  if (length(x) != length(shape)) {
    stop("The lengths of x (", length(x), ") and shape (",
         length(shape), ") must be equal.")
  }

  z <- fd(x)
  xk <- x[length(x)]
  a0 <- sum(shape)

  ld <- sum(shape * log(z)) - a0 * log(scale) - sum(lgamma(shape)) - xk / scale

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname dMGamma
rMGamma <- function(n, shape, scale) {

  k <- length(shape)
  x <- matrix(nrow = k, ncol = n)
  for (i in 1:k) {
    x[i, ] <- stats::rgamma(n, shape[i], scale = scale)
  }

  apply(x, 2, cumsum)

}

setClass("MGammaParameter",
         representation = representation(shape = "numeric", scale = "numeric"),
         prototype = prototype(shape = c(1, 1), scale = 1,
                               name = gettext("Parameter of a MGamma distribution")),
         contains = "Parameter"
)

#' @title Multivariate Gamma Distribution S4 Class
#'
#' @slot shape numeric. The shape parameter vector.
#' @slot scale numeric. The scale parameter vector.
#'
#' @return An object of class `MGamma`.
#' @export
#'
#' @inherit dMGamma examples
MGamma <- setClass("MGamma",
                   slots = list(shape = "numeric", scale = "numeric"),
                   prototype = prototype(
                     r = function(n) {
                       rMGamma(n, shape = c(1, 1), scale = 1)
                     },
                     d = function(x, log = FALSE) {
                       dMGamma(x, shape = c(1, 1), scale = 1, log = log)
                     },
                     param = new("MGammaParameter"),
                     .logExact = TRUE,
                     .lowerExact = TRUE
                   ),
                   contains = "AbscontDistribution"
)

# Access methods
setMethod("shape", "MGammaParameter",
          function(object) object@shape)

setMethod("scale", signature = c(x = "MGammaParameter"),
          function(x, center = TRUE, scale = TRUE) x@scale)

# Replace Methods
setReplaceMethod("shape", "MGammaParameter",
                 function(object, value){ object@shape <- value; object})

setReplaceMethod("scale", signature = c(object = "MGammaParameter"),
                 function(object, value){ object@scale <- value; object})

setValidity("MGammaParameter", function(object){
  if (any(shape(object) <= 0)) {
    stop("shape has to be positive")
  } else if (object@scale <= 0) {
    stop("scale has to be positive")
  } else {
    return(TRUE)
  }
})

# wrapped access methods
setMethod("shape", "MGamma",
          function(object) { distr::shape(distr::param(object)) })

setMethod("scale", signature = c(x = "MGamma"),
          function(x, center = TRUE, scale = TRUE) {
            distr::param(x)@scale
          })

# wrapped replace methods
setMethod("shape<-", "MGamma",
          function(object, value) { new("MGamma", shape = value(object)) })

setMethod("scale<-", "MGamma",
          function(object, value) {
            new("MGamma", shape = distr::shape(object), scale = value)
          })

setMethod("initialize", "MGamma",
          function(.Object, shape = c(1, 1), scale = 1) {
            .Object@img <- new("Reals")
            .Object@param <- new("MGammaParameter", shape = shape, scale = scale)
            .Object@r <- function(n) {}
            .Object@d <- function(x, log = FALSE) {}
            body(.Object@r) <- substitute(
              { rMGamma(n, shape = shapeSub, scale = scaleSub) },
              list(shapeSub = shape, scaleSub = scale)
            )
            body(.Object@d) <- substitute(
              { dMGamma(x, shape = shapeSub, scale = scaleSub, log = log) },
              list(shapeSub = shape, scaleSub = scale)
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
          signature  = c(prm = "numeric", x = "matrix", distr = "MGamma"),
          definition = function(prm, x, distr) {

  k <- length(prm)
  sum(apply(x,
            MARGIN = 2,
            FUN = dMGamma,
            shape = prm[1:(k - 1)], scale = prm[k], log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "MGamma"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz - log(b))

  - sum(a) * log(b) - sum(lgamma(a)) - xk / b + sum((a - 1) * logz)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "MGamma"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz -log(b))

  db <- - xk / par ^ 2
  da <- 1 / (par * trigamma(a))

  - sum(da) * log(b) - sum(a) * db / b - sum(digamma(a) * da) +
    xk * db / b ^ 2 + sum(logz * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MLE                    ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  k <- nrow(x)
  logz <- rowMeans(log(fd(x)))
  xk <- mean(x[k, ])
  tx <- c(logz, xk)

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))[1:k]),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  b <- xk / par
  a <- idigamma(logz - log(b))

  par <- c(a, b)

  names(par) <- c(paste0("shape", 1:k), "scale")
  par

})

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)
  k <- length(a)
  a0 <- sum(a)

  D <- solve(rbind(cbind(diag(trigamma(a)), matrix(1 / b, nrow = k, ncol = 1)),
                   c(rep(1 / b, k), a0 / b ^ 2)))

  D <- as.matrix(nearPD(D))
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ME                     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  z <- fd(x)
  mz <- rowMeans(z)

  scale <- mean(rowVar(z) / mz)
  shape <- mz / scale

  c(shape = shape, scale = scale)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  # Preliminaries
  a <- distr::shape(distr)
  b <- distr::scale(distr)
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A11 <- (Matrix(2 + 1 / a, k, 1) %*% Matrix(a, 1, k) /k + diag(1, k, k)) / b
  A21 <- - Matrix(1 / a, k, 1) %*% Matrix(a, 1, k) / (k * b ^ 2)
  A12 <- - (2 + 1 / a) / k
  A22 <- 1 / (a * k * b)
  A <- cbind(rbind(A11, A21), Matrix(c(A12, A22), 2 * k, 1))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- 2 * a * (a + 1) * b ^ 4 * ( 2 * a + 3)
  B12 <- 2 * a * (a + 1) * b ^ 3
  B <- rbind(cbind(diag(B11), diag(B12)),
             cbind(diag(B12), diag(B22)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(Matrix::t(A) %*% B %*% A)
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAME                   ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  z <- t(fd(x))

  scale <- mean(diag(stats::cov(z, log(z))))
  shape <- colMeans(z) / scale

  c(shape = shape, scale = scale)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  # Preliminaries
  a <- distr::shape(distr)
  b <- distr::scale(distr)
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A21 <- Matrix(a, k, 1) %*% Matrix(a, 1, k) / k
  A31 <- - Matrix(1, k, 1) %*% Matrix(a, 1, k) / (k * b)
  A12 <- - (digamma(a) + log(b)) / k
  A22 <- - a * b / k
  A32 <- rep(1 / k, k)
  A11 <- (- Matrix(A12, k, 1) %*% Matrix(a, 1, k) + diag(1, k, k)) / b
  A <- cbind(rbind(A11, A21, A31), Matrix(c(A12, A22, A32), 3 * k, 1))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- trigamma(a)
  B33 <- a * (a + 1) * b ^ 2 * (trigamma(a + 2) + (digamma(a + 2) + log(b)) ^ 2) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b)) ^ 2
  B12 <- rep(b, k)
  B13 <- a * (a + 1) * b ^ 2 * (digamma(a + 2) + log(b)) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b))
  B23 <- a * b * (trigamma(a + 1) + (digamma(a + 1) + log(b)) ^ 2) -
    a * b * (digamma(a) + log(b)) * (digamma(a + 1) + log(b))
  B <- rbind(cbind(diag(B11), diag(B12), diag(B13)),
             cbind(diag(B12), diag(B22), diag(B23)),
             cbind(diag(B13), diag(B23), diag(B33)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(Matrix::t(A) %*% B %*% A)
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dirichlet ME           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname me
setMethod("me2",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  w <- gendir(x)
  shape <- me(w, Dirichlet())
  xk <- mean(x[nrow(x), ])
  scale <- xk / sum(shape)

  c(shape, scale = scale)

})

#' @rdname acov_me
setMethod("acov_me2",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)
  a0 <- sum(a)

  acov <- acov_me(Dirichlet(shape = a), comp = TRUE)
  AD <- acov$A
  BD <- acov$B

  A12 <- - matrix(rowSums(AD) * b / a0)

  A <- rbind(cbind(AD, A12),
             c(rep(0, ncol(AD)), 1 / a0))

  B <- rbind(cbind(BD, c(rep(0, ncol(BD)))),
             c(rep(0, nrow(BD)), a0 * b ^ 2))

  D <- nearPD(Matrix::t(A) %*% B %*% A)

  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dirichlet SAME         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same2",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  w <- gendir(x)
  shape <- same(w, Dirichlet())
  xk <- mean(x[nrow(x), ])
  scale <- xk / sum(shape)

  c(shape, scale = scale)

})

#' @rdname acov_same
setMethod("acov_same2",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)
  a0 <- sum(a)

  acov <- acov_same(Dirichlet(shape = a), comp = TRUE)
  AD <- acov$A
  BD <- acov$B

  A12 <- - matrix(rowSums(AD) * b / a0)

  A <- rbind(cbind(AD, A12),
             c(rep(0, ncol(AD)), 1 / a0))

  B <- rbind(cbind(BD, c(rep(0, ncol(BD)))),
             c(rep(0, nrow(BD)), a0 * b ^ 2))

  D <- nearPD(Matrix::t(A) %*% B %*% A)

  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})
