# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Score-Adjusted Moment Estimator
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  same(x, distr)[ , 1]

})

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr) {

  x <- as.matrix(x)
  mx <- colMeans(x)
  mlx <- colMeans(log(x))
  mxlx <- colMeans(x * log(x))
  my <- 1 - mx
  mly <- colMeans(log(1 - x))
  myly <- colMeans((1 - x) * log(1 - x))

  sx <- mxlx - mx * mlx
  sy <- myly - my * mly

  a <- mx / (sx + sy)
  b <- my / (sx + sy)

  dn <- list(prm = c("shape1", "shape2"), sam = 1:ncol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  p1a <- trigamma(a)
  p1b <- trigamma(b)
  m1  <- matrix(c(a ^ 2, prd, prd, b ^ 2), nrow = 2, ncol = 2)
  m2  <- matrix(c(prd, th2 - prd, th2 - prd, prd), nrow = 2, ncol = 2)

  (s2 * th2 * (p1a + p1b) + 1) * m1 - m2 / (th + 1)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  same(x, distr)[ , 1]

})

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr) {

  mx  <- colMeans(x)
  mlx <- colMeans(log(x))
  mxlx <- colMeans(x * log(x))
  b <- mxlx - mx * mlx
  a  <- mx / b

  dn <- list(prm = c("shape", "scale"), sam = 1:ncol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)

  c1 <- 1 + a * trigamma(a + 1)
  c2 <- 1 + a * trigamma(a)

  v11 <- a ^ 2 * c1
  v21 <- - a * b * c1
  v22 <- b ^ 2 * c2

  matrix(c(v11, v21, v21, v22), 2, 2)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dirichlet
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  diag(B33) <- (Ddigamma(a + 2, a0 + 2) ^ 2 + Dtrigamma(a + 2, a0 + 2))  * a * (a + 1) / (a0 * (a0 + 1)) - (Ddigamma(a + 1, a0 + 1) * a / a0) ^ 2
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
    return(as.matrix(D))
  } else {
    return(list(A = as.matrix(A), B = as.matrix(B)))
  }

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate Gamma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname same
setMethod("same",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  w <- gendir(x)
  shape <- same(w, Dirichlet())
  xk <- mean(x[nrow(x), ])
  scale <- xk / sum(shape)

  c(shape, scale = scale)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)
  a0 <- sum(a)

  acov <- acov_same(Dirichlet(shape = a), comp = TRUE)
  AD <- acov$A
  BD <- acov$B

  A12 <- matrix(rowSums(AD) * b / a0)

  A <- rbind(cbind(AD, A12),
             c(rep(0, ncol(AD)), 1 / a0))

  B <- rbind(cbind(BD, c(rep(0, ncol(BD)))),
             c(rep(0, nrow(BD)), a0 * b ^ 2))

  D <- nearPD(Matrix::t(A) %*% B %*% A)

  as.matrix(D)

})
