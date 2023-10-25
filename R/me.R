# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Moment Estimator                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta                   ----

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  me(x, distr)[ , 1]

})

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  d  <- (m - m2) / (m2 - m ^ 2)
  a  <- d * m
  b  <- d * (1 - m)
  dn <- list(prm = c("shape1", "shape2"), sam = 1:ncol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  prd <- a * b
  th  <- a + b
  th2 <- th ^ 2
  s2  <- prd / (th2 * (th + 1))
  s4  <- s2 ^ 2
  m3  <- 2 * (b - a) * s2 / (th * (th + 2))
  m4  <- 3 * prd * (prd * (th + 2) + 2 * (b - a) ^2) / ((th ^ 4) * (th + 1) * (th + 2) * (th + 3))
  d   <- (th + 1) ^ 2 * (th + 2) ^ 2 * s2
  e   <- (th + 1) ^ 3 * (m4 - s4 - m3 ^ 2 / s2) / s2

  s11 <- (a * (a + 1)) ^ 2 / d + a * e / b
  s22 <- (b * (b + 1)) ^ 2 / d + b * e / a
  s12 <- - a * (a + 1) * b * (b + 1) / d + e
  matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma                  ----

#' @rdname me
setMethod("me",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr) {

  x  <- as.matrix(x)
  me(x, distr)[ , 1]

})

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr) {

  m  <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  s2 <- m2 - m ^ 2
  a  <- m ^ 2 / s2
  b  <- s2 / m

  dn <- list(prm = c("shape", "scale"), sam = 1:ncol(x))
  matrix(c(a, b), nrow = 2, byrow = TRUE, dimnames = dn)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)

  s11 <- 2 * a * (a + 1)
  s22 <- b ^ 2 * (2 * a + 3) / a
  s12 <- - 2 * b * (a + 1)
  matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dirichlet              ----

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
  diag(B22) <- (a * (a + 1) * (a + 2) * (a + 3)) / (a0 * (a0 + 1) * (a0 + 2) * (a0 + 3)) - (a * (a + 1) / (a0 * (a0 + 1))) ^ 2
  B22 <- nearPD(B22)

  B <- rbind(cbind(B11, B12),
             cbind(Matrix::t(B12), B22))
  B <- nearPD(B)

  if (!comp) {
    D <- nearPD(Matrix::t(A) %*% B %*% A)
    return(as.matrix(D))
  } else {
    return(list(A = as.matrix(A), B = as.matrix(B)))
  }

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate Gamma     ----

#' @rdname me
setMethod("me",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr) {

  w <- gendir(x)
  shape <- me(w, Dirichlet())
  xk <- mean(x[nrow(x), ])
  scale <- xk / sum(shape)

  c(shape, scale = scale)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {

  a <- distr::shape(distr)
  b <- distr::scale(distr)
  a0 <- sum(a)

  acov <- acov_me(Dirichlet(shape = a), comp = TRUE)
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
