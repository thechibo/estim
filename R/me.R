#-------------------------------------------------------------------------------
# Moment Estimator
#-------------------------------------------------------------------------------

#' Moment Estimator
#'
#' @inherit mle params return examples
#' @export
setGeneric("me", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("me") })

#' @inherit acov_mle title params return examples
#' @export
setGeneric("acov_me", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_me") })

#' @rdname me
setMethod("me",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

  me(x, get_distr_class(distr), ...)

})

# Beta         ----

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

# Gamma        ----

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

# Dirichlet    ----

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
          definition = function(distr) {


})

# Matrix Gamma ----

#' @rdname me
setMethod("me",
          signature  = c(x = "array", distr = "MGamma"),
          definition = function(x, distr) {

  x2 <- x
  for (k in 1:dim(x)[3]) {
    x2[, , k] <- x2[, , k] %*% x2[, , k]
  }

  m  <- apply(x, FUN = mean, MAR = 1:2)
  m2  <- apply(x2, FUN = mean, MAR = 1:2)
  m2minv <- m2 %*% solve(m)
  m2minv[lower.tri(m2minv)] <- m2minv[upper.tri(m2minv)]

  Sigma <- 2 * m2minv - 2 * m
  diag(Sigma) <- diag(Sigma) - sum(diag(Sigma)) / (nrow(Sigma) + 1)
  shape <- mean(diag(m %*% solve(Sigma)))

  list(shape = shape, Sigma = Sigma)

})

#' @rdname acov_me
setMethod("acov_me",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {


})
