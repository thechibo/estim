#-------------------------------------------------------------------------------
# Score-Adjusted Moment Estimator
#-------------------------------------------------------------------------------

#' Score-Adjusted Moment Estimator
#'
#' @inherit mle params return examples
#' @export
setGeneric("same", signature = c("x", "distr"),
           function(x, distr, ...) { standardGeneric("same") })

#' @inherit acov_mle title params return examples
#' @export
setGeneric("acov_same", signature = c("distr"),
           function(distr, ...) { standardGeneric("acov_same") })

#' @rdname same
setMethod("same",
          signature  = c(x = "ANY", distr = "character"),
          definition = function(x, distr, ...) {

  same(x, get_distr_class(distr), ...)

})

# Beta         ----

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

# Gamma         ----

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

# Dirichlet    ----

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
          definition = function(distr) {


})

# Matrix Gamma ----

#' @rdname same
setMethod("same",
          signature  = c(x = "array", distr = "MGamma"),
          definition = function(x, distr) {

  ldx  <- apply(x, FUN = function(x) {log(det(x))}, MAR = 3)
  xldx <- x
  for (k in 1:dim(x)[3]) {
    xldx[, , k] <- x[, , k] * ldx[k]
  }

  m  <- apply(x, FUN = mean, MAR = 1:2)
  mldx <- mean(ldx)
  mxldx  <- apply(xldx, FUN = mean, MAR = 1:2)

  Sigma <- mxldx - m * mldx
  shape <- mean(diag(m %*% solve(Sigma)))
  list(shape = shape, Sigma = Sigma)

})

#' @rdname acov_same
setMethod("acov_same",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {


})

