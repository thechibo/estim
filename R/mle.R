# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Maximum Likelihood Estimator                                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#' @return matrix. The asymptotic covariance matrix of the estimator.
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Beta"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  x  <- as.matrix(x)
  mle(x, distr)[ , 1]

})

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Beta"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  dn <- list(prm = c("shape1", "shape2"), sam = seqcol(x))
  y  <- matrix(0, nrow = 2, ncol = ncol(x), dimnames = dn)

  for (j in seqcol(x)) {
    y[, j] <- optim(par = do.call(par0, list(x = x[ , j], distr = distr)),
                    fn = ll,
                    x = x[ , j],
                    distr = distr,
                    method = method,
                    lower = lower,
                    upper = upper,
                    control = list(fnscale = -1))$par
  }

  y

})

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "Beta"),
          definition = function(distr) {

  a <- distr::shape1(distr)
  b <- distr::shape2(distr)

  p1a  <- trigamma(a)
  p1b  <- trigamma(b)
  p1   <- trigamma(a + b)
  D    <- 1 / (p1a * p1b - (p1a + p1b) * p1)

  matrix(D * c(p1b - p1, p1, p1, p1a - p1), nrow = 2, ncol = 2)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gamma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "numeric", distr = "Gammad"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  x  <- as.matrix(x)
  mle(x, distr)[ , 1]

})

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Gammad"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = c(1e-5, 1e-5),
                                upper = c(Inf, Inf)) {

  dn <- list(prm = c("shape", "scale"), sam = seqcol(x))
  y  <- matrix(0, nrow = 2, ncol = ncol(x), dimnames = dn)

  for (j in seqcol(x)) {
    y[, j] <- optim(par = do.call(par0, list(x = x[ , j], distr = distr)),
                    fn = ll,
                    x = x[ , j],
                    distr = distr,
                    method = method,
                    lower = lower,
                    upper = upper,
                    control = list(fnscale = -1))$par
  }

  y

})

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "Gammad"),
          definition = function(distr) {

  # Ye and Chen (2017)

  a <- distr::shape(distr)
  b <- distr::param(distr)@scale

  psi1a <- trigamma(a)
  d <- 1 / (a * psi1a - 1)

  v11 <- a
  v21 <- - b
  v22 <- b ^ 2 * psi1a

  matrix(d * c(v11, v21, v21, v22), nrow = 2, ncol = 2)

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dirichlet
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Dirichlet"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  tx  <- rowMeans(log(x))

  a0 <- optim(par = sum(do.call(par0, list(x = x, distr = distr))),
              fn = lloptim,
              gr = dlloptim,
              tx = tx,
              distr = distr,
              method = method,
              lower = lower,
              upper = upper,
              control = list(fnscale = -1))$par

  shape <- idigamma(digamma(a0) + tx)

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
  as.matrix(nearPD(D))

})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate Gamma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "MGamma"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf) {

  k <- nrow(x)
  logz <- rowMeans(log(fd(x)))
  xk <- mean(x[k, ])

  par <- optim(par = sum(do.call(par0, list(x = x, distr = distr))[1:k]),
               fn = lloptim,
               gr = dlloptim,
               tx = c(logz, xk),
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

  as.matrix(nearPD(D))

})
