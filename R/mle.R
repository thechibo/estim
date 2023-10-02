#-------------------------------------------------------------------------------
# Maximum Likelihood Estimator
#-------------------------------------------------------------------------------

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
#' @return numeric or list. The estimator produced by the sample. In case the
#' distribution parameters hold specific structures (e.g. the Matrix Gamma
#' parameters are a positive real number and a positive definite matrix), a
#' list is returned instead.
#'
#' @importClassesFrom distr Beta Gammad
#' @export
#'
#' @examples \dontrun{
#' # Distribution
#' d_beta <- Beta(shape1 = 1, shape2 = 1.5)
#'
#' # Simulation
#' set.seed(2)
#' x <- r(d_beta)(50)
#'
#' # Estimation
#' me(x, "Beta")
#' same(x, "Beta")
#' mle(x, "Beta")
#'
#' # Asymptotic Covariance Matrix
#' acov_me(d_beta)
#' acov_same(d_beta)
#' acov_mle(d_beta)
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

# Beta         ----

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

  dn <- list(prm = c("shape1", "shape2"), sam = 1:ncol(x))
  y  <- matrix(0, nrow = 2, ncol = ncol(x), dimnames = dn)

  for (j in 1:ncol(x)) {
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

# Gamma        ----

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

  dn <- list(prm = c("shape", "scale"), sam = 1:ncol(x))
  y  <- matrix(0, nrow = 2, ncol = ncol(x), dimnames = dn)

  for (j in 1:ncol(x)) {
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

# Dirichlet    ----

#' @rdname mle
setMethod("mle",
          signature  = c(x = "matrix", distr = "Dirichlet"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = rep(1e-5, times = nrow(x)),
                                upper = rep(Inf, times = nrow(x))) {

  shape <- optim(par = do.call(par0, list(x = x, distr = distr)),
                fn = ll,
                x = x,
                distr = distr,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  names(shape) <- paste0("shape", seq_along(shape))
  shape

})

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "Dirichlet"),
          definition = function(distr) {


})

# Matrix Gamma ----

#' @rdname mle
setMethod("mle",
          signature  = c(x = "array", distr = "MGamma"),
          definition = function(x, distr,
                                par0 = same,
                                method = "L-BFGS-B",
                                lower = NULL,
                                upper = NULL) {

  p <- dim(x)[1]

  if (is.null(lower)) {
    lower <- c(0.5 * (p - 1), rep(0.01, times = p), rep(-Inf, times = 0.5 * (p ^ 2 - p)))
  }
  if (is.null(upper)) {
    upper <- rep(Inf, times = 0.5 * (p ^ 2 + p + 4))
  }

  # Initialize
  par0 <- do.call(par0, list(x = x, distr = distr))

  if (!is_pd(par0$Sigma)) {
    par0$Sigma <- matrix(0, nrow = p, ncol = p)
    diag(par0$Sigma) <- rep(1, times = p)
    warning("Setting Sigma to Identity matrix for initialization.")
  }

  par0 <- c(par0$shape, mat_to_vec(par0$Sigma))

  prm <- optim(par = par0,
                fn = ll,
                x = x,
                distr = distr,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  list(shape = prm[1], Sigma = vec_to_mat(prm[2:length(prm)]))

})

#' @rdname acov_mle
setMethod("acov_mle",
          signature  = c(distr = "MGamma"),
          definition = function(distr) {



})
