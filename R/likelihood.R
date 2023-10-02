#-------------------------------------------------------------------------------
# Log-Likelihood
#-------------------------------------------------------------------------------

#' Log-Likelihood
#'
#' @param prm numeric. A vector of the distribution parameters.
#' @param x numeric. A sample under estimation.
#' @param distr A subclass of `Distribution`. The distribution of interest.
#' @param ... extra arguments.
#'
#' @return Numeric. The value of the log-likelihood function.
#' @export
setGeneric("ll", signature = c("prm", "x", "distr"),
           function(prm, x, distr, ...) { standardGeneric("ll") })

# Beta ---

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "numeric", distr = "Beta"),
          definition = function(prm, x, distr) {

  sum(dbeta(x = x, shape1 = prm[1], shape2 = prm[2], log = TRUE))

})

# Bias Corrected log-likelihood
# (Firth, 1993, Cribari-Neto and Vasconcellos, 2010)
#ll = function(prm, x) {
#  p1a = trigamma(prm[1])
#  p1b = trigamma(prm[2])
#  p1  = trigamma(sum(prm))
#  d   = p1a * p1b - p1 * (p1a + p1b)
#  ld = log((length(x) ^ 2) * d) / 2
#  sum(do.call(dbeta, c(list(x = x, log = TRUE), prm))) + ld
#}

# Gamma ---

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "numeric", distr = "Gammad"),
          definition = function(prm, x, distr) {

  sum(dgamma(x = x, shape = prm[1], scale = prm[2], log = TRUE))

})

# Dirichlet ---

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "matrix", distr = "Dirichlet"),
          definition = function(prm, x, distr) {

  sum(apply(x, MAR = 2, FUN = dDirichlet, prm = prm, log = TRUE))

})

# Matrix Gamma ---

#' @rdname ll
setMethod("ll",
          signature  = c(prm = "numeric", x = "array", distr = "MGamma"),
          definition = function(prm, x, distr) {

  Sigma <- vec_to_mat(prm[2:length(prm)])
  sum(apply(x, MAR = 3, FUN = dMGamma, shape = prm[1], Sigma = Sigma, log = TRUE))

})
