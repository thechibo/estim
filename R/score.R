# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Score Equations                                                           ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setGeneric("lloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("lloptim") })

setGeneric("dlloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("dlloptim") })

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dirichlet              ----

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate Gamma     ----

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

  - sum(da) * log(b) - sum(a) * db / b - sum(digamma(a) * da) + xk * db / b ^ 2 + sum(logz * da)

})
