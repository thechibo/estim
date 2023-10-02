#-------------------------------------------------------------------------------
# Tools
#-------------------------------------------------------------------------------

get_distr_class <- function(distr) {

  distr <- stringr::str_to_title(distr)

  if (distr == "Gamma") {
    distr <- "Gammad"
  } else if (distr == "Mgamma") {
    distr <- "MGamma"
  }

  new(Class = distr)

}

trigammaInverse <- function(x) {
  if (!is.numeric(x))
    stop("Non-numeric argument to mathematical function")
  if (length(x) == 0)
    return(numeric(0))
  omit <- is.na(x)
  if (any(omit)) {
    y <- x
    if (any(!omit))
      y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x < 0)
  if (any(omit)) {
    y <- x
    y[omit] <- NaN
    warning("NaNs produced")
    if (any(!omit))
      y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x > 1e+07)
  if (any(omit)) {
    y <- x
    y[omit] <- 1/sqrt(x[omit])
    if (any(!omit))
      y[!omit] <- Recall(x[!omit])
    return(y)
  }
  omit <- (x < 1e-06)
  if (any(omit)) {
    y <- x
    y[omit] <- 1/x[omit]
    if (any(!omit))
      y[!omit] <- Recall(x[!omit])
    return(y)
  }
  y <- 0.5 + 1/x
  iter <- 0
  repeat {
    iter <- iter + 1
    tri <- trigamma(y)
    dif <- tri * (1 - tri/x)/psigamma(y, deriv = 2)
    y <- y + dif
    if (max(-dif/y) < 1e-08)
      break
    if (iter > 50) {
      warning("Iteration limit exceeded")
      break
    }
  }
  y
}

gammap <- function(x, p) {
  g <- x
  for (i in seq_along(x)) {
    g[i] <- pi ^ (p * (p - 1) / 4) * prod(gamma(x[i] + (1 - 1:p) / 2))
  }
  g
}

lgammap <- function(x, p) {
  g <- x
  for (i in seq_along(x)) {
  g[i] <- (p * (p - 1) / 4) * log(pi)  + sum(lgamma(x[i] + (1 - 1:p) / 2))
  }
  g
}

vec_to_mat <- function(prm) {

  p <- 0.5 * (- 1 + sqrt(1 + 8 * length(prm)))
  d <- prm[1:p]
  L <- diag(p)
  L[lower.tri(L)] <- prm[(p+1):length(prm)]

  L %*% diag(d) %*% t(L)
}

mat_to_vec <- function(Sigma) {

  x <- fastmatrix::ldl(Sigma)
  d <- x$d
  L <- x$lower
  c(d, L[lower.tri(L)])
}

is_pd <- function(x) {
  LaplacesDemon::is.positive.definite(x)
}
