# Tools ------------------------------------------------------------------------
#_______________________________________________________________________________

# Turn distribution name from character to an S4 class
get_distr_class <- function(distr) {

  distr <- paste(toupper(substr(distr, 1, 1)),
                 substr(tolower(distr), 2, nchar(distr)), sep = "")

  if (distr == "Gamma") {
    distr <- "Gammad"
  } else if (distr == "Mgamma") {
    distr <- "MGamma"
  }

  new(Class = distr)

}

# Gamma Function ----

idigamma <- function(x) {
  distr::igamma(x)
}

# Digamma Difference
Ddigamma <- function(x, y) {
  digamma(x) - digamma(y)
}

# Trigamma Difference
Dtrigamma <- function(x, y) {
  trigamma(x) - trigamma(y)
}

# p-variate Gamma function
gammap <- function(x, p) {
  g <- x
  for (i in seq_along(x)) {
    g[i] <- pi ^ (p * (p - 1) / 4) * prod(gamma(x[i] + (1 - 1:p) / 2))
  }
  g
}

# logarithm of p-variate Gamma function
lgammap <- function(x, p) {
  g <- x
  for (i in seq_along(x)) {
  g[i] <- (p * (p - 1) / 4) * log(pi)  + sum(lgamma(x[i] + (1 - 1:p) / 2))
  }
  g
}

# Matrix Algebra ----

Matrix <- function(...) {
  Matrix::Matrix(...)
}

nearPD <- function(x) {
  Matrix::nearPD(x)$mat
}

vec_to_mat <- function(prm) {

  p <- 0.5 * (- 1 + sqrt(1 + 8 * length(prm)))
  d <- prm[1:p]
  L <- diag(p)
  L[lower.tri(L)] <- prm[(p + 1):length(prm)]

  L %*% diag(d) %*% t(L)

}

mat_to_vec <- function(Sigma) {

  x <- Matrix::Cholesky(A, perm = FALSE)
  D <- Matrix::expand1(x, "D")
  L <- Matrix::expand1(x, "L1")
  c(diag(D), L[lower.tri(L)])

}

is_symmetric <- function(x) {
  sum(x == t(x)) == (nrow(x)^2)
}

is_pd <- function(x) {

  if (!is_symmetric(x)) {
    stop("x is not a symmetric matrix.")
  }

  eigs <- eigen(x, symmetric = TRUE)$values
  if (any(is.complex(eigs))) {
    return(FALSE)
  }

  if (all(eigs > 0)) {
    pd <- TRUE
  } else {
    pd <- FALSE
  }

  pd

}

# Multivariate Gamma ----

fd <- function(x) {
  if (is.matrix(x)) {
    return(rbind(x[1,], diff(x)))
  } else if (is.vector(x)) {
    return(c(x[1], diff(x)))
  } else {
    stop("x must be an atomic vector or a matrix.")
  }
}

gendir <- function(x) {
  z <- fd(x)
  if (is.matrix(x)) {
    return(t(t(z) / x[nrow(x), ]))
  } else if (is.vector(x)) {
    return(z / x[length(x)])
  } else {
    stop("x must be an atomic vector or a matrix.")
  }
}
