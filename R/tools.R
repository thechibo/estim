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

# Calculate the Trigamma Inverse function
trigammaInverse <- function(x) {

  if (!is.numeric(x)) {
    stop("Non-numeric argument to mathematical function.")
  }

  if (length(x) == 0) {
    return(numeric(0))
  }

  omit <- is.na(x)
  if (any(omit)) {
    y <- x
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }

  omit <- (x < 0)
  if (any(omit)) {
    y <- x
    y[omit] <- NaN
    warning("NaNs produced")
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }

  omit <- (x > 1e+07)
  if (any(omit)) {
    y <- x
    y[omit] <- 1 / sqrt(x[omit])
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }

  omit <- (x < 1e-06)
  if (any(omit)) {
    y <- x
    y[omit] <- 1 / x[omit]
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }

  y <- 0.5 + 1 / x
  iter <- 0

  while (max(- dif / y) > 1e-08) {
    iter <- iter + 1
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2)
    y <- y + dif
    if (iter > 50) {
      warning("Iteration limit exceeded in Trigamma Inverse calculation.")
      break
    }
  }

  y

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
  L[lower.tri(L)] <- prm[(p+1):length(prm)]

  L %*% diag(d) %*% t(L)

}

mat_to_vec <- function(Sigma) {

  # check Matrix::Cholesky
  x <- fastmatrix::ldl(Sigma)
  d <- x$d
  L <- x$lower
  c(d, L[lower.tri(L)])

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

# Statistics ----

s2 <- function(x) {
  ((length(x) - 1) / length(x)) * var(x)
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
