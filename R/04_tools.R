# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools                                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## General                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

loading_bar <- function(total) {
  frm <- "Processing [:bar] :percent | Remaining: :eta | Elapsed: :elapsedfull"
  progress::progress_bar$new(format = frm, total = total, clear = FALSE)
}

seqcol <- function(x) {
  seq_along(x[1, ])
}

seqrow <- function(x) {
  seq_along(x[ , 1])
}

is_pos <- function(x) {
  all(is.finite(x)) && all(x > 0)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Structures             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get the available moment methods for a distribution
get_moment_methods <- function(x) {

  # All available moments
  mom <- c("mean", "median", "mode", "var", "sd", "skew", "kurt",
           "entro", "finf")

  # Get class methods
  df_meth <- attr(methods(class = class(x)), "info")
  meth <- df_meth[df_meth$from == "estim", ]$generic

  mom[mom %in% meth]

}

# Turn distribution name from character to an S4 class
get_distr_class <- function(distr) {

  distr <- paste(toupper(substr(distr, 1, 1)),
                 substr(tolower(distr), 2, nchar(distr)), sep = "")

  if (distr == "Gamma") {
    distr <- "Gam"
  } else if (distr == "Mgamma") {
    distr <- "MGamma"
  }

  new(Class = distr)

}

# Turn an S4 object to a list
s4_to_list <- function(object) {

  # Get the slot names
  names <- methods::slotNames(class(object))

  # Initialize an empty list to store slot values
  y <- list()

  # Loop through the slot names and extract slot values
  for (name in names) {
    y[[name]] <- methods::slot(object, name)
  }

  y

}

# Get the parameters of a distribution
get_params <- function(D, list = TRUE) {
  params <- s4_to_list(D)
  params["name"] <- NULL
  params["ncp"] <- NULL

  if (!list) {
    params <- unlist(params)
  }

  params
}

# Get the unknown parameters of a distribution
get_unknown_params <- function(D, list = TRUE) {

  params <- get_params(D)

  if (is(D, "Binom") || is(D, "Nbinom") || is(D, "Multinom")) {
    params$size <- NULL
  }

  if (!list) {
    params <- unlist(params)
  }

  params
}

# Update the distribution parameters
update_params <- function(D, prm, i) {

  # Position of parameter (e.g. the third element of the shape parameter vector)
  if (is.null(prm$pos)) {
    prm$pos <- 1
  }

  params <- s4_to_list(D)
  params[[prm$name]][prm$pos] <- prm$val[i]

  do.call("new", c(params, Class = class(D)))

}

array_to_df <- function(x) {

  dn <- dimnames(x)
  names(dn) <- names(dimnames(x))

  df <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE)
  df$Value <- as.vector(x)

  df

}

set1of1 <- function(x, i) {
  x
}

set1of2 <- function(x, i) {
  x[i, ]
}

set1of3 <- function(x, i) {
  x[i, , ]
}

set2of3 <- function(x, i) {
  x[, i, ]
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Statistics             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Biased variance
bvar <- function(x) {
  ((length(x) - 1) / length(x)) * var(x)
}

# Biased standard deviation
bsd <- function(x) {
  sqrt(bvar(x))
}

rowVar <- function(x) {
  rowMeans(x ^ 2) - rowMeans(x) ^ 2
}

colVar <- function(x) {
  colMeans(x ^ 2) - colMeans(x) ^ 2
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Gamma Function         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Polygamma Functions
#'
#' @description
#' This set of functions revolve around the polygamma functions, i.e. the
#' derivatives of the gamma function.
#'
#' @param x,y numeric. The points to evaluate the function.
#' @param p integer. The p-variate Gamma function.
#' @param log logical. Should the logarithm of the result be returned?
#' @param ... extra arguments passed to `optim()`.
#'
#' @describeIn idigamma inverse digamma function.
#'
#' @return numeric. The evaluated function.
#'
#' @export
#'
#' @details
#' These functions are needed for the beta and gamma distribution families (and
#' their multivariate analogs, e.g. the Dirichlet). They appear in the
#' estimation and the asymptotic variance-covariance matrix of the MLE and the
#' SAME.
#'
#' The `idigamma()` function implements the inverse of the digamma function
#' \eqn{\psi}. It is a numerical approximation based on the L-BFGS-U
#' quasi-Newton algorithm. Specifically, `idigamma()` makes a call to `optim()`
#' in order to to solve the equation \eqn{\psi(x) = y}; more accurately, to find
#' the minimum of \eqn{f(x) = \log\Gamma(x) - xy}, whose derivative is
#' \eqn{f'(x) = \psi(x) - y}. The optimization is restricted within the tight
#' bounds derived by Batir (2017). The function is vectorized.
#'
#' @references
#' Necdet Batir (2017), INEQUALITIES FOR THE INVERSES OF THE POLYGAMMA FUNCTIONS
#' https://arxiv.org/pdf/1705.06547
#'
#' Oikonomidis, I. & Trevezas, S. (2023), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @examples
#' idigamma(2)
#' Ddigamma(2, 3)
#' Dtrigamma(2, 3)
#' gammap(1:3, 3)
idigamma <- function(x, ...) {

  unlist(lapply(x, FUN = function(x) {

    # Tight bounds derived in https://arxiv.org/pdf/1705.06547
    l <- 1 / log(1 + exp(- x))
    u <- exp(x) + 0.5

    # Quasi-Newton optimization
    optim(par = (l + u) / 2,
          fn = function(y) {lgamma(y) - x * y},
          gr = function(y) {digamma(y) - x},
          method = "Brent",
          lower = l,
          upper = u, ...)$par
  }))

}

#' @describeIn idigamma digamma difference function.
#' @export
Ddigamma <- function(x, y) {
  digamma(x) - digamma(y)
}

#' @describeIn idigamma trigamma difference function.
#' @export
Dtrigamma <- function(x, y) {
  trigamma(x) - trigamma(y)
}

#' @describeIn idigamma p-variate gamma function
#' @export
gammap <- function(x, p, log = FALSE) {

  g <- x

  for (i in seq_along(x)) {
    g[i] <- (p * (p - 1) / 4) * log(pi)  + sum(lgamma(x[i] + (1 - 1:p) / 2))
  }

  if (!log) { g <- exp(g) }
  g

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Matrix Algebra         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# mat_to_vec <- function(Sigma) {
#
#   x <- Matrix::Cholesky(Sigma, perm = FALSE)
#   D <- Matrix::expand1(x, "D")
#   L <- Matrix::expand1(x, "L1")
#   c(diag(D), L[lower.tri(L)])
#
# }

is_symmetric <- function(x) {
  sum(x == t(x)) == (nrow(x) ^ 2)
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

is_integer <- function(x) {
  identical(x, round(x))
}

is_natural <- function(x) {
  is_integer(x) && (x > 0)
}

inv2x2 <- function(x) {
  det <- x[1, 1] * x[2, 2] - x[1, 2] * x[2, 1]

  if (det == 0) {
    return("The matrix is singular, its inverse does not exist.")
  } else {
    inv <- matrix(c(x[2, 2], -x[2, 1], -x[1, 2], x[1, 1]), nrow = 2)
    inv <- inv / det
    dimnames(inv) <- dimnames(x)
    return(inv)
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multivariate Gamma     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

tdiff <- function(x) {
  t(diff(t(x)))
}

fd <- function(x) {
  if (is.matrix(x)) {
    return(cbind(x[, 1], tdiff(x)))
  } else if (is.vector(x)) {
    return(c(x[1], diff(x)))
  } else {
    stop("x must be an atomic vector or a matrix.")
  }
}

gendir <- function(x) {
  z <- fd(x)
  if (is.matrix(x)) {
    return(z / x[, ncol(x)])
  } else if (is.vector(x)) {
    return(z / x[length(x)])
  } else {
    stop("x must be an atomic vector or a matrix.")
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generalized Gamma      ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

stacy_to_prentice <- function(a, b, k) {
  c("mu" = log(a) + digamma(k) / b,
    "sigma" = 1 / (b * sqrt(k)),
    "Q" = 1 / sqrt(k))
}

prentice_to_stacy <- function(mu, sigma, q) {
  k <- 1 / q ^ 2
  b <- q / sigma
  a <- exp(mu - digamma(k) / b)
  c("a" = a, "b" = b, "k" = k)
}
