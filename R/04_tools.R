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
get_params <- function(D) {
  params <- s4_to_list(D)
  params["name"] <- NULL
  params["ncp"] <- NULL
  unlist(params)
}

# Get the parameters of a distribution
get_params_list <- function(D) {
  params <- s4_to_list(D)
  params["name"] <- NULL
  params["ncp"] <- NULL
  params
}

get_unknown_params <- function(D) {
  prm <- get_params(D)
  if (is(D, "Binom")) {
    return(prm["prob"])
  } else if (is(D, "Nbinom")) {
    return(prm["prob"])
  } else if (is(D, "Multinom")) {
    return(prm[-1])
  } else {
    return(prm)
  }
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

# .x <- seq(1e-70, exp(25), length = 1e7)
# .y <- digamma(.x)
# idigamma_approx <- approxfun(.y, .x)
# rm(.x, .y)

## an extensive grid of x-values
.xg <- sort(c(10^(-70:-1),qexp(unique(pmin(seq(0,1,length=5e3)+1e-10,1-1e-10))),qcauchy(seq(0.999,1-1e-10,length=5e3))))
.dxg <- digamma(.xg)
igamma <- approxfun(.dxg,.xg)
rm(.xg,.dxg)

#' @title Polygamma Functions
#'
#' @description
#' This set of functions revolve around the polygamma functions.
#'
#' @param x,y numeric. The points to evaluate the function.
#' @param p integer. The p-variate Gamma function.
#' @param log logical. Should the logarithm of the result be returned?
#'
#' @describeIn idigamma inverse digamma function.
#'
#' @return numeric. The evaluated function.
#'
#' @export
#'
#' @examples
#' idigamma(2)
#' Ddigamma(2, 3)
#' Dtrigamma(2, 3)
#' gammap(1:3, 3)
idigamma <- function(x) {
  igamma(x)#idigamma_approx(x)
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
