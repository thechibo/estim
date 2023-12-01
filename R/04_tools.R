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
  params <- methods::slot(D, "param")
  params <- s4_to_list(params)
  params["name"] <- NULL
  params["ncp"] <- NULL
  unlist(params)
}

get_unknown_params <- function(D) {
  prm <- get_params(D)
  if (is(D, "Binom")) {
    return(prm["prob"])
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

  params <- methods::slot(D, "param")
  slot(params, prm$name)[prm$pos] <- prm$val[i]

  x <- c(s4_to_list(params), Class = class(D))
  x["name"] <- NULL

  do.call("new", x)

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Gamma Function         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#' @importFrom distr igamma
#' @export
#'
#' @seealso [distr::igamma()]
#'
#' @examples \dontrun{
#' idigamma(2)
#' Ddigamma(2, 3)
#' Dtrigamma(2, 3)
#' gammap(1:3, 3)
#' }
idigamma <- function(x) {
  distr::igamma(x)
}

#' @describeIn idigamma digamma difference function.
Ddigamma <- function(x, y) {
  digamma(x) - digamma(y)
}

#' @describeIn idigamma trigamma difference function.
Dtrigamma <- function(x, y) {
  trigamma(x) - trigamma(y)
}

#' @describeIn idigamma p-variate gamma function
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

inv2x2 <- function(x) {
  det <- x[1, 1] * x[2, 2] - x[1, 2] * x[2, 1]

  if (det == 0) {
    return("The matrix is singular, its inverse does not exist.")
  } else {
    inv <- matrix(c(x[2, 2], -x[2, 1], -x[1, 2], x[1, 1]), nrow = 2)
    inv <- inv / det
    return(inv)
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multivariate Gamma     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

fd <- function(x) {
  if (is.matrix(x)) {
    return(rbind(x[1, ], diff(x)))
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
