#-------------------------------------------------------------------------------
# Distributions
#-------------------------------------------------------------------------------

# Generics ----

setGeneric("shape", function(object) standardGeneric("shape"))
setGeneric("shape<-", function(object, value) standardGeneric("shape<-"))

# Dirichlet      ----

#' @title The Dirichlet Distribution
#'
#' @param x numeric. The quantile vector.
#' @param shape numeric. The parameter vector.
#' @param log logical. If TRUE, probabilities p are given as log(p).
#' @param n numeric. The number of observations.
#'
#' @return `dDirichlet` returns the evaluated density function.
#'         `rDirichlet` performs Monte-Carlo simulation.
#'
#' @export
#'
#' @examples \dontrun{
#' # Classic R Stats Format
#' dDirichlet(c(0.3, 0.7), shape = c(2, 3))
#' set.seed(1)
#' rDirichlet(10, shape = c(2, 3))
#'
#' # S4 Distribution Class
#' d_dirichlet <- Dirichlet(shape = c(2, 3))
#' d(d_dirichlet)(c(0.3, 0.7))
#' set.seed(1)
#' r(d_dirichlet)(10)
#' }
dDirichlet <- function(x, shape, log = FALSE) {

  if (length(x) != length(shape)) {
    stop("The lengths of x (", length(x), ") and shape (",
         length(shape), ") must be equal.")
  }

  ld <- log(gamma(sum(shape))) - sum(log(gamma(shape))) + sum((shape - 1) * log(x))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname dDirichlet
rDirichlet <- function(n, shape) {

  k <- length(shape)
  x <- matrix(nrow = n, ncol = k)
  for (j in 1:k) {
    x[, j] <- stats::rgamma(n, shape[j], 1)
  }

  t(x / rowSums(x))

}

setClass("DirichletParameter",
 representation = representation(shape = "numeric"),
 prototype = prototype(shape = c(1, 1),
                       name = gettext("Parameter of a Dirichlet distribution")),
 contains = "Parameter"
)

#' @title Dirichlet Distribution S4 Class
#'
#' @slot shape numeric. The parameter vector.
#'
#' @return An object of class `Dirichlet`.
#' @export
#'
#' @inherit dDirichlet examples
Dirichlet <- setClass("Dirichlet",
 slots = list(shape = "numeric"),
 prototype = prototype(
   r = function(n) {
     rDirichlet(n, shape = c(1, 1))
   },
   d = function(x, log = FALSE) {
     dDirichlet(x, shape = c(1, 1), log = log)
   },
   param = new("DirichletParameter"),
   .logExact = TRUE,
   .lowerExact = TRUE
 ),
 contains = "AbscontDistribution"
)

## Access methods
setMethod("shape", "DirichletParameter", function(object) object@shape)

## Replace Methoden
setReplaceMethod("shape", "DirichletParameter",
                 function(object, value){ object@shape <- value; object})

setValidity("DirichletParameter", function(object){
  if (any(shape(object) <= 0)) {
    stop("shape has to be positive")
  } else {
    return(TRUE)
  }
})

## wrapped access methods
setMethod("shape", "Dirichlet",
          function(object) { shape(distr::param(object)) })

## wrapped replace methods
setMethod("shape<-", "Dirichlet",
          function(object, value) { new("Dirichlet", shape = value(object)) })

setMethod("initialize", "Dirichlet",
  function(.Object, shape = c(1, 1)) {
    .Object@img <- new("Reals")
    .Object@param <- new("DirichletParameter", shape = shape)
    .Object@r <- function(n){}
    .Object@d <- function(x, log = FALSE){}
    .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
    .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){}
    body(.Object@r) <- substitute(
      { rDirichlet(n, shape = shapeSub) },
      list(shapeSub = shape)
    )
    body(.Object@d) <- substitute(
      { dDirichlet(x, shape = shapeSub, log = log) },
      list(shapeSub = shape)
    )
    body(.Object@p) <- substitute(
      { pDirichlet(q, shape = shapeSub, lower.tail = lower.tail,
                   log.p = log.p) },
      list(shapeSub = shape)
    )
    body(.Object@q) <- substitute(
      { qDirichlet(p, shape = shapeSub, lower.tail = lower.tail,
                   log.p = log.p) },
      list(shapeSub = shape)
    )
    .Object@.withSim   <- FALSE
    .Object@.withArith <- FALSE
    .Object
})

# Multivariate Gamma ----

#' @title The Multivariate Gamma Distribution
#'
#' @param x numeric. The quantile vector.
#' @param shape numeric. The shape parameter vector.
#' @param scale numeric. The scale parameter vector.
#' @param log logical. If TRUE, probabilities p are given as log(p).
#' @param n numeric. The number of observations.
#'
#' @return `dMGamma` returns the evaluated density function.
#'         `rMGamma` performs Monte-Carlo simulation.
#'
#' @export
#'
#' @examples \dontrun{
#' # Classic R Stats Format
#' dMGamma(c(4, 6), shape = c(2, 3), scale = 2)
#' set.seed(1)
#' rMGamma(10, shape = c(2, 3), scale = 2)
#'
#' # S4 Distribution Class
#' d_MGamma <- MGamma(shape = c(2, 3), scale = 2)
#' d(d_MGamma)(c(4, 6))
#' set.seed(1)
#' r(d_MGamma)(10)
#' }
dMGamma <- function(x, shape, scale, log = FALSE) {

  if (length(x) != length(shape)) {
    stop("The lengths of x (", length(x), ") and shape (",
         length(shape), ") must be equal.")
  }

  z <- fd(x)
  a0 <- sum(shape)

  ld <- - a0 * log(scale) - sum(lgamma(shape)) - x[length(x)] / scale + sum(shape * log(z))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname dMGamma
rMGamma <- function(n, shape, scale) {

  k <- length(shape)
  x <- matrix(nrow = k, ncol = n)
  for (i in 1:k) {
    x[i, ] <- stats::rgamma(n, shape[i], scale = scale)
  }

  apply(x, 2, cumsum)

}

setClass("MGammaParameter",
 representation = representation(shape = "numeric", scale = "numeric"),
 prototype = prototype(shape = c(1, 1), scale = 1,
                       name = gettext("Parameter of a MGamma distribution")),
 contains = "Parameter"
)

#' @title Multivariate Gamma Distribution S4 Class
#'
#' @slot shape numeric. The shape parameter vector.
#' @slot scale numeric. The scale parameter vector.
#'
#' @return An object of class `MGamma`.
#' @export
#'
#' @inherit dMGamma examples
MGamma <- setClass("MGamma",
  slots = list(shape = "numeric", scale = "numeric"),
  prototype = prototype(
    r = function(n) {
      rMGamma(n, shape = c(1, 1), scale = 1)
    },
    d = function(x, log = FALSE) {
      dMGamma(x, shape = c(1, 1), scale = 1, log = log)
    },
    param = new("MGammaParameter"),
    .logExact = TRUE,
    .lowerExact = TRUE
  ),
  contains = "AbscontDistribution"
)

## Access methods
setMethod("shape", "MGammaParameter",
          function(object) object@shape)

setMethod("scale", signature = c(x = "MGammaParameter"),
          function(x, center = TRUE, scale = TRUE) x@scale)

## Replace Methods
setReplaceMethod("shape", "MGammaParameter",
                 function(object, value){ object@shape <- value; object})

setReplaceMethod("scale", signature = c(object = "MGammaParameter"),
                 function(object, value){ object@scale <- value; object})

setValidity("MGammaParameter", function(object){
  if (any(shape(object) <= 0)) {
    stop("shape has to be positive")
  } else if (object@scale <= 0) {
    stop("scale has to be positive")
  } else {
    return(TRUE)
  }
})

## wrapped access methods
setMethod("shape", "MGamma",
          function(object) { shape(distr::param(object)) })

setMethod("scale", signature = c(x = "MGamma"),
          function(x, center = TRUE, scale = TRUE) {
            distr::param(x)@scale
          })

## wrapped replace methods
setMethod("shape<-", "MGamma",
          function(object, value) { new("MGamma", shape = value(object)) })

setMethod("scale<-", "MGamma",
          function(object, value) {
            new("MGamma", shape = shape(object), scale = value)
          })

setMethod("initialize", "MGamma",
  function(.Object, shape = c(1, 1), scale = 1) {
    .Object@img <- new("Reals")
    .Object@param <- new("MGammaParameter", shape = shape, scale = scale)
    .Object@r <- function(n){}
    .Object@d <- function(x, log = FALSE){}
    .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
    .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){}
    body(.Object@r) <- substitute(
      { rMGamma(n, shape = shapeSub, scale = scaleSub) },
      list(shapeSub = shape, scaleSub = scale)
    )
    body(.Object@d) <- substitute(
      { dMGamma(x, shape = shapeSub, scale = scaleSub, log = log) },
      list(shapeSub = shape, scaleSub = scale)
    )
    body(.Object@p) <- substitute(
      { pMGamma(q, shape = shapeSub, scale = scaleSub, lower.tail = lower.tail,
                   log.p = log.p) },
      list(shapeSub = shape, scaleSub = scale)
    )
    body(.Object@q) <- substitute(
      { qMGamma(p, shape = shapeSub, scale = scaleSub, lower.tail = lower.tail,
                   log.p = log.p) },
      list(shapeSub = shape, scaleSub = scale)
    )
    .Object@.withSim   <- FALSE
    .Object@.withArith <- FALSE
    .Object
})
