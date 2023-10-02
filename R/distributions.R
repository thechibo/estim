#-------------------------------------------------------------------------------
# Distributions
#-------------------------------------------------------------------------------

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

setGeneric("shape", function(object) standardGeneric("shape"))
setGeneric("shape<-", function(object, value) standardGeneric("shape<-"))

## Access methods
setMethod("shape", "DirichletParameter", function(object) object@shape)

## Replace Methoden
setReplaceMethod("shape", "DirichletParameter",
                 function(object, value){ object@shape <- value; object})

setValidity("DirichletParameter", function(object){
  if(any(shape(object) <= 0))
      stop("shape has to be positive")
  else return(TRUE)
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

# Matrix Gamma   ----

#' @title The Dirichlet Distribution
#'
#' @param X matrix. The random matrix.
#' @param shape numeric. The shape parameter.
#' @param Sigma matrix. The matrix parameter.
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
#' X <- matrix(c(2, 1, 3, 2), 2, 2)
#' dMGamma(X, shape = 3, Sigma = diag(2))
#' set.seed(1)
#' rMGamma(1, shape = 3, Sigma = diag(2))
#'
#' # S4 Distribution Class
#' X <- matrix(c(2, 1, 3, 2), 2, 2)
#' d_mgamma <- MGamma(shape = 3, Sigma = diag(2))
#' d(d_mgamma)(X)
#' set.seed(1)
#' r(d_mgamma)(1)
#' }
dMGamma <- function(X, shape = 1, Sigma = diag(2), log = FALSE) {

  p <- nrow(Sigma)
  ldetS <- log(det(Sigma))
  ldetX <- log(det(X))

  ld <- - shape * ldetS - lgammap(shape, p) + (shape - (p + 1) / 2) * ldetX - sum(diag(solve(Sigma) %*% X))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

#' @rdname dMGamma
#' @importFrom matrixsampling rmatrixgamma
rMGamma <- function(n, shape = 1, Sigma = diag(2)) {

  matrixsampling::rmatrixgamma(n, nu = shape, theta = 1, Sigma = Sigma, checkSymmetry = TRUE)

}

setClass("MGammaParameter",
         representation = representation(shape = "numeric", Sigma = "matrix"),
         prototype = prototype(shape = 2, Sigma = diag(2),
                               name = gettext("Parameter of a MGamma distribution")),
         contains = "Parameter"
)

#' @title Matrix Gamma Distribution S4 Class
#'
#' @slot shape numeric. The shape parameter.
#' @slot Sigma matrix. The matrix parameter.
#'
#' @return An object of class `MGamma`.
#' @export
#'
#' @inherit dMGamma examples
MGamma <- setClass("MGamma",
 slots = list(shape = "numeric",
              Sigma = "matrix"),
 prototype = prototype(
   r = function(n){
     rMGamma(n, shape = 2, Sigma = diag(2))
   },
   d = function(x, log = FALSE){
     dMGamma(x, shape = 2, Sigma = diag(2), log = log)
   },
   param = new("MGammaParameter"),
   .logExact = TRUE,
   .lowerExact = TRUE
 ),
 contains = "AbscontDistribution"
)

setGeneric("Sigma", function(object) standardGeneric("Sigma"))
setGeneric("Sigma<-", function(object, value) standardGeneric("Sigma<-"))

## Access methods
setMethod("shape", "MGammaParameter", function(object) object@shape)
setMethod("Sigma", "MGammaParameter", function(object) object@Sigma)

## Replace Methoden
setReplaceMethod("shape", "MGammaParameter",
                 function(object, value){ object@shape <- value; object})
setReplaceMethod("Sigma", "MGammaParameter",
                 function(object, value){ object@Sigma <- value; object})

setValidity("MGammaParameter", function(object){
  if(shape(object) <= 0) {
    stop("shape has to be positive")
  } else if (!matrixcalc::is.positive.definite(Sigma(object))){
    stop("Sigma has to be a positive definite matrix")
  } else {
    return(TRUE)
  }
})

## wrapped access methods
setMethod("shape", "MGamma", function(object) shape(distr::param(object)))
setMethod("Sigma", "MGamma", function(object) Sigma(distr::param(object)))

## wrapped replace methods
setMethod("shape<-", "MGamma",
          function(object, value) new("MGamma", shape = value(object)))
setMethod("Sigma<-", "MGamma",
          function(object, value) new("MGamma", Sigma = value(object)))

setMethod("initialize", "MGamma",
          function(.Object, shape = 2, Sigma = diag(2)) {
  .Object@img <- new("Reals")
  .Object@param <- new("MGammaParameter", shape = shape, Sigma = Sigma)
  .Object@r <- function(n){}
  .Object@d <- function(x, log = FALSE){}
  .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
  .Object@q <- function(p, lower.tail = TRUE, log.p = FALSE){}
  body(.Object@r) <- substitute(
    { rMGamma(n, shape = shapeSub, Sigma = SigmaSub) },
    list(shapeSub = shape, SigmaSub = Sigma)
  )
  body(.Object@d) <- substitute(
    { dMGamma(x, shape = shapeSub, Sigma = SigmaSub, log = log) },
    list(shapeSub = shape, SigmaSub = Sigma)
  )
  body(.Object@p) <- substitute(
    { pMGamma(q, shape = shapeSub, Sigma = SigmaSub, lower.tail = lower.tail,
                   log.p = log.p) },
    list(shapeSub = shape, SigmaSub = Sigma)
  )
  body(.Object@q) <- substitute(
    { qMGamma(p, shape = shapeSub, Sigma = SigmaSub, lower.tail = lower.tail,
                   log.p = log.p) },
    list(shapeSub = shape, SigmaSub = Sigma)
  )
  .Object@.withSim   <- FALSE
  .Object@.withArith <- FALSE
  .Object
})
