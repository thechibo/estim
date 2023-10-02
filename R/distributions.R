#-------------------------------------------------------------------------------
# Distributions
#-------------------------------------------------------------------------------

# Dirichlet      ----

if(!isGeneric("shape"))
  setGeneric("shape", function(object) standardGeneric("shape"))
if(!isGeneric("shape<-"))
  setGeneric("shape<-", function(object, value) standardGeneric("shape<-"))

dDirichlet <- function(x, prm, log = FALSE) {

  ld <- log(gamma(sum(prm))) - sum(log(gamma(prm))) + sum((prm-1) * log(x))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

rDirichlet <- function(n, shape) {

  k <- length(shape)
  x <- matrix(nrow = n, ncol = k)
  for (j in 1:k) {
    x[, j] <- rgamma(n, shape[j], 1)
  }

  t(x / rowSums(x))

}

#' @export
setClass("DirichletParameter",
         representation = representation(shape = "numeric"),
         prototype = prototype(shape = 1,
                               name = gettext("Parameter of a Dirichlet distribution")),
         contains = "Parameter"
)

#' @export
setClass("Dirichlet",
         prototype = prototype(
           r = function(n){
             rDirichlet(n, shape = 1)
           },
           d = function(x, log = FALSE){
             dDirichlet(x, prm = c(1, 1), log = log)
           },
           param = new("DirichletParameter"),
           .logExact = TRUE,
           .lowerExact = TRUE
         ),
         contains = "AbscontDistribution"
)

Dirichlet <- function(shape = 1) {
  new("Dirichlet", shape = shape)
}

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
setMethod("shape", "Dirichlet", function(object) shape(param(object)))

## wrapped replace methods
setMethod("shape<-", "Dirichlet",
          function(object, value) new("Dirichlet", shape = value(object)))

setMethod("initialize", "Dirichlet",
          function(.Object, shape = 1) {
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

dMGamma <- function(X, shape, Sigma, log = FALSE) {

  p <- nrow(Sigma)
  ldetS <- log(det(Sigma))
  ldetX <- log(det(X))

  ld <- - shape * ldetS - lgammap(shape, p) + (shape - (p + 1) / 2) * ldetX - sum(diag(solve(Sigma) %*% X))

  if (!log) {
    ld <- exp(ld)
  }

  ld

}

rMGamma <- function(n, shape, Sigma) {

  matrixsampling::rmatrixgamma(n, nu = shape, theta = 1, Sigma = Sigma, checkSymmetry = TRUE)

}

setGeneric("Sigma", function(object) standardGeneric("Sigma"))
setGeneric("Sigma<-", function(object, value) standardGeneric("Sigma<-"))

#' @export
setClass("MGammaParameter",
         representation = representation(shape = "numeric", Sigma = "matrix"),
         prototype = prototype(shape = 2, Sigma = matrix(c(1, 0, 0, 1), 2, 2),
                               name = gettext("Parameter of a MGamma distribution")),
         contains = "Parameter"
)

#' @export
setClass("MGamma",
         prototype = prototype(
           r = function(n){
             rMGamma(n, shape = 2, Sigma = matrix(c(1, 0, 0, 1), 2, 2))
           },
           d = function(x, log = FALSE){
             dMGamma(x, shape = 2, Sigma = matrix(c(1, 0, 0, 1), 2, 2), log = log)
           },
           param = new("MGammaParameter"),
           .logExact = TRUE,
           .lowerExact = TRUE
         ),
         contains = "AbscontDistribution"
)

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
setMethod("shape", "MGamma", function(object) shape(param(object)))
setMethod("Sigma", "MGamma", function(object) Sigma(param(object)))

## wrapped replace methods
setMethod("shape<-", "MGamma",
          function(object, value) new("MGamma", shape = value(object)))
setMethod("Sigma<-", "MGamma",
          function(object, value) new("MGamma", Sigma = value(object)))

setMethod("initialize", "MGamma",
          function(.Object, shape = 2, Sigma = matrix(c(1, 0, 0, 1), 2, 2)) {
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

MGamma <- function(shape = 2, Sigma = matrix(c(1, 0, 0, 1), 2, 2)) {
  new("MGamma", shape = shape, Sigma = Sigma)
}
