# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribution Calculus                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculus               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Distribution Calculus
#' @name calculus
#'
#' @param x,e1,e2 objects of subclass `Distribution`.
#' @param na.rm logical. Should missing values be removed?
#' @param ... extra arguments.
#'
#' @return All calculations return Distribution objects (specifically, objects
#' of a class that is a subclass of `Distribution`), accordingly to the property
#' at hand.
#'
#' @examples
#' # -----------------------------------------------------
#' # Distribution Calculus Example
#' # -----------------------------------------------------
#'
#' library(estimators)
#'
#' # Normal location - scale transformation
#' x <- Norm(mean = 2, sd = 3)
#' y <- 3 * x + 1 # Norm(mean = 7, sd = 9)
#'
#' # Addition of two independent Normal random variables
#' x1 <- Norm(mean = 1, sd = 3)
#' x2 <- Norm(mean = 2, sd = 4)
#' x3 <- x1 + x2 # Norm(mean = 3, sd = 5)
NULL

#' @rdname calculus
setMethod("+", signature = c(e1 = "Norm", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1@mean + e2@mean,
                 sd = sqrt(e1@sd ^ 2 + e2@sd ^ 2))
          })

#' @rdname calculus
setMethod("+", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 + e2@mean, sd = e2@sd)
          })

#' @rdname calculus
setMethod("+", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            e2 + e1
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "Norm", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1@mean - e2@mean,
                 sd = sqrt(e1@sd ^ 2 + e2@sd ^ 2))
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 - e2@mean, sd = e2@sd)
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            Norm(mean = e1@mean - e2, sd = e2@sd)
          })

#' @rdname calculus
setMethod("*", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 * e2@mean, sd = e1 * e2@sd)
          })

#' @rdname calculus
setMethod("*", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            e2 * e1
          })

#' @rdname calculus
setMethod("/", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            (1 / e2) * e1
          })

#' @rdname calculus
setMethod("sum", signature = c(x = "Norm", na.rm = "logical"),
          function(x, ..., na.rm = FALSE) {
            d <- list(x, ...)
            m <- unlist(lapply(d, FUN = function(x) {x@mean}))
            s <- unlist(lapply(d, FUN = function(x) {x@sd}))
            Norm(mean = sum(m), sd = sqrt(sum(s ^ 2)))
          })

#' @rdname calculus
setMethod("exp", signature = c(x = "Norm"),
          function(x) {
            Lnorm(meanlog = x@mean, sdlog = x@sd)
          })
